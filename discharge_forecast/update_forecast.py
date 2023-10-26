import requests
import numpy as np
from datetime import datetime, timedelta
import pandas as pd
import xarray as xr
import os
import subprocess as sbp
from glob import glob
import sys

from discharge_forecast.config import *
from discharge_forecast import secrets

def process_temperature_forecasts(data):
    
    temp_values = np.array(
        [[datetime.strptime(entry['time'], "%Y-%m-%dT%H:%M:%SZ"),entry['data']['instant']['details']['air_temperature']] for entry in data['properties']['timeseries']]
    )

    temp_values_24h = [0]
    tt_start = temp_values[0,0]
    time_vec = [tt_start + timedelta(days=1)]
    i = 0
    count = [0]
    for tt,val in temp_values:
        # print(tt,val)
        temp_values_24h[i] += val
        count[i] += 1        
            
        if (tt.hour == 6) and (tt != tt_start):
            temp_values_24h.append(val)
            time_vec.append(tt+timedelta(days=1))
            count.append(1)
            i += 1
    if time_vec[0].hour != 6:
        time_vec[0] = datetime(time_vec[0].year,time_vec[0].month,time_vec[0].day,6)
    if time_vec[0] == time_vec[1]:
        time_vec[0] = time_vec[0] - timedelta(days=1)

    st_24 = xr.DataArray.from_series(pd.Series([tt/cc for tt,cc in zip(temp_values_24h,count)],index = pd.Index(time_vec,name='date')))

    mask = xr.DataArray.from_series(pd.Series([0 if cc >4 else 1 for cc in count],index = pd.Index(time_vec,name='date'))) # mask average for days with at less than 5 forecasts
    if count[0] < 25:
        mask[0] = 1

    return st_24,mask


def process_precipitation_forecasts(data):
    
    # precipitation_amount_6h = data['properties']['timeseries'][80]['data']['next_6_hours']['details']['precipitation_amount']
    precip_values = np.array(
        [[datetime.strptime(entry['time'], "%Y-%m-%dT%H:%M:%SZ"),float(entry['data']['next_6_hours']['details']['precipitation_amount'])] for entry in data['properties']['timeseries'] if 'next_6_hours' in entry['data'].keys() and datetime.strptime(entry['time'], "%Y-%m-%dT%H:%M:%SZ").hour%6 == 0]
    )
    
    # Accumulate precip to daily sum (6AM to 6AM)
    precip_values_24h_accum = [0]
    tt_start = precip_values[0,0]
    time_vec = [tt_start + timedelta(days=1)]
    i = 0
    count = [0]
    for tt,val in precip_values:
        if tt.hour != 6 or tt == tt_start:
            precip_values_24h_accum[i] += val
            count[i] += 1
        else:
            precip_values_24h_accum.append(val)
            # if i != 0 or tt != tt_start:
            time_vec.append(tt+timedelta(days=1))
            count.append(1)
            i += 1


    if time_vec[0].hour != 6:
        time_vec[0] = datetime(time_vec[0].year,time_vec[0].month,time_vec[0].day,6)
    if time_vec[0] == time_vec[1]:
        time_vec[0] = time_vec[0] - timedelta(days=1)

    tp_24 = xr.DataArray.from_series(pd.Series(precip_values_24h_accum,index = pd.Index(time_vec,name='date')))

    mask = xr.DataArray.from_series(pd.Series([0 if cc == 4 else 1 for cc in count],index = pd.Index(time_vec,name='date'))) # mask daily sums for days with at less than 4 forecasts
    
    return tp_24,mask


if __name__ == '__main__':

    today = datetime.today()

    for catchments_from in ['smaakraft','nve']:
        
        desired_fc = '{0:s}/data/regular_downloads/metno/{1:s}/metno_{2:s}T06:00:00Z.csv'.format(proj_base,catchments_from,today.strftime('%Y-%m-%d'))

        existing_fc = sorted(glob('{0:s}/data/regular_downloads/metno/{1:s}/metno_{2:s}*.csv'.format(proj_base,catchments_from,today.strftime('%Y-%m-%d'))))

        fail = 'Failed'

        if catchments_from == 'nve':
            catchprop_nve = pd.read_csv(proj_base + '/results/catchment_properties/nve/catchprop_nveapi.csv')

        if os.path.exists(desired_fc):
            print('Forecast already exists. Skipping. {:}'.format(today))
        else:
            try:
                
                # get file that contains all seNorge grid cells that have centroid inside a catchment
                csv_file_path = '{0:s}/results/catchment_properties/{1:s}/sN_allcoords_{1:s}.csv'.format(proj_base,catchments_from)
                # sn_old = pd.read_csv('{0:s}/data_for_ole/smaakraft_catchments/seNorge_daily_smaakraft.csv'.format(path))

                # Read the CSV file into a pandas dataframe:
                data_frame = pd.read_csv(csv_file_path).rename(columns={'catchid':'stat_id'})

                stations_list = [[stat,data_frame['elev'][id]] for (id,stat) in data_frame['stat_id'].drop_duplicates().items()]
                stations = pd.DataFrame(stations_list,columns=['stat_id','elev'])

                st_24_6to6 = []
                tp_24_6to6 = []
                mask_st = []
                mask_tp = []
                stations_ = {}
                elevation = []
                stat_coll = [] # collect all stations that downloads were executed for

                for _,stat in stations.iterrows():

                    # from NVE, onyl take catchments smaller than 70km2
                    # otherwise downloads will take forever...
                    # think about downsampling (instead of taking point each km)
                    # to be able to include larger catchments
                    if catchments_from == 'nve':
                        if catchprop_nve[catchprop_nve.stat_id == stat.stat_id].area_total.any():
                            total_area = catchprop_nve[catchprop_nve.stat_id == stat.stat_id].area_total.to_list()[0]
                            if total_area > 100:
                                print('skipping {0:}, too large ({1:})'.format(stat.stat_id,total_area))
                                continue
                        else:
                            print('skipping {:} due to lack of area information'.format(stat.stat_id))
                            continue
                    
                    LATS = data_frame[data_frame['stat_id'] == stat['stat_id']]['lat'].to_list()
                    LONS = data_frame[data_frame['stat_id'] == stat['stat_id']]['lon'].to_list()

                    temp_24_6to6 = []; st_mask = []
                    precip_24_6to6 = []; tp_mask = []
                    elev = []

                    for cc,(latitude,longitude) in enumerate(zip(LATS,LONS)):
                        
                        # latitude and longitude of the location you want to get the forecast for
                        latitude = np.round(latitude,4)
                        longitude = np.round(longitude,4)

                        # Construct the API endpoint URL
                        url = 'https://api.met.no/weatherapi/locationforecast/2.0/compact?lat={lat}&lon={lon}'.format(lat=latitude,lon=longitude)

                        # Set headers with your API key
                        headers = {
                            'User-Agent': secrets.metno_api_user,
                            # 'Authorization': f'Bearer {api_key}'
                        }

                        # Make the API request
                        response = requests.get(url, headers=headers)

                        # Check if the request was successful
                        if response.status_code in [200, 203]:
                            # Print the forecast data
                            # print(response.json())
                            data = response.json()

                            data_first_step = datetime.strptime(data['properties']['timeseries'][0]['time'],'%Y-%m-%dT%H:%M:%SZ')

                            filename = '{0:s}/data/regular_downloads/metno/{1:s}/metno_{2:s}.csv'.format(
                                proj_base,catchments_from,data['properties']['timeseries'][0]['time']
                            )

                            if cc == 0 and data_first_step.hour > 6:
                                print('Warning: first time step in forecast is {:} (after 6AM)!'.format(data_first_step))

                            # if the file doesn't exist, check if there is one closer to 6AM initialization and stop if that is the case
                            if existing_fc:
                                des_time = datetime.strptime(desired_fc.split('_')[-1].split('.')[0],'%Y-%m-%dT%H:%M:%SZ')
                                tdiff_ex = [abs((datetime.strptime(exfc.split('_')[-1].split('.')[0],'%Y-%m-%dT%H:%M:%SZ')) - des_time) for exfc in existing_fc]
                                tdiff_new = data_first_step - des_time
                                log = [1 if tdx < tdiff_new  else 0 for tdx in tdiff_ex]
                                if 1 in log:
                                    print('Forecast closer to 6AM already exists')
                                    fail = 'Chose not'
                                    sys.exit()
                            # stop execution if the forecast already exists:
                            if os.path.exists(filename):
                                print('Forecast already exists')
                                fail = 'Chose not'
                                sys.exit()

                            temp_24_6to6_gp,msk_t = process_temperature_forecasts(data)
                            precip_24_6to6_gp,msk_p = process_precipitation_forecasts(data)

                            elev.append(data['geometry']['coordinates'][-1])

                            temp_24_6to6.append(temp_24_6to6_gp)
                            st_mask.append(msk_t)
                            precip_24_6to6.append(precip_24_6to6_gp)
                            tp_mask.append(msk_p)

                        else:
                            print('Error: {0:d} at {1:s} for gp {2:d}'.format(response.status_code,stat['stat_id'],cc))

                        
                    st_24_6to6.append(xr.concat(temp_24_6to6,dim='loc').mean(dim='loc'))
                    tp_24_6to6.append(xr.concat(precip_24_6to6,dim='loc').mean(dim='loc'))
                    mask_st.append(xr.concat(st_mask,dim='loc').sum(dim='loc'))
                    mask_tp.append(xr.concat(tp_mask,dim='loc').sum(dim='loc'))
                    stations_[stat.stat_id] = np.array(elev).mean()
                    elevation.append(np.array(elev).mean())
                    stat_coll.append(stat.stat_id)

                tp = xr.DataArray(xr.concat(tp_24_6to6,dim='catchname').assign_coords(catchname=stat_coll),name='prec')
                st = xr.DataArray(xr.concat(st_24_6to6,dim='catchname').assign_coords(catchname=stat_coll),name='st')
                flag_st = xr.where(xr.DataArray(xr.concat(mask_st,dim='catchname').assign_coords(catchname=stat_coll),name='st_flag'),1.,0)
                flag_tp = xr.where(xr.DataArray(xr.concat(mask_tp,dim='catchname').assign_coords(catchname=stat_coll),name='prec_flag'),1,0)

                stations_ = pd.DataFrame({'stat_id':stat_coll,'elev':elevation})

                DF = xr.merge([tp,st,flag_tp,flag_st]).to_dataframe()
                # add elevation info:
                DF = DF.sort_index(level=[1]).reset_index(level=[0]).sort_index().reset_index('catchname')
                DF = DF.merge(stations_.rename(columns={'stat_id':'catchname'}).set_index('catchname'), left_on='catchname', right_on='catchname')
                DF['elev'] = DF['elev'].round(0)
                DF['date'] = DF['date'].dt.strftime('%Y-%m-%d')

                DF.to_csv(filename,index=False)
                if existing_fc:
                    sbp.call(['rm','-r'] + existing_fc)

            except Exception as e:
                print(e)
                print('{0:s} to retrieve latest forecast {1:}'.format(fail,datetime.today()))