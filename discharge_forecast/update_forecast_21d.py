import requests
import numpy as np
from datetime import datetime, timedelta, date
import pandas as pd
import xarray as xr
import os
import sys

from discharge_forecast.config import *

def json_extract(data,key,data_type='next_24_hours'):
    values = np.array(
        [[datetime.strptime(entry['time'], "%Y-%m-%dT%H:%M:%SZ"),entry['data'][data_type]['details'][key]] for entry in data['properties']['timeseries']]
    )
    return values


def collect_daily_percentiles(data,variable,shortname=None,percentiles=[10,50,90]):

    extr = [json_extract(data,variable) if pctl == 50 else json_extract(data,'{0:s}_percentile_{1:d}'.format(variable,pctl)) for pctl in percentiles]

    data_coll = pd.DataFrame(
        {
            '{0:s}_{1:d}'.format(shortname,key): np.array(vals[:,1],dtype=float) for key, vals in zip(percentiles, extr)
        },index=extr[0][:,0]
    )

    return data_coll

if __name__ == '__main__':

    today = date.today()

    for catchments_from in ['smaakraft','nve']:
        
        desired_fc = '{0:s}/data/regular_downloads/metno_21d/{1:s}/metno_{2:}.csv'.format(proj_base,catchments_from,today)

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

                temp = []; prec = []
                stations_ = {}
                elevation = []
                stat_coll = [] # collect all stations that downloads were executed for

                for _,stat in stations.iterrows():

                    # from NVE, only take catchments smaller than 70km2
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

                    temp_ = []; precip_ = []; elev = []

                    for cc,(latitude,longitude) in enumerate(zip(LATS,LONS)):
                        
                        # latitude and longitude of the location you want to get the forecast for
                        latitude = np.round(latitude,4)
                        longitude = np.round(longitude,4)

                        # Construct the API endpoint URL (staging version at https://staging.forti.met.no/api/subseasonal/v2/complete?)
                        url = 'https://api.met.no/weatherapi/subseasonal/1.0/complete?lat={lat}&lon={lon}'.format(lat=latitude,lon=longitude)

                        # Set headers with your API key
                        headers = {
                            'User-Agent': 'owul@norceresearch.no',
                            # 'Authorization': f'Bearer {api_key}'
                        }

                        # Make the API request
                        response = requests.get(url, headers=headers)

                        # Check if the request was successful
                        if response.status_code in [200, 203]:
                            # Print the forecast data
                            # print(response.json())
                            data = response.json()

                            # get first time step (should always be 'tomorrow'):
                            data_first_step = datetime.strptime(data['properties']['timeseries'][0]['time'],'%Y-%m-%dT%H:%M:%SZ')
                            # make sure that the forecast is the latest:
                            if not (data_first_step.date() == today + timedelta(days=1)):
                                print('Forecast is old. First date is not tomorrow')
                                sys.exit()

                            filename = '{0:s}/data/regular_downloads/metno_21d/{1:s}/metno_{2:}.csv'.format(
                                proj_base,catchments_from,today
                            )

                            if os.path.exists(filename):
                                print('Forecast already exists')
                                fail = 'Chose not'
                                sys.exit()

                            st_df = collect_daily_percentiles(data,'air_temperature_mean',shortname='st')
                            tp_df = collect_daily_percentiles(data,'precipitation_amount',shortname='prec')

                            st_ds = xr.Dataset.from_dataframe(st_df)
                            prec_ds = xr.Dataset.from_dataframe(tp_df)

                            temp_.append(st_ds)
                            precip_.append(prec_ds)

                            elev.append(data['geometry']['coordinates'][-1])

                        else:
                            print('Error: {0:d} at {1:s} for gp {2:d}'.format(response.status_code,stat['stat_id'],cc))

                        
                    temp.append(xr.concat(temp_,dim=pd.Index(np.arange(len(temp_)),name='loc')).mean(dim='loc'))
                    prec.append(xr.concat(precip_,dim=pd.Index(np.arange(len(temp_)),name='loc')).mean(dim='loc'))
                    stations_[stat.stat_id] = np.array(elev).mean()
                    elevation.append(np.array(elev).mean())
                    stat_coll.append(stat.stat_id)

                stations_ = pd.DataFrame({'stat_id':stat_coll,'elev':elevation})
                
                st = xr.concat(temp,dim=pd.Series(stations_['stat_id'].to_list(),name='catchname'))
                tp = xr.concat(prec,dim=pd.Series(stations_['stat_id'].to_list(),name='catchname'))
                
                DF = xr.merge([tp,st]).to_dataframe()
                DF = DF.sort_index(level=[1]).reset_index(level=[0]).sort_index().reset_index('catchname').rename(columns={'index':'date'})
                DF = DF.merge(stations_.rename(columns={'stat_id':'catchname'}).set_index('catchname'), left_on='catchname', right_on='catchname')
                DF['elev'] = DF['elev'].round(0)
                # shift data plus a day to match the time stamps in seNorge where the date refers to the end of the 24h window!
                DF.date = DF.date + pd.Timedelta(days=1)
                DF['date'] = DF['date'].dt.strftime('%Y-%m-%d')

                DF.to_csv(filename,index=False)
                
            except Exception as e:
                print(e)
                print('{0:s} to retrieve latest forecast {1:}'.format(fail,datetime.today()))