import xarray as xr
import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree
from glob import glob
import re
import os

from discharge_forecast.config import *

# location of ecmwf s2s hindcasts:
ecmwf_dir = '/projects/NS9873K/etdu/raw/s2s/mars/ecmwf/hindcast/sfc/6hourly/'

# find all hindcast files for which all 6 files (3 streams, cf and pf) exist:
resol = ['0.25x0.25','vr_0.5x0.5','0.5x0.5'] # all resolution keys
fctype = ['cf','pf'] # all forecast type keys (control and perturbed)

def find_complete_hc_set(vrbl,dir=ecmwf_dir):

    all_tp = glob('{0:s}{1:s}/*'.format(dir,vrbl))
    dts = [re.search('_(\d{4}-\d{2}-\d{2})_',tpf).groups()[0] for tpf in all_tp]
    unique_dts = list(set(dts))
    missing = []
    full = []

    if vrbl == 'tp':
        res = resol
    else:
        res = [resol[0],resol[-1]]

    for udt in unique_dts:
        ud_list = glob('{0:s}{2:s}/*_{1:s}_*'.format(dir,udt,vrbl))
        if len(ud_list) < len(res)*len(fctype):
            missing.append(udt)
        else:
            MV = [re.search('(CY\d{2}R\d{1})',udl).groups()[0] for udl in ud_list]
            # make sure model version is unique for initialization date
            if len(list(set(MV))) == 1:
                file_exist = [os.path.exists('{0:s}{1:s}/{1:s}_{2:s}_{3:s}_{4:s}_{5:s}.nc'.format(dir,vrbl,MV[0],RES,udt,FCTYPE)) for RES in res for FCTYPE in fctype]
                # make sure all required files exist:
                if np.array(file_exist).sum() == len(res)*len(fctype):
                    full.append(udt)
                else:
                    print('sth is weird')
                    missing.append(udt)
            else:
                print('multiple model versions found for same date...')

    return sorted(full),sorted(missing)


def mv_for_date(date,var,dir = ecmwf_dir):
    flist = glob('{0:s}{1:s}/{1:s}_*_{2:s}_{3:s}_{4:s}.nc'.format(dir,var,resol[0],date,fctype[0]))
    if len(flist) == 1:
        return re.search('CY\d{2}R\d{1}',flist[0])[0]

def get_ctchmnt_latlons(ctchmnts):

    csv_file_path = proj_base + '/results/catchment_properties/{0:s}/{0:s}_prop_and_clim.csv'.format(ctchmnts)

    # Read the CSV file into a pandas dataframe:
    data_frame = pd.read_csv(csv_file_path)

    df_lat_lon = data_frame[['stat_id','latitude','longitude']].set_index('stat_id')

    return df_lat_lon


def retrieve_closest_gp(date,var,dir=ecmwf_dir,ctchmnts='smaakraft'):
    """
    """

    mv = mv_for_date(date,var,dir=dir)

    if var == 'tp':
        res = resol
    else:
        res = [resol[0],resol[-1]]

    # get catchment coordinates:
    ctch_lat_lon = get_ctchmnt_latlons(ctchmnts)
    # transform catchment center coords to radians
    xcoords = np.deg2rad(
        np.stack([ctch_lat_lon.latitude,ctch_lat_lon.longitude],axis=-1)
    )

    # load overlap stream:
    ds_lowres =  xr.open_dataarray('{0:s}{1:s}/{1:s}_{2:s}_{3:s}_{4:s}_{5:s}.nc'.format(dir,var,mv,res[1],date,fctype[0])).isel(time=0)
    # stack coordinates:
    ds_stack = ds_lowres.stack(point=['latitude','longitude'])
    # transform to radians
    ycoords = np.deg2rad(
        np.stack([ds_stack.latitude,ds_stack.longitude],axis=-1)
    )        

    tree    = BallTree(ycoords,metric="haversine")
    # find indices:
    k_indices = tree.query(xcoords,return_distance=False).squeeze()

    # now go through all files and use k_indices to retrieve closest grid point to catchments:
    DA = []
    for FCTYPE in fctype:
        DA_1type = []
        for RES in res:
            with xr.open_dataarray('{0:s}{1:s}/{1:s}_{2:s}_{3:s}_{4:s}_{5:s}.nc'.format(dir,var,mv,RES,date,FCTYPE)) as da:
                if RES == '0.25x0.25':
                    ds_lowres_extended = ds_lowres.reindex({'time':da.time})
                    da = da.interp_like(ds_lowres_extended)
                da_stck = da.stack(point=['latitude','longitude'])
            DA_1type.append(da_stck.isel(point=k_indices))
        DA_1type = xr.concat(DA_1type,dim='time') # will have two values for the overlap date
        if FCTYPE == 'cf':
            DA_1type = DA_1type.expand_dims({'number':1}).assign_coords({'number': [0]})
        DA.append(DA_1type)

    DA_out = xr.concat(DA,dim='number').assign_coords(station=('point',ctch_lat_lon.index.to_list()))
    DA_out.attrs['model_version'] = mv
    return DA_out

def precip_6to6(dataarray,time_offset=np.timedelta64(6,'h'),timedim='time'):
    """
    compute daily precip sums from 6AM to 6AM from the 6hrly accumulated precip (since initialization)
    date of daily sums is given as the date of the end point of accumulation (i.e. values are representative of the 24hrs prior to the date tag)
    """
    
    dataarray_diff = dataarray.diff(dim=timedim)
    # find the indices of all unique time values (first occurrence)
    dataarray_unique = dataarray_diff.isel(time=np.unique(dataarray_diff.time,return_index=True)[-1])
    # compute 24H sums using the given offset and label values with the end date of the aggregation period
    return dataarray_unique.resample(time='24H', offset=time_offset,label='right').sum()

def temp_6to6(dataarray):
    """
    compute daily temperature averages (6AM to 6AM) from instantaneous temperature values
    ASSUMES 6HRLY DATA AS INPUT!
    uses each 6AM value twice!
    date of daily averages is given as the date of the end point of averaging (i.e. values are representative of the 24hrs prior to the date tag)
    """
    
    dataarray_res = dataarray.rolling(time=5).mean()
    
    return dataarray_res.sel(time=dataarray_res.time.dt.hour == 6).dropna(dim='time')


if __name__ == '__main__':
    savedir = proj_base + '/catchment_hindcast'
    for VAR in ['tp','t2m']:
        print(VAR)
        full,missing = find_complete_hc_set(VAR)
        for date in full:
            print(date)
            filename = '{0:s}/{1:s}/{1:s}_{2:s}.nc'.format(savedir,VAR,date)
            DA = retrieve_closest_gp(date,VAR)
            if VAR == 'tp':
                DA_6to6 = precip_6to6(DA)
            elif VAR == 't2m':
                DA_6to6 = temp_6to6(DA)

            # 'flatten' the time axis of the data 
            DA_stacked = DA_6to6.stack(valid_time=['hdate','time'])
            vd_timestamp = []
            init = []
            for vd in DA_stacked.valid_time.values:
                init.append(pd.Timestamp(int(str(vd[0])[:4]),int(str(vd[0])[4:6]),int(str(vd[0])[6:])))
                valid_year = int(str(vd[0])[:4]) if (vd[-1].month - int(str(vd[0])[4:6]) >= 0) else int(str(vd[0])[:4]) + 1
                vd_timestamp.append(pd.Timestamp(valid_year,vd[-1].month,vd[-1].day,vd[-1].hour))
            DA_6to6_final = DA_stacked.drop_vars({'hdate', 'valid_time', 'time'}).assign_coords(valid_time=('valid_time',vd_timestamp)).assign_coords(init_time=('valid_time',init))

            DA_6to6_final = DA_6to6_final.reset_index('point').set_index({'point':'station'}).rename({'point':'station'})
            DA_6to6_final.to_netcdf(filename)
