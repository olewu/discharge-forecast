from glob import glob
import re
import os

import xarray as xr
import pandas as pd
import numpy as np

from discharge_forecast.post_processing import Qmap

from discharge_forecast.config import *

saveraw = False

dir_seno = os.path.join(proj_base,'data/historical_data/senorge/smaakraft/seNorge_daily_smaakraft.csv')
dir_hc = os.path.join(proj_base,'results/catchment_hindcast/')
target_dir = os.path.join(proj_base,'results/catchment_hindcast/postprocessed/qmap_cv_4seas/')

# load senorge data in catchment
seno_catch_df = pd.read_csv(dir_seno)

seno_catch_df.date = pd.to_datetime(seno_catch_df.date)   

# convert to xarray.dataset
seno_catch = seno_catch_df.set_index(['date','catchname']).to_xarray()
seno_catch = seno_catch.rename({'date':'valid_time','catchname':'station'})
seno_catch = seno_catch.assign_coords(valid_time = seno_catch.valid_time + pd.Timedelta(hours=6))

lead_times = np.arange(1,46)

# collect hindcast data and load one leadtime:
t2m_hc_files = sorted(glob(dir_hc+'t2m/t2m*_????-??-??.nc'))
for LT in lead_times:

    t2m_path = dir_hc + 't2m/t2m_{0:0>2d}d.nc'.format(LT)
    tp_path = dir_hc + 'tp/tp_{0:0>2d}d.nc'.format(LT)

    if os.path.exists(t2m_path) and os.path.exists(tp_path):
        da_t2m = xr.open_dataarray(t2m_path)
        da_tp = xr.open_dataarray(tp_path)
        xda_lt = xr.Dataset(
            dict(t2m = da_t2m,
            tp = da_tp)
        )
        da_t2m.close()
        da_tp.close()
    else:
        xda_t2m_lt = []
        xda_tp_lt = []
        for fi in t2m_hc_files:

            init_date_str = re.search('_(\d{4}-\d{2}-\d{2}).nc',fi).groups()[0]
            lt_date = pd.Timestamp(init_date_str) + pd.Timedelta(days=LT)

            with xr.open_dataarray(fi) as xda_t2m:
                xda_t2m_lt.append(xda_t2m.sel(valid_time=np.logical_and(xda_t2m.valid_time.dt.month == lt_date.month, xda_t2m.valid_time.dt.day == lt_date.day)))

            # corresponding precip hindcasts:
            with xr.open_dataarray(fi.replace('t2m','tp')) as xda_tp:
                xda_tp_lt.append(xda_tp.sel(valid_time=np.logical_and(xda_tp.valid_time.dt.month == lt_date.month, xda_tp.valid_time.dt.day == lt_date.day)))

        xda_lt = xr.merge([xr.concat(xda_tp_lt,dim='valid_time'),xr.concat(xda_t2m_lt,dim='valid_time')])

        xda_lt.t2m = xda_lt.t2m - 273.15
        xda_lt.tp = xda_lt.tp * 1000

        # save the raw forecasts split by lead time for evaluation:
        if saveraw:
            basepath = '/'.join(fi.split('/')[:-2])
            t2m_path = basepath + '/t2m/t2m_{:0>2d}d.nc'.format(LT)
            tp_path = basepath + '/tp/tp_{:0>2d}d.nc'.format(LT)
            xda_lt.t2m.to_netcdf(t2m_path)
            xda_lt.tp.to_netcdf(tp_path)

    # split hindcasts
    yrs = sorted(list(set(xda_lt.valid_time.dt.year.values)))

    ssns = ['MAM','JJA','SON','DJF']

    t2m_pred = []
    tp_pred = []


    # Post-process catchment by catchment:
    for ctchname in xda_lt.station.values:
        hc_lt_catchsel = xda_lt.sel(station=ctchname)
        # get groups to post-process by season:
        hc_seas_groups = hc_lt_catchsel.groupby(hc_lt_catchsel.valid_time.dt.season)
        obs_sel = seno_catch.sel(station=ctchname)
        
        t2m_pp = []
        tp_pp = []
        
        for seas in ssns:
            hc_subs_seas = hc_lt_catchsel.isel(valid_time=hc_seas_groups.groups[seas])
            for withhold_year in yrs:
                # Split into training and testing subset:
                hc_train_subs = hc_subs_seas.sel(valid_time=hc_subs_seas.valid_time.dt.year!=withhold_year)
                hc_test_subs = hc_subs_seas.sel(valid_time=hc_subs_seas.valid_time.dt.year==withhold_year)

                # Get matching observation set from seno:
                obs_train = obs_sel.sel(valid_time = hc_train_subs.valid_time)
                obs_test = obs_sel.sel(valid_time = hc_test_subs.valid_time)

                # Perform quantile mapping:
                
                # Temperature
                qm_t2m = Qmap(obs_train.st.values,hc_train_subs.t2m.values)
                qm_t2m.fit(fill_val=obs_train.st.values.min()-.5)
                prediction_t2m = xr.DataArray(
                    qm_t2m.predict(hc_test_subs.t2m.values),
                    dims=('number','valid_time'),
                    coords=hc_test_subs.coords
                )
                # print('\n\tTemperature {0:s} in {1:s}\n'.format(ctchname,seas))
                # qm_t2m.evaluate_prediction(obs_test.st)

                # Precipitation
                qm_tp = Qmap(obs_train.prec.values,hc_train_subs.tp.values)
                qm_tp.fit(fill_val=0)
                prediction_tp = xr.DataArray(
                    qm_tp.predict(hc_test_subs.tp.values),
                    dims=('number','valid_time'),
                    coords=hc_test_subs.coords
                )
                # print('\n\tPrecipitation {0:s} in {1:s}\n'.format(ctchname,seas))
                # qm_tp.evaluate_prediction(obs_test.prec)

                t2m_pp.append(prediction_t2m)
                tp_pp.append(prediction_tp)

        t2m_pred.append(xr.concat(t2m_pp,dim='valid_time'))
        tp_pred.append(xr.concat(tp_pp,dim='valid_time'))

    t2m_pred = xr.concat(t2m_pred,dim='station')
    tp_pred = xr.concat(tp_pred,dim='station')

    t2m_pred.to_netcdf(target_dir + 't2m/t2m_{:0>2d}d.nc'.format(LT))
    tp_pred.to_netcdf(target_dir + 'tp/tp_{:0>2d}d.nc'.format(LT))
