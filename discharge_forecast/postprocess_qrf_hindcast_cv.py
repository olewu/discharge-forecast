from glob import glob
import re
import os

import xarray as xr
import pandas as pd
import numpy as np

from quantile_forest import RandomForestQuantileRegressor
from scipy.stats import skew,kurtosis

from discharge_forecast.config import *

method = 'qrf'
MINLEAF = 20
NTREE = 500
ret_quantiles = np.arange(.05,1,.05)

dir_seno = os.path.join(proj_base,'data/historical_data/senorge/smaakraft/seNorge_daily_smaakraft.csv')
dir_hc = os.path.join(proj_base,'results/catchment_hindcast/')
target_dir = os.path.join(proj_base,'results/catchment_hindcast/postprocessed/{2:s}_cv_minleaf{1:d}_ntree{0:d}/'.format(MINLEAF,NTREE,method))

os.makedirs(target_dir + 't2m/',exist_ok=True)
os.makedirs(target_dir + 'tp/',exist_ok=True)

# load senorge data in catchment
seno_catch_df = pd.read_csv(dir_seno)

seno_catch_df.date = pd.to_datetime(seno_catch_df.date)   

# convert to xarray.dataset
seno_catch = seno_catch_df.set_index(['date','catchname']).to_xarray()
seno_catch = seno_catch.rename({'date':'valid_time','catchname':'station'})
seno_catch = seno_catch.assign_coords(valid_time = seno_catch.valid_time + pd.Timedelta(hours=6))

lead_times = np.arange(2,21)

for LT in lead_times:

    print(LT)

    t2m_path = dir_hc + 't2m/t2m_{0:0>2d}d.nc'.format(LT)
    tp_path = dir_hc + 'tp/tp_{0:0>2d}d.nc'.format(LT)

    da_t2m = xr.open_dataarray(t2m_path)
    da_tp = xr.open_dataarray(tp_path)
    xda_lt = xr.Dataset(
        dict(t2m = da_t2m,
        tp = da_tp)
    )
    da_t2m.close()
    da_tp.close()

    # split hindcasts
    yrs = sorted(list(set(xda_lt.valid_time.dt.year.values)))

    t2m_pred = []
    tp_pred = []


    # Post-process catchment by catchment:
    for ctchname in xda_lt.station.values:

        print(ctchname)

        hc_lt_catchsel = xda_lt.sel(station=ctchname)
        # get groups to post-process by season:
        obs_sel = seno_catch.sel(station=ctchname)
        
        t2m_pp = []
        tp_pp = []
        for withhold_year in yrs:
            # Prepare predictors for quantile regression forest:
        
            mon = hc_lt_catchsel.valid_time.dt.month.values
            t2m = hc_lt_catchsel.t2m.values.T
            tp = hc_lt_catchsel.tp.values.T

            skew_t2m = skew(t2m,axis=-1)
            skew_t2m[np.isnan(skew_t2m)] = 0
            kurtosis_t2m = kurtosis(t2m,axis=-1)
            kurtosis_t2m[np.isnan(kurtosis_t2m)] = 0

            skew_tp = skew(tp,axis=-1)
            skew_tp[np.isnan(skew_tp)] = 0
            kurtosis_tp = kurtosis(tp,axis=-1)
            kurtosis_tp[np.isnan(kurtosis_tp)] = 0


            t2m_predictors = np.concatenate([
                mon[:,np.newaxis],
                t2m[:,0][:,np.newaxis], # control member
                t2m.mean(-1)[:,np.newaxis],
                t2m.std(-1)[:,np.newaxis],
                skew_t2m[:,np.newaxis],
                kurtosis_t2m[:,np.newaxis],
                np.percentile(t2m,[10,50,90],axis=-1).T,
                np.sum(t2m < -15,axis=-1)[:,np.newaxis]/t2m.shape[-1],
                np.sum(t2m < -5,axis=-1)[:,np.newaxis]/t2m.shape[-1],
                np.sum(t2m > 0,axis=-1)[:,np.newaxis]/t2m.shape[-1],
                np.sum(t2m > 15,axis=-1)[:,np.newaxis]/t2m.shape[-1],
                np.sum(t2m > 5,axis=-1)[:,np.newaxis]/t2m.shape[-1],
            ],
            axis=-1)


            tp_predictors = np.concatenate([
                mon[:,np.newaxis],
                tp[:,0][:,np.newaxis], # control member
                tp.mean(-1)[:,np.newaxis],
                tp.std(-1)[:,np.newaxis],
                skew_tp[:,np.newaxis],
                kurtosis_tp[:,np.newaxis],
                np.percentile(tp,[30,50,90],axis=-1).T,
                np.sum(tp > .3,axis=-1)[:,np.newaxis]/tp.shape[-1],
                np.sum(tp > 1,axis=-1)[:,np.newaxis]/tp.shape[-1],
                np.sum(tp > 3,axis=-1)[:,np.newaxis]/tp.shape[-1],
                np.sum(tp > 5,axis=-1)[:,np.newaxis]/tp.shape[-1],
                np.sum(tp < .05,axis=-1)[:,np.newaxis],
            ],
            axis=-1)

            # split into training and test set

            t2m_train = xr.DataArray(t2m_predictors,dims=('valid_time','predix'),coords=dict(valid_time=hc_lt_catchsel.valid_time)).sel(valid_time=hc_lt_catchsel.valid_time.dt.year!=withhold_year)
            t2m_test = xr.DataArray(t2m_predictors,dims=('valid_time','predix'),coords=dict(valid_time=hc_lt_catchsel.valid_time)).sel(valid_time=hc_lt_catchsel.valid_time.dt.year==withhold_year)

            tp_train = xr.DataArray(tp_predictors,dims=('valid_time','predix'),coords=dict(valid_time=hc_lt_catchsel.valid_time)).sel(valid_time=hc_lt_catchsel.valid_time.dt.year!=withhold_year)
            tp_test = xr.DataArray(tp_predictors,dims=('valid_time','predix'),coords=dict(valid_time=hc_lt_catchsel.valid_time)).sel(valid_time=hc_lt_catchsel.valid_time.dt.year==withhold_year)


            # Get matching observation set from seno:
            obs_train = obs_sel.sel(valid_time = t2m_train.valid_time)
            obs_test = obs_sel.sel(valid_time = t2m_test.valid_time)

            qrf_t2m = RandomForestQuantileRegressor(n_estimators=NTREE,min_samples_leaf=MINLEAF,n_jobs=2)
            qrf_t2m.fit(t2m_train.values,obs_train.st.values)
            
            qrf_tp = RandomForestQuantileRegressor(n_estimators=NTREE,min_samples_leaf=MINLEAF,n_jobs=2)
            qrf_tp.fit(tp_train.values,obs_train.prec.values)

            prediction_t2m = xr.DataArray(
                qrf_t2m.predict(t2m_test.values,quantiles=list(ret_quantiles)),
                dims=('valid_time','quantile'),
                coords=dict(
                    valid_time = t2m_test.valid_time,
                    quantile = ret_quantiles,
                )
            )

            prediction_tp = xr.DataArray(
                qrf_tp.predict(tp_test.values,quantiles=list(np.arange(.05,1,.05))),
                dims=('valid_time','quantile'),
                coords=dict(
                    valid_time = tp_test.valid_time,
                    quantile = ret_quantiles,
                )
            )

            t2m_pp.append(prediction_t2m)
            tp_pp.append(prediction_tp)

        t2m_pred.append(xr.concat(t2m_pp,dim='valid_time'))
        tp_pred.append(xr.concat(tp_pp,dim='valid_time'))

    t2m_pred = xr.concat(t2m_pred,dim='station').assign_coords({'station':xda_lt.station})
    tp_pred = xr.concat(tp_pred,dim='station').assign_coords({'station':xda_lt.station})

    t2m_pred.to_netcdf(target_dir + 't2m/t2m_{:0>2d}d.nc'.format(LT))
    tp_pred.to_netcdf(target_dir + 'tp/tp_{:0>2d}d.nc'.format(LT))
