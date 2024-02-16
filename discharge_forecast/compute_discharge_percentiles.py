import pandas as pd
from glob import glob
from functools import reduce
import numpy as np
from datetime import date

from discharge_forecast.config import * 


def compute_percentiles(fc_start,percentiles,catchments):
    pcols = ['p{0:d}'.format(p) for p in percentiles]
    pcols[0] = 'min'; pcols[-1] = 'max'

    fcens_filelist = sorted(glob(proj_base + '/results/discharge_forecast/{1:s}/daily21d_{0:s}/daily21d_{0:s}_ens*'.format(fc_start,catchments)))
    droplist = ['utm_x','utm_y','lon','lat','area','month','year','day']
    df_list = []
    for f in fcens_filelist:
        DF = pd.read_csv(f).drop(columns=droplist)
        df_list.append(DF[pd.to_datetime(DF.date) > pd.Timestamp(fc_start)])
    fcens = reduce(lambda  left,right: pd.merge(left,right,on=['date','stat_id'], how='outer'), df_list)
    
    q_col_names = [col for col in fcens.columns if 'q_sim_mm' in col]
    q_col_cumecs_names = [col for col in fcens.columns if 'q_sim_cumecs' in col]

    fcens.loc[:,pcols] = np.percentile(fcens[q_col_names],percentiles,axis=1).T
    fcens = fcens.drop(columns=q_col_names + q_col_cumecs_names)

    fcens.to_csv(proj_base + '/results/discharge_forecast/{1:s}/daily21d_{0:s}/daily21d_{0:s}_empiricalpercentiles.csv'.format(fc_start,catchments),index=False)


if __name__ == '__main__':
    percentiles = np.arange(0,101,5)
    init_date = date.today().strftime('%Y-%m-%d')
    for catchments in ['smaakraft', 'nve']:
        compute_percentiles(init_date,percentiles,catchments)