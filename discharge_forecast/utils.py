import pandas as pd
import os

from discharge_forecast.config import proj_base

def concatenate_sn_fc(fc_init_date,ts_start,catchments_from='smaakraft',db_path=proj_base,save=False):
    """
    concatenate station forecasts obtained from met.no and senorge data
    from a specified start point in the past (relative to forecast init)
    up to forecast initialization.
    fc_init_date,ts_start should be pandas.Timestamp objects!
    """

    assert type(fc_init_date) == pd.Timestamp, 'pass fc_init_date as pandas.Timestamp'

    assert type(ts_start) == pd.Timestamp, 'pass ts_start as pandas.Timestamp'

    assert fc_init_date >= pd.Timestamp(2023,9,1), 'no forecasts archived before Sept 1, 2023'

    start_time_str = '{:0>4d}-{:0>2d}-{:0>2d}'.format(ts_start.year,ts_start.month,ts_start.day)
    
    # forecast part:
    init_time_str = '{:0>4d}-{:0>2d}-{:0>2d}T06:00:00Z'.format(fc_init_date.year,fc_init_date.month,fc_init_date.day)
    dffc = pd.read_csv(db_path + '/data/regular_downloads/metno/{1:s}/metno_{0:s}.csv'.format(init_time_str,catchments_from))

    # updated SN part:
    if ts_start < pd.Timestamp(2022,12,31):
        drange_sn = pd.date_range(pd.Timestamp(2022,12,31),fc_init_date)
    else:
        drange_sn = pd.date_range(ts_start,fc_init_date)

    sn_filelist = ['{0:s}/data/regular_downloads/senorge/{2:s}/seNorge_{1:s}.csv'.format(db_path,drsn.strftime('%Y%m%d'),catchments_from) for drsn in drange_sn]

    # interrupt if the latest seNorge file does not exist yet!
    if not os.path.exists(sn_filelist[-1]):
        print('Latest seNorge file does not exist yet. Interrupting merge operation.')
        return

    df = pd.concat((pd.read_csv(f) for f in sn_filelist), ignore_index=True)

    # more SN if necessary:
    if ts_start < pd.Timestamp(2022,12,31):
        old_sn = pd.read_csv(db_path + '/data/historical_data/senorge/{0:s}/seNorge_daily_{0:s}.csv'.format(catchments_from))
        sn_arch_req = old_sn[pd.to_datetime(old_sn.date) >= ts_start].rename(columns={'catchid':'catchname'})

    hydro_inp = pd.concat(
        [
            sn_arch_req.drop(columns=['st_min','st_max']),
            df.drop(columns=['st_min','st_max']),
            dffc.drop(columns=['prec_flag','st_flag'])
        ]
    ).sort_values(['catchname','date']).reset_index(drop=True)

    if save:
        svpth = db_path + '/results/forecast_input/{2:s}/fc_init_{0:s}_merge_sn_{1:s}.csv'.format(init_time_str,start_time_str,catchments_from)
        hydro_inp.to_csv(svpth,index=False)
        print('saving to {0:s}'.format(svpth))

    return hydro_inp

