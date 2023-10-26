import pandas as pd
from glob import glob

from discharge_forecast.config import *

def merge_senorge_daily(data_dir):
    tmp_path = '{:}/tmp'.format(data_dir)
    outfile_daily = '{:}/historical_data/senorge/smaakraft/seNorge_daily_smaakraft.csv'.format(data_dir)
    outfile_coords = '{:}/historical_data/senorge/smaakraft/sN_allcoords_smaakraft.csv'.format(data_dir)
    outfile_lonlat = '{:}/historical_data/senorge/smaakraft/seNorge_lonlat_smaakraft.csv'.format(data_dir)

    # daily values:
    all_daily = sorted(glob(tmp_path + '/seNorge_daily_smaakraft_*.csv'))
    df_daily = pd.concat(map(pd.read_csv, all_daily), ignore_index=True)
    df_daily.to_csv(outfile_daily,index=False)

    # coordinate values:
    all_coords = sorted(glob(tmp_path + '/sN_allcoords_smaakraft_*.csv'))
    # take only one file (the last) since they're all the same but get generated for every year
    df_coords = pd.read_csv(all_coords[-1], index_col=[0])
    df_coords.to_csv(outfile_coords,index=False)

    # latlon values:
    all_latlon = sorted(glob(tmp_path + '/seNorge_lonlat_smaakraft_*.csv'))
    # take only one file (the last) since they're all the same but get generated for every year
    df_latlon = pd.read_csv(all_latlon[-1], index_col=[0])
    df_latlon.to_csv(outfile_lonlat,index=False)


if __name__ == '__main__':
    
    # input data directory; {data_dir}/tmp/ needs to include all the years to be merged
    data_dir = os.path.join(proj_base,'data')

    merge_senorge_daily(data_dir)    