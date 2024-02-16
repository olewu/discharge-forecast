import pandas as pd
from glob import glob

from discharge_forecast.config import *



def merge_senorge_daily(data_dir, ctchmnt_class):
    tmp_path = '{:}/tmp'.format(data_dir)
    outfile_daily = '{1:}/historical_data/senorge/{0:s}/seNorge_daily_{0:s}_1960-1989.csv'.format(ctchmnt_class,data_dir)
    outfile_coords = '{1:}/historical_data/senorge/{0:s}/sN_allcoords_{0:s}_1960-1989.csv'.format(ctchmnt_class,data_dir)
    outfile_lonlat = '{1:}/historical_data/senorge/{0:s}/seNorge_lonlat_{0:s}_1960-1989.csv'.format(ctchmnt_class,data_dir)

    # daily values:
    all_daily = sorted(glob(tmp_path + '/seNorge_daily_{0:s}_*.csv'.format(ctchmnt_class)))
    df_daily = pd.concat(map(pd.read_csv, all_daily), ignore_index=True)
    df_daily.to_csv(outfile_daily,index=False)

    if ctchmnt_class == 'smaakraft':
        # coordinate values:
        all_coords = sorted(glob(tmp_path + '/sN_allcoords_{0:s}_*.csv'.format(ctchmnt_class)))
        # take only one file (the last) since they're all the same but get generated for every year
        df_coords = pd.read_csv(all_coords[-1], index_col=[0])
        df_coords.to_csv(outfile_coords,index=False)

        # latlon values:
        all_latlon = sorted(glob(tmp_path + '/seNorge_lonlat_{0:s}_*.csv'.format(ctchmnt_class)))
        # take only one file (the last) since they're all the same but get generated for every year
        df_latlon = pd.read_csv(all_latlon[-1], index_col=[0])
        df_latlon.to_csv(outfile_lonlat,index=False)


if __name__ == '__main__':
    
    # input data directory; {data_dir}/tmp/ needs to include all the years to be merged
    data_dir = os.path.join(proj_base,'data')

    merge_senorge_daily(data_dir, ctchmnt_class = 'nve')    