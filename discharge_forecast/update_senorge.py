import pandas as pd
import xarray as xr
from datetime import datetime, timedelta

from discharge_forecast.config import *

# load the required grid points:
def get_indices(catchments_from):

    csv_coord = '{0:s}/data/historical_data/senorge/{1:s}/sN_allcoords_{1:s}.csv'.format(proj_base,catchments_from)
    df_coord = pd.read_csv(csv_coord)

    sidx = (df_coord['xcoor'] - 1).to_list() # subtract 1 from index since the table is generated in R, which uses 1 as the first index!
    sidy = (df_coord['ycoor'] - 1).to_list()

    if catchments_from == 'smaakraft':
        ID = df_coord['stat_id'].to_list()
    elif catchments_from == 'nve':
        ID = df_coord['catchid'].to_list()

    return sidx, sidy, ID

def load_senorge(time='latest'):
    """
    'time' needs to be string, either in ['current','latest'] or a time given as %Y%m%d (days and months zero padded)
    """

    if time in ['current','latest']:
        # Current date:
        date = datetime.today().strftime('%Y%m%d')
    else:
        date = time
        
    # Open latest seNorge nc file on thredds:
    seno_td = xr.open_dataset('https://thredds.met.no/thredds/dodsC/senorge/seNorge_2018/Latest/seNorge2018_{:}.nc'.format(date))

    return seno_td

def organize_dataframe(DF,organized_df):

    DF_out = DF.rename(columns={'time':'date'})
    DF_out['date'] = DF_out['date'].dt.strftime('%Y-%m-%d')
    DF_out['elev'] = organized_df[['catchname','elev']].drop_duplicates().set_index('catchname')
    DF_out = DF_out.reset_index(level=['catchname'])
    try:
        DF_out = DF_out[organized_df.reset_index().columns.to_list()].round(2)
    except:
        DF_out = DF_out[[col for col in organized_df.reset_index().columns.to_list() if col not in ['st_min','st_max']]].round(2)    
    DF_out = DF_out.set_index('catchname').reindex(organized_df['catchname'].drop_duplicates()).rename_axis('catchname').reset_index('catchname')

    return DF_out


def subselect_senorge_to_catchments(catchments_from,time='latest',append_latest=False):

    selection_index_x, selection_index_y, stat_id = get_indices(catchments_from)

    senorge = load_senorge(time)

    ind_x = xr.DataArray(selection_index_x,dims=['Y'])

    rr_ctch = senorge.rr.squeeze()[selection_index_y,ind_x]
    tg_ctch = senorge.tg.squeeze()[selection_index_y,ind_x]
    RR_ctch = rr_ctch.assign_coords({'catchname':xr.DataArray(stat_id,dims=['Y'])}).groupby('catchname').mean()
    TG_ctch = tg_ctch.assign_coords({'catchname':xr.DataArray(stat_id,dims=['Y'])}).groupby('catchname').mean()
    try:

        tn_ctch = senorge.tn.squeeze()[selection_index_y,ind_x]
        tx_ctch = senorge.tx.squeeze()[selection_index_y,ind_x]

        TN_ctch = tn_ctch.assign_coords({'catchname':xr.DataArray(stat_id,dims=['Y'])}).groupby('catchname').mean()
        TX_ctch = tx_ctch.assign_coords({'catchname':xr.DataArray(stat_id,dims=['Y'])}).groupby('catchname').mean()

        DS = xr.Dataset(data_vars = {'prec':RR_ctch,'st':TG_ctch,'st_min':TN_ctch,'st_max':TX_ctch})
    except:
        DS = xr.Dataset(data_vars = {'prec':RR_ctch,'st':TG_ctch})

    # first_file = '{0:s}/catchment_data_latest/{1:s}/seNorge_20230101.csv'.format(path,catchments_from)
    # if os.path.exists(first_file):
    #     sn_old = pd.read_csv(first_file)
    # else:
    sn_old = pd.read_csv('{0:s}/data/historical_data/senorge/{1:s}/seNorge_daily_{1:s}.csv'.format(proj_base,catchments_from),index_col=0)

    dataframe = DS.to_dataframe()

    if catchments_from == 'nve':
        sn_old = sn_old.rename(columns={'catchid':'catchname'})

    organized_dataframe = organize_dataframe(dataframe,sn_old)
    
    df_filename = '{0:s}/data/regular_downloads/senorge/{1:s}/seNorge_{2:s}.csv'.format(
        proj_base,catchments_from,senorge.time.dt.strftime('%Y%m%d').values[0]
    )

    print('writing {:}'.format(df_filename))

    organized_dataframe.to_csv(
        df_filename,
        index=False
    )

    if append_latest:
        file_latest = '{0:s}/data/regular_downloads/senorge/{1:s}/seNorge_2022-12-31_today.csv'.format(
            proj_base,catchments_from,
        )
        organized_dataframe.to_csv(
            file_latest,
            mode = 'a',
            header = not os.path.exists(file_latest),
            index=False
        )

if __name__ == '__main__':

    for catchments_from in ['smaakraft','nve']:
        subselect_senorge_to_catchments(catchments_from)
        # update yesterday's senorge file to include min and max temp (valid from 18 to 18, so only included later)
        subselect_senorge_to_catchments(catchments_from,time=(datetime.today()-timedelta(days=1)).strftime('%Y%m%d'))