from discharge_forecast.utils import concatenate_sn_fc
import pandas as pd

def update_merge(catchments_from,fc_init_date = None,ts_start = None, save = False):
    if fc_init_date is None:
        fc_init_date = pd.Timestamp('today')
    if ts_start is None:
        ts_start = pd.Timestamp(fc_init_date.year-3,1,1)

    merge = concatenate_sn_fc(fc_init_date,ts_start,save=save,catchments_from=catchments_from)

    return merge


if __name__ == '__main__':
    for catchments_from in ['smaakraft','nve']:
        update_merge(catchments_from,save=True)
