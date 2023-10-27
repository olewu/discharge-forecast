import requests
import pandas as pd
import os
import time

from discharge_forecast.config import *
from discharge_forecast.hidden_secrets import sildre_api_key

def discharge_dl(catch_id_list,
                 date = pd.Timestamp('today'),
                 base_url='https://hydapi.nve.no/api/v1/Observations?StationID={station}&Parameter={parameter}&ResolutionTime={resolution_time}&ReferenceTime={reference_time}',
                 slow_down=True,
                 ):

    curr_date_str = date.strftime('%Y-%m-%d')
    # define time period that should be downloaded (yesterday 6AM to today 6AM, inclusive)
    ref_time = '{0:s}/{1:s}'.format(
        (date-pd.Timedelta(days=1)).strftime('%Y-%m-%dT06:00'),
        date.strftime('%Y-%m-%dT06:00')
    )

    # 
    out_file = os.path.join(proj_base,'data/regular_downloads/sildre_nve/nve/nve_{date}.csv'.format(date=curr_date_str))

    failed_catch = []

    for catchment_id in catch_id_list:

        url = base_url.format(
            station = catchment_id,
            parameter = 1001,
            resolution_time = 60,
            reference_time = ref_time, #'2023-10-10T06:00/2023-10-11T06:00',
        )

        headers = {
            "Accept": "application/json",
            "X-API-Key": sildre_api_key, 
        }

        response = requests.get(url, headers=headers)

        if response.status_code == 200:
            data = response.json()  # Process the response data as needed
        elif response.status_code == 429:
            # catch station IDs for requests that were denied (usually 1-5 per daily request)
            # will most likely be available if request is resend
            failed_catch.append(catchment_id)
        else:
            print(f"Request for station {catchment_id} failed with status code {response.status_code}")

        # only compute if data exist for station:
        if data['data'][0]['observations']:
            # average hourly data from 6AM to 6AM:
            data_df = pd.DataFrame(data['data'][0]['observations'])
            data_df.time = pd.to_datetime(data_df.time)

            mean_discharge = pd.DataFrame(
                {
                    'stationName' : data['data'][0]['stationName'],
                    'catchname' : catchment_id,
                    data['data'][0]['parameterNameEng'] : data_df.value.mean(),
                    'quality' : data_df.quality.mean(),
                    'date' : data_df.time[len(data_df)-1].strftime('%Y-%m-%d')
                },
                index=[data_df.time[len(data_df)-1]]
            )
        else:
            print('no data for {0:s} ({1:s})'.format(data['data'][0]['stationName'],catchment_id))
            mean_discharge = pd.DataFrame(
                {
                    'stationName' : data['data'][0]['stationName'],
                    'catchname' : catchment_id,
                    data['data'][0]['parameterNameEng'] : 'NaN',
                    'quality' : 0,
                    'date' : data_df.time[len(data_df)-1].strftime('%Y-%m-%d')
                },
                index=[data_df.time[len(data_df)-1]]
            )

        mean_discharge.to_csv(out_file,mode='a',header = not os.path.exists(out_file),index=False)

        # pause the loop for a split second to avoid sending too many requests per second (might result in being black-listed?)
        # this slows down the download but it still finishes in ~60s for the 230 nve catchments
        # with slow_down = False takes ~45s
        if slow_down:
            time.sleep(.05)

    return failed_catch, out_file

if __name__ == '__main__':

    # run the above function for latest observation (daily mean yesterday 6AM to today 6AM)
    nve_prop = pd.read_csv(os.path.join(proj_base,'results/catchment_properties/nve/nve_prop_and_clim.csv'))
    failed_catch, out_file = discharge_dl(nve_prop.stat_id)

    # try again for the catchments where data request was denied (error 429):
    start_time = pd.Timestamp('now')
    while failed_catch:
        failed_catch,_ = discharge_dl(failed_catch)
        
        # time out after 5 minutes:
        time_delta = pd.Timestamp('now') - start_time
        if time_delta.total_seconds() >= 5*60:
            print('request timed out, stations {:s} are still missing'.format(', '.join(failed_catch)))
            break