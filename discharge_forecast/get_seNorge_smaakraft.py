import numpy as np
import cartopy.crs as ccrs
import pandas as pd
import xarray as xr
from shapely.geometry import Point
import geopandas as gpd
from netCDF4 import Dataset, num2date
from datetime import datetime

from glob import glob
import re
import os
import sys

from discharge_forecast.config import *

data_path = os.path.join(proj_base,'historical_data/senorge/smaakraft')

# can run this from the dir where script is saved with a year as input
if len(sys.argv) > 1:
    YEAR = [sys.argv[1]]
else: # if no year is passed, go through all years from 1990 to 2023
    YEAR = range(1990,2023)

# files to write (same as in Thea's get_seNorge_smaakraft.R):
# the code appends data to each of these files, so if they already exist, it might create duplicates for already existing catchments (move old files first for a fresh start) 
# dataframe with daily catchment averages from seNorge:
daily_out = os.path.join(data_path,'seNorge_daily_smaakraft.csv')
# dataframe of longitudes and latitudes of the catchment centers
lola_out = os.path.join(data_path,'seNorge_lonlat_smaakraft.csv')
# dataframe with all information on each grid point entering the catchment averages: 
allcord_path = os.path.join(data_path,'sN_allcoords_smaakraft.csv')

# find all catchment shape files:
all_catchments_dir = os.path.join(proj_base,'data/catchment_shape_files/*')
sk_all_catchments_dir = sorted(glob(all_catchments_dir))

# open seNorge geoinfo file (constant over time; contains elevation, easting, northing):
sn_geo = xr.open_dataset('https://thredds.met.no/thredds/dodsC/senorge/geoinfo/seNorge2_dem_UTM33.nc')

# iterate through desired years in seNorge (can go back to 1957)
for year in YEAR:
    # open seNorge
    seno_files_thredds = 'https://thredds.met.no/thredds/dodsC/senorge/seNorge_2018/Archive/seNorge2018_{0:d}.nc'.format(year)
    NC = Dataset(seno_files_thredds)

    # read out coordinate variables:
    longitude = NC.variables['longitude'][:]
    latitude = NC.variables['latitude'][:]
    X = NC.variables['X'][:]
    Y = NC.variables['Y'][:]
    TIME = NC.variables['time']
    # transform time units to string
    time_dt = [datetime(tobj.year,tobj.month,tobj.day) for tobj in num2date(TIME[:],TIME.units,calendar=TIME.calendar)]
    time_pd = pd.to_datetime(time_dt)

    # loop through all catchments:
    for sk_catch_dir in sk_all_catchments_dir:


        # extract catchment name and replace Roman numerals with Arabic
        catchn = sk_catch_dir.split('/')[-1]
        repl = {' I': ' 1', ' II': ' 2', ' III': ' 3'}
        catchn = re.sub(
            '|'.join(repl.keys()),
            lambda match: repl[match.string[match.start():match.end()]],
            catchn
        )

        # change to lowercase and substitute special Norwegian characters and white spaces
        replaces = {'ø': 'oe', 'å': 'aa', 'æ': 'ae', ' ': '_'}
        ctchmnt = re.sub(
            '|'.join(replaces.keys()),
            lambda match: replaces[match.string[match.start():match.end()]],
            catchn.lower()
        )

        # skip catchment if data already exists in the specified daily data csv
        if os.path.exists(daily_out):
            if ctchmnt in pd.read_csv(daily_out).catchname.to_list():
                print('skipping {0:s}, daily senorge data already exists in file'.format(ctchmnt))
                continue

        print(catchn,sk_catch_dir)

        # open catchment shape file as geopandas dataframe (gpdf):
        ctchm_file = '{0:s}/NedbfeltF_v3.shp'.format(sk_catch_dir)
        gdf = gpd.read_file(ctchm_file)

        # initialize projection
        proj = ccrs.PlateCarree()
        # Reproject the GeoDataFrame to the desired projection
        desired_crs = proj.proj4_init
        if gdf.crs is not None and gdf.crs != desired_crs:
            gdf = gdf.to_crs(desired_crs)


        # retrieve geometry from gpdf
        polygon_geometry = gdf['geometry'][0]

        # Get the bounding box of the polygon to cut down the senorge grid
        lat_l = [0 for i in range(2)]
        lon_l = lat_l.copy()
        [lon_l[0],lat_l[0],lon_l[1],lat_l[1]] = polygon_geometry.bounds
        lims = np.array([(X,Y) for X in lon_l for Y in lat_l])

        # cut the bounding box out of senorge:
        log_lat = np.logical_and(lat_l[0] < latitude, latitude < lat_l[-1])
        log_lon = np.logical_and(lon_l[0] < longitude, longitude < lon_l[-1])
        log = np.logical_and(log_lat,log_lon)

        # all lon/lat values inside bounding box 
        lons_sn = longitude[log]#.flatten()
        lats_sn = latitude[log]#.flatten()
        grid_points_sn = [Point(lon,lat) for lon,lat in zip(lons_sn,lats_sn)]

        # keep only values that lie inside the catchment shape:
        sn_points_within_polygon = [point for point in grid_points_sn if polygon_geometry.contains(point)]
        if sn_points_within_polygon:
            # use the above to generate the allcoords dataframe that contains all constant grid point
            # information for the catchments (such as lat/lon, elevation, ...)
            xcoor, ycoor = [], []
            xcoor_py, ycoor_py = [], []
            vect_ind, vect_ind_py = [], []
            for pt in sn_points_within_polygon:
                if pt.y in latitude and pt.x in longitude:
                    # get x and y indices for the desired lon/lat points:
                    ind_xy = np.where(np.logical_and(latitude == pt.y,longitude == pt.x)) # note that dimensions of the array are (lat,lon)
                    if len(ind_xy[0]) == 1 and len(ind_xy[1]) == 1:
                        xcoor_py.append(ind_xy[1][0])
                        xcoor.append(ind_xy[1][0] + 1) # add one to be consisten with R indexing which starts at 1 instead 0!
                        ycoor_py.append(ind_xy[0][0])
                        ycoor.append(ind_xy[0][0] + 1) # add one to be consisten with R indexing which starts at 1 instead 0!
                    else:
                        print('warning')

                    # get index of flattened field:
                    ind_vct = np.where(np.logical_and(latitude.flatten() == pt.y,longitude.flatten() == pt.x))
                    if len(ind_vct[0]) == 1:
                        vect_ind_py.append(ind_vct[0][0])
                        vect_ind.append(ind_vct[0][0] + 1) # add one to be consisten with R indexing which starts at 1 instead 0!
                    else:
                        print('warning')
        else: # case where no grid point lies inside catchment (often small intakes)
            centr = [[geom.centroid.x,geom.centroid.y] for geom in gdf.geometry][0]
        
            vect_ind_py = np.array([np.argmin((longitude-centr[0])**2 + (latitude-centr[1])**2)])
            lon_nearest,lat_nearest = longitude.flatten()[vect_ind_py],latitude.flatten()[vect_ind_py]
            sn_points_within_polygon = [Point(lon_nearest,lat_nearest)]    
            ycoor_py,xcoor_py = np.where(np.logical_and(longitude == lon_nearest,latitude==lat_nearest))
            xcoor = xcoor_py + 1
            ycoor = ycoor_py + 1
            vect_ind = vect_ind_py + 1

        elev = sn_geo.elevation.values.flatten()[vect_ind_py]
        utmx = sn_geo.easting.values.flatten()[xcoor_py]
        utmy = sn_geo.northing.values.flatten()[ycoor_py]
        plons = [pt.x.round(6) for pt in sn_points_within_polygon]
        plats = [pt.y.round(6) for pt in sn_points_within_polygon]

        # write to dataframe
        ctchment_coord_df = pd.DataFrame(
            dict(
                vect_ind = vect_ind,
                xcoor = xcoor,
                ycoor = ycoor,
                lon = plons,
                lat = plats,
                elev = elev,
                stat_id = ctchmnt,
                utmx = utmx,
                utmy = utmy
            )
        )

        ctchment_coord_df.to_csv(allcord_path, mode='a', header=not os.path.exists(allcord_path))

        # define slice to only load a small subset of the senorge data that contains the catchment
        bxx = slice(np.array(xcoor_py).min(),np.array(xcoor_py).max()+1)
        bxy = slice(np.array(ycoor_py).min(),np.array(ycoor_py).max()+1)

        # select the catchment grid points
        tg_catch = NC.variables['tg'][:,bxy,bxx][:,np.array(ycoor_py)-bxy.start,np.array(xcoor_py)-bxx.start]
        rr_catch = NC.variables['rr'][:,bxy,bxx][:,np.array(ycoor_py)-bxy.start,np.array(xcoor_py)-bxx.start]
        tn_catch = NC.variables['tn'][:,bxy,bxx][:,np.array(ycoor_py)-bxy.start,np.array(xcoor_py)-bxx.start]
        tx_catch = NC.variables['tx'][:,bxy,bxx][:,np.array(ycoor_py)-bxy.start,np.array(xcoor_py)-bxx.start]

        sn_daily = pd.DataFrame(
            {
                'date': time_pd.strftime('%Y-%m-%d'),
                'catchname' : ctchmnt,
                'prec' : np.round(np.nanmean(rr_catch-1),2),
                'st' : np.round(np.nanmean(tg_catch-1),2),
                'st_min' : np.round(np.nanmean(tn_catch-1),2),
                'st_max' : np.round(np.nanmean(tx_catch-1),2),
                'elev' : np.round(np.nanmean(elev),2)
            }
        )

        sn_daily.to_csv(daily_out, mode='a', header=not os.path.exists(daily_out))

        # get central lon/lat:
        centr = [[geom.centroid.x,geom.centroid.y] for geom in gdf.geometry][0]
        lonlat = pd.DataFrame(
            dict(
                lon = [centr[0]],
                lat = [centr[-1]],
            )
        )
        # write to file
        lonlat.to_csv(lola_out, mode='a', header = not os.path.exists(lola_out))