from datetime import datetime
import os

import netCDF4
import pandas as pd
import xarray as xr

from data import GetGFSData
from charts import CalculateCharts

URL = 'http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml'
dataset = 'Latest Collection for GFS Quarter Degree Forecast'

variables = ['Relative_humidity_isobaric',
             'Temperature_isobaric',
             'u-component_of_wind_isobaric',
             'v-component_of_wind_isobaric',
             'Best_4_layer_Lifted_Index_surface',
             'Geopotential_height_isobaric',
             'Precipitable_water_entire_atmosphere_single_layer',
             'Pressure_reduced_to_MSL_msl',
             'Vertical_velocity_pressure_isobaric']

gfs = GetGFSData(variables)
data = gfs.get()
path = 'data.nc'
if os.path.exists(path):
    os.remove(path)
data.to_netcdf('data.nc')

data = xr.open_dataset('data.nc')

for time_step in range(len(data['time3'])):
    charts = CalculateCharts(data, time_step)
    charts.clouds_humidity()
    charts.showers_heat_humidity()
    charts.rain()
    charts.thunderstorm_showers()
    charts.storms()
    charts.hail()
    charts.instability()