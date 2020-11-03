from datetime import datetime
import os
import time

import netCDF4
import pandas as pd
import xarray as xr

from data import GetGFSData
from charts import CalculateCharts
import warnings
warnings.filterwarnings("ignore")

if __name__ == '__main__':
    gfs = GetGFSData()
    data = gfs.get()
    PATH = 'data/data.nc'
    if os.path.exists(PATH):
        os.remove(PATH)
    data.to_netcdf(PATH)

    data = xr.open_dataset(PATH)

    start_time = time.time()
    print('Plotting charts')

    for time_step in range(len(data['time'])):
        charts = CalculateCharts(data, time_step)
        charts.clouds_humidity()
        charts.showers_heat_humidity()
        charts.rain()
        charts.thunderstorm_showers()
        charts.storms()
        charts.hail()
        charts.instability()

    elapsed_time = time.time() - start_time
    print(f'Process done in {elapsed_time} seconds')
