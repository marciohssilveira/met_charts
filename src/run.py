from datetime import datetime
import os
import sys
src_dir = os.path.join(os.getcwd(), '..', 'src')
sys.path.append(src_dir)

data_dir = os.path.join(os.getcwd(), '..', 'data')
sys.path.append(data_dir)

import time

import netCDF4
import pandas as pd
import xarray as xr

from d01_data.data import GetGFSData
from d03_visualisation.charts import CalculateCharts
import warnings
warnings.filterwarnings("ignore")

if __name__ == '__main__':
    gfs = GetGFSData() # instantiate the class containing the functions to download data
    data = gfs.get() # fetch data
    gfs.save(data) # save data to a local file in PATH (will replace the previous file)
    PATH = 'data/01_raw/data.nc'

    data = xr.open_dataset(PATH)

    start_time = time.time()

    for time_step in range(len(data[list(dict(data.dims).keys())[-1]])):
        print('Plotting charts')
        charts = CalculateCharts(data, time_step)
        charts.clouds_humidity()
        charts.showers_heat_humidity()
        charts.rain()
        charts.thunderstorm_showers()
        charts.storms()
        charts.hail()
        charts.instability()

    elapsed_time = time.time() - start_time
    print('Process done in {} seconds'.format(elapsed_time))
