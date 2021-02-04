import warnings
from d03_visualisation.charts import CalculateCharts
from d03_visualisation.soundings import Sounding
from d01_data.data import GetGFSData
import xarray as xr
import time
import os
import datetime
warnings.filterwarnings("ignore")


# Downloading data
PATH = 'data/01_raw/data.nc'
delta = datetime.datetime.fromtimestamp(os.path.getmtime(PATH)) - datetime.datetime.today()

start_time = time.time()
print('Checking data...')
# if not os.path.exists(PATH):
if delta > datetime.timedelta(0.5):
    print('Downloading new data...')
    gfs = GetGFSData()  # instantiate the class containing the functions to download data
    data = gfs.get(n_hours=34)  # fetch data
    gfs.save(data)
    print('Done!')

data = xr.open_dataset(PATH)
elapsed_time = time.time() - start_time
print(f'Data obtained in {elapsed_time} seconds')


# Plotting graphics
print('Plotting images...')
start_time = time.time()
for time_step in range(len(data[list(dict(data.dims).keys())[-1]])):
    sounding = Sounding(data, time_step)
    sounding.plot_skewt()
    charts = CalculateCharts(data, time_step)
    charts.clouds_humidity()
    charts.showers_heat_humidity()
    charts.rain()
    charts.thunderstorm_showers()
    charts.storms()
    charts.hail()
    charts.instability()

elapsed_time = time.time() - start_time
print(f'Process done in {elapsed_time} seconds!')

