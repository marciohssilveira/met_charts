import warnings
from d03_visualisation.charts import CalculateCharts
from d03_visualisation.soundings import Sounding
from d03_visualisation.tables import Tables
from d01_data.data import GetGFSData
import xarray as xr
import time
import os
import datetime
warnings.filterwarnings("ignore")

start_time = time.time()

# Obtaining data
PATH = 'data/01_raw/data.nc'

# Checking if local data is outdated
print('Checking local data...')
if not os.path.exists(PATH):
    gfs = GetGFSData()  # instantiate the class containing the functions to download data
    print('Fetching new data.')
    data = gfs.get()  # fetch data
    print('Saving new data.')
    gfs.save(data)  # use the gfs.save() function to export to a netcdf format

    elapsed_time = time.time() - start_time
    print(f'Data obtained in {elapsed_time} seconds')
else:
    delta = datetime.datetime.fromtimestamp(os.path.getmtime(PATH)) - datetime.datetime.today()
    if delta > datetime.timedelta(0.5):
        print('Outdated file found!\nDownloading new data...')
        gfs = GetGFSData()  # instantiate the class containing the functions to download data
        print('Fetching new data.')
        data = gfs.get()  # fetch data
        print('Saving new data.')
        gfs.save(data)  # use the gfs.save() function to export to a netcdf format

        elapsed_time = time.time() - start_time
        print(f'Data obtained in {elapsed_time} seconds')
    else:
        print('Recent data found!')


# Open data.nc using xarray
data = xr.open_dataset(PATH)

# Plotting graphics
print('Plotting images...')
start_time = time.time()

# Looping through the duration of the data in use
for time_step in range(len(data[list(dict(data.dims).keys())[-1]])):
    # Calculate and plot soundings
    sounding = Sounding(data, time_step)
    sounding.plot_skewt()
    # Calculate and plot charts
    charts = CalculateCharts(data, time_step)
    charts.clouds_humidity()
    charts.showers_heat_humidity()
    charts.rain()
    charts.thunderstorm_showers()
    charts.storms()
    charts.hail()
    charts.instability()
    # # Calculate tables
    # tables = Tables(data, time_step)
    # print(tables.create_table())


elapsed_time = time.time() - start_time
print(f'Process done in {elapsed_time} seconds!')

