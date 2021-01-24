import warnings
from d03_visualisation.charts import CalculateCharts
from d03_visualisation.sounding import Sounding
from d01_data.data import GetGFSData
import xarray as xr
import time
import os

PATH = 'data/01_raw/data.nc'
warnings.filterwarnings("ignore")

if not os.path.exists(PATH):
    gfs = GetGFSData()  # instantiate the class containing the functions to download data
    try:  # Sometimes the available memory on the machine won't support handling more than 24 hours of data
        #    so I will use this workaround until i find a better way to deal with in memory NECDF4 data
        print('Trying 34 hours of data')
        data = gfs.get(n_hours=34)  # fetch data
    except:
        print('Failed... \nTrying 1 hours of data')
        data = gfs.get(n_hours=1)
    # save data to a local file in PATH (will replace the previous file)
    gfs.save(data)
    print('Done!')

data = xr.open_dataset(PATH)

start_time = time.time()

print('Plotting images')
for time_step in range(len(data[list(dict(data.dims).keys())[-1]])):
    charts = CalculateCharts(data, time_step)
    charts.clouds_humidity()
    charts.showers_heat_humidity()
    charts.rain()
    charts.thunderstorm_showers()
    charts.storms()
    charts.hail()
    charts.instability()
    sounding = Sounding(data, time_step)
    sounding.plot_skewt()

elapsed_time = time.time() - start_time
print('Process done in {} seconds'.format(elapsed_time))

print('Plotting soundings')
