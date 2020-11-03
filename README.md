# Meteorological charts
Inspired by the work of [Pinheiro et al.](https://www.scielo.br/pdf/rbmet/v29n2/a06v29n2.pdf), I created some personalized weather charts to help me with my forecasts.

[data.py](https://github.com/marciohssilveira/met_charts/blob/master/src/data.py) downloads data from the latest run of the [Global Forecast System (GFS)](https://thredds.ucar.edu/thredds/ncss/grib/NCEP/GFS/Global_0p5deg/Best/dataset.html);

[variables.py](https://github.com/marciohssilveira/met_charts/blob/master/src/variables.py) extracts the variables;

[indices.py](https://github.com/marciohssilveira/met_charts/blob/master/src/indices.py) calculate thermodynamic indices;

[charts.py](https://github.com/marciohssilveira/met_charts/blob/master/src/charts.py) perform the calculations; and

[plot_charts.py](https://github.com/marciohssilveira/met_charts/blob/master/src/plot_charts.py) plot the charts and organizes them in the folder [/img](https://github.com/marciohssilveira/met_charts/tree/master/img).

## Usage

The method ```get()``` inside the class ```GetGFSData``` in [data.py](https://github.com/marciohssilveira/met_charts/blob/master/src/data.py) downloads the GFS data for the predefined area of interest, converts into [xarray](http://xarray.pydata.org/en/stable/) and saves it in a big NetCDF file (data/data.nc).


```python
import matplotlib.pyplot as plt
import sys
sys.path.append('./src')
import xarray as xr
PATH = './data/data.nc'
data = xr.open_dataset(PATH)
```

Then we calculate the desired fields using the methods inside ```CalculateCharts``` in [charts.py](https://github.com/marciohssilveira/met_charts/blob/master/src/charts.py) 


```python
from charts import CalculateCharts
# Instantiate the class with the data and the timestep
time_step = range(len(data['time']))[0]
charts = CalculateCharts(data, time_step)

```


```python
charts.clouds_humidity()
plt.show()
```


    
![svg](https://github.com/marciohssilveira/met_charts/blob/master/img/umidade_00.jpg)
    



```python
charts.showers_heat_humidity()
plt.show()
```


    
![svg](https://github.com/marciohssilveira/met_charts/blob/master/img/pancadas_00.jpg)
    



```python
charts.rain()
plt.show()
```


    
![svg](https://github.com/marciohssilveira/met_charts/blob/master/img/chuva_00.jpg)
    



```python
charts.thunderstorm_showers()
plt.show()
```


    
![svg](https://github.com/marciohssilveira/met_charts/blob/master/img/trovoadas_00.jpg)
    



```python
charts.storms()
plt.show()
```


    
![svg](https://github.com/marciohssilveira/met_charts/blob/master/img/tempestades_00.jpg)
    



```python
charts.hail()
plt.show()
```


    
![svg](https://github.com/marciohssilveira/met_charts/blob/master/img/granizo_00.jpg)
    



```python
charts.instability()
plt.show()
```


    
![svg](https://github.com/marciohssilveira/met_charts/blob/master/img/instabilidade_00.jpg)
    



```python

```
