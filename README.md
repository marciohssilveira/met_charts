# Meteorological charts
Inspired by the work of [Pinheiro et al.](https://www.scielo.br/pdf/rbmet/v29n2/a06v29n2.pdf), I created some personalized weather charts to help me with my forecasts.

[data.py](https://github.com/marciohssilveira/met_charts/blob/master/src/data.py) downloads data from the latest run of the [Global Forecast System (GFS)](https://thredds.ucar.edu/thredds/ncss/grib/NCEP/GFS/Global_0p5deg/Best/dataset.html);

The method ```get()``` inside the class ```GetGFSData``` in [data.py](https://github.com/marciohssilveira/met_charts/blob/master/src/data.py) downloads the GFS data for the predefined area of interest, converts into [xarray](http://xarray.pydata.org/en/stable/) and saves it in a big NetCDF file (data/data.nc).

[variables.py](https://github.com/marciohssilveira/met_charts/blob/master/src/variables.py) extracts the variables;

[indices.py](https://github.com/marciohssilveira/met_charts/blob/master/src/indices.py) calculate thermodynamic indices;

[charts.py](https://github.com/marciohssilveira/met_charts/blob/master/src/charts.py) perform the calculations; and

[run.py](https://github.com/marciohssilveira/met_charts/blob/master/src/plot_charts.py) plot the charts and organizes them in the folder [/images](https://github.com/marciohssilveira/met_charts/tree/master/img).

## Examples
    
![svg](https://github.com/marciohssilveira/met_charts/blob/master/images/charts/umidade_00.jpg)
    
![svg](https://github.com/marciohssilveira/met_charts/blob/master/images/charts/pancadas_00.jpg)
    
![svg](https://github.com/marciohssilveira/met_charts/blob/master/images/charts/chuva_00.jpg)
       
![svg](https://github.com/marciohssilveira/met_charts/blob/master/images/charts/trovoadas_00.jpg)
    
![svg](https://github.com/marciohssilveira/met_charts/blob/master/images/charts/tempestades_00.jpg)
        
![svg](https://github.com/marciohssilveira/met_charts/blob/master/images/charts/granizo_00.jpg)
      
![svg](https://github.com/marciohssilveira/met_charts/blob/master/images/charts/instabilidade_00.jpg)

![svg](https://github.com/marciohssilveira/met_charts/blob/master/images/soundings/area_1/SBBR/sounding_SBBR_00.jpg)
