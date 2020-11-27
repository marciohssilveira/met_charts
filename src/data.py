import datetime as dt
import time
import os

import xarray as xr
from siphon.catalog import TDSCatalog
from xarray.backends import NetCDF4DataStore


class GetGFSData:
    def __init__(self):
        self.variables = ['Relative_humidity_isobaric',
                          'Temperature_isobaric',
                          'u-component_of_wind_isobaric',
                          'v-component_of_wind_isobaric',
                          'Best_4_layer_Lifted_Index_surface',
                          'Geopotential_height_isobaric',
                          'Precipitable_water_entire_atmosphere_single_layer',
                          'Pressure_reduced_to_MSL_msl',
                          'Vertical_velocity_pressure_isobaric']
        self.URL = 'http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml'
        self.dataset = 'Latest Collection for GFS Quarter Degree Forecast'

    def get(self):
        """
        :param coordinates: tuple like (lon, lat)
        :param variables: chosen list of variables based on the variables list for the dataset
        :param n_hours: number of hours for the prediction
        :return: a subset of the netCDF4 dataset based on the given coordinates and variables
        """
        ###########################################
        # First we construct a TDSCatalog instance using the url
        gfs_catalog = TDSCatalog(self.URL)

        # We see this catalog contains three datasets.
        # print(gfs_catalog.datasets)

        gfs_subset = gfs_catalog.datasets[self.dataset].subset()

        ###########################################
        # Define sub_point to proceed with the query
        query = gfs_subset.query()

        ###########################################
        # Then we construct a query asking for data corresponding to desired latitude and longitude and
        # for the time interval. We also ask for NetCDF version 4 data and choose the variables.
        # This request will return all vertical levels for a single point and for the time interval.
        # Note the string representation of the query is a properly encoded query string.
        # lonlat_box(west, east, south, north)
        query.lonlat_box(north=0, south=-40, east=-25, west=-70)
        now = dt.datetime.utcnow()
        n_hours = 24  # Sometimes, depending on the available RAM, it will give an error because it is trying to download a large .nc file
        query.time_range(now, now + dt.timedelta(hours=n_hours))
        query.accept('netcdf4')

        ###########################################
        # We'll be pulling out the variables we want to use in the future,
        # as well as the values of pressure levels.
        # To get the name of the correct variable we look at the 'variables' attribute on.
        # The last of the variables listed in `coordinates` is the pressure dimension.
        # print(gfs_subset.variables)

        query.variables(*self.variables)

        ###########################################
        # We now request data from the server using this query.
        print('Downloading new data...')
        start_time = time.time()

        raw_data = gfs_subset.get_data(query)

        elapsed_time = time.time() - start_time
        print(f'Process done in {elapsed_time} seconds')

        # We need the datastore so that we can open the existing netcdf dataset we downloaded
        dataset = xr.open_dataset(NetCDF4DataStore(raw_data))
        return dataset

    def save(self, data):
        """
        Takes the result of the get method, which is a xarray.Dataset 
        and stores it in a local file for debugging purposes
        """
        PATH = 'data/data.nc'
        if os.path.exists(PATH):
            os.remove(PATH)
        data.to_netcdf(PATH)


if __name__ == "__main__":
    pass
