import metpy.calc as mpcalc
import numpy as np
from metpy.units import units
import xarray as xr


class ExtractVariables:
    def __init__(self, data):
        self.data = data
        self.time_step = 0

    def coordinates(self):
        """
        Extracts the lon-lat 2d meshgrid of the data for future plotting usage
        """
        lon = self.data['lon'].values
        lat = self.data['lat'].values
        lon_2d, lat_2d = np.meshgrid(lon, lat)
        return lon_2d, lat_2d

    def temperature(self, level):
        """
        Receives the integer value of the desired vertical pressure level
        to extract from data the temperature values in degree Celsius.
        In return you will have a xarray.Dataset for the desired time and level.
        """
        # Obtaining the index for the given pressure level
        index_level = np.where(
            np.array(self.data["isobaric6"]) == level * 100)[0][0]
        tair = self.data["Temperature_isobaric"][self.time_step][index_level] - 273.15
        tair.attrs["units"] = "degree_Celsius"
        return tair

    def relative_humidity(self, level):
        """
        Receives the integer value of the desired vertical pressure level
        to extract from data the relative humidity values in %.
        In return you will have a xarray.Dataset for the desired time and level.
        """
        # Obtaining the index for the given pressure level
        index_level = np.where(
            np.array(self.data["isobaric"]) == level * 100)[0][0]
        rhum = self.data["Relative_humidity_isobaric"][self.time_step][index_level]
        rhum.attrs["units"] = "percent"
        return rhum

    def dew_point(self, level):
        """
        Receives the integer value of the desired vertical pressure level
        to extract from datathe dew point temperature values in degree Celsius.
        It will use previously defined self.temperature() function 
        to extract temperature and use it in a metpy function to obtain dewpoint values.
        In return you will have a xarray.Dataset for the desired time and level.
        """
        # Obtaining the index for the given pressure level
        index_level = np.where(
            np.array(self.data["isobaric"]) == level * 100)[0][0]
        rhum = self.data["Relative_humidity_isobaric"][self.time_step][index_level]
        rhum.attrs["units"] = "percent"
        # Obtain air temperature data from the previously defined function
        tair = self.temperature(level)
        # Calculate dewpoints using a metpy function
        dewp = mpcalc.dewpoint_from_relative_humidity(tair, rhum)
        # Convert the result of the metpy function into a xarray.Dataset
        dewp = xr.DataArray(dewp)
        dewp.attrs["units"] = "degree_Celsius"
        return dewp

    def wind(self, level):
        """
        Receives the integer value of the desired vertical pressure level
        to extract from data a tuple containing values of wind direction and speed.
        In return you will have a pint.Quantity array for the desired time and level.
        """
        # Obtaining the index for the given pressure level
        index_level = np.where(
            np.array(self.data["isobaric"]) == level * 100)[0][0]
        uwnd = self.data["u-component_of_wind_isobaric"][self.time_step][index_level]
        vwnd = self.data["v-component_of_wind_isobaric"][self.time_step][index_level]
        # Calculate wind speed and direction using metpy functions
        wind_spd = mpcalc.wind_speed(uwnd, vwnd)
        wind_dir = mpcalc.wind_direction(uwnd, vwnd)
        return wind_dir, wind_spd

    def wind_components(self, level):
        """
        Receives the integer value of the desired vertical pressure level
        to extract from data a tuple containing u and v components of wind.
        In return you will have a xarray.Dataset for the desired time and level.
        """
        # Obtaining the index for the given pressure level
        index_level = np.where(
            np.array(self.data["isobaric"]) == level * 100)[0][0]
        uwnd = self.data["u-component_of_wind_isobaric"][self.time_step][index_level]
        vwnd = self.data["v-component_of_wind_isobaric"][self.time_step][index_level]
        return uwnd, vwnd

    def geopotential_height(self, level):
        """
        Receives the integer value of the desired vertical pressure level
        to extract from data the geopotential values in gpm.
        In return you will have a xarray.Dataset for the desired time and level.
        """
        # Obtaining the index for the given pressure level
        index_level = np.where(
            np.array(self.data["isobaric6"]) == level * 100)[0][0]
        hgpt = self.data["Geopotential_height_isobaric"][self.time_step][index_level]
        hgpt.attrs["units"] = "gpm"
        return hgpt

    def mean_sea_level_pressure(self):
        """
        Extract values of mean sea level pressure.
        In return you will have a xarray.Dataset for the desired time and level.
        """
        mslp = self.data["Pressure_reduced_to_MSL_msl"][self.time_step] / 100
        mslp.attrs["units"] = "hPa"
        return mslp
