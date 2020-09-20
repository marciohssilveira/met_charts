import metpy.calc as mpcalc
import numpy as np
from metpy.units import units
import xarray as xr


class ExtractVariables:
    def __init__(self, data):
        self.data = data
        self.time_step = 0

        # Obtaining the index for the pressure levels of interest
        # isobaric variables
        index_lvl_1000 = np.where(np.array(self.data["isobaric"]) == 100000)[0][0]
        index_lvl_850 = np.where(np.array(self.data["isobaric"]) == 85000)[0][0]
        index_lvl_700 = np.where(np.array(self.data["isobaric"]) == 70000)[0][0]
        index_lvl_500 = np.where(np.array(self.data["isobaric"]) == 50000)[0][0]
        index_lvl_300 = np.where(np.array(self.data["isobaric"]) == 30000)[0][0]
        index_lvl_250 = np.where(np.array(self.data["isobaric"]) == 25000)[0][0]
        # isobaric6 variables
        index6_lvl_1000 = np.where(np.array(self.data["isobaric6"]) == 100000)[0][0]
        index6_lvl_850 = np.where(np.array(self.data["isobaric6"]) == 85000)[0][0]
        index6_lvl_700 = np.where(np.array(self.data["isobaric6"]) == 70000)[0][0]
        index6_lvl_500 = np.where(np.array(self.data["isobaric6"]) == 50000)[0][0]
        index6_lvl_300 = np.where(np.array(self.data["isobaric6"]) == 30000)[0][0]
        index6_lvl_250 = np.where(np.array(self.data["isobaric6"]) == 25000)[0][0]

        # Extracting data for the variables and levels of interest
        # Keeping them available throughout the class
        self.tair_850 = (self.data["Temperature_isobaric"][self.time_step][index6_lvl_850] - 273.15)
        self.tair_850.attrs["units"] = "degree_Celsius"
        self.rhum_850 = self.data["Relative_humidity_isobaric"][self.time_step][index_lvl_850]
        self.rhum_850.attrs["units"] = "percent"
        uwnd_850 = self.data["u-component_of_wind_isobaric"][self.time_step][index_lvl_850]
        vwnd_850 = self.data["v-component_of_wind_isobaric"][self.time_step][index_lvl_850]

        self.tair_700 = (self.data["Temperature_isobaric"][self.time_step][index6_lvl_700] - 273.15)
        self.tair_700.attrs["units"] = "degree_Celsius"
        self.rhum_700 = self.data["Relative_humidity_isobaric"][self.time_step][index_lvl_700]
        self.rhum_700.attrs["units"] = "percent"

        self.tair_500 = (self.data["Temperature_isobaric"][self.time_step][index6_lvl_500] - 273.15)
        self.tair_500.attrs["units"] = "degree_Celsius"
        uwnd_500 = self.data["u-component_of_wind_isobaric"][self.time_step][index_lvl_500]
        vwnd_500 = self.data["v-component_of_wind_isobaric"][self.time_step][index_lvl_500]

        # Calculate wind speed and direction using metpy functions
        self.wind_spd_850 = mpcalc.wind_speed(uwnd_850, vwnd_850)
        self.wind_dir_850 = mpcalc.wind_direction(uwnd_850, vwnd_850)
        self.wind_spd_500 = mpcalc.wind_speed(uwnd_500, vwnd_500)
        self.wind_dir_500 = mpcalc.wind_direction(uwnd_500, vwnd_500)
        
        # Calculate dewpoints using a metpy function
        dewp_850 = mpcalc.dewpoint_from_relative_humidity(self.tair_850, self.rhum_850)
        dewp_700 = mpcalc.dewpoint_from_relative_humidity(self.tair_700, self.rhum_700)
        # The results of those metpy functions are not xarray.Dataset objects
        self.dewp_850 = xr.DataArray(dewp_850)
        self.dewp_850.attrs["units"] = "degree_Celsius"
        self.dewp_700 = xr.DataArray(dewp_700)
        self.dewp_700.attrs["units"] = "degree_Celsius"