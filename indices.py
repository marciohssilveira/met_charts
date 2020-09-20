import metpy.calc as mpcalc
import numpy as np
from metpy.units import units
import xarray as xr


class CalculateIndices:
    def __init__(self, data):
        self.data = data
        self.time_step = 0

        # Obtaining the index for the pressure levels of interest
        # isobaric variables
        index_lvl_850 = np.where(np.array(self.data["isobaric"]) == 85000)[0][0]
        index_lvl_700 = np.where(np.array(self.data["isobaric"]) == 70000)[0][0]
        index_lvl_500 = np.where(np.array(self.data["isobaric"]) == 50000)[0][0]
        # isobaric6 variables
        index6_lvl_850 = np.where(np.array(self.data["isobaric6"]) == 85000)[0][0]
        index6_lvl_700 = np.where(np.array(self.data["isobaric6"]) == 70000)[0][0]
        index6_lvl_500 = np.where(np.array(self.data["isobaric6"]) == 50000)[0][0]

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

    def k(self):
        # Calculate K-index
        k_index = (
            (self.tair_850.values - self.tair_500.values)
            + self.dewp_850.values
            - (self.tair_700.values - self.dewp_700.values)
        )

        # Smooth the data
        k_index = mpcalc.smooth_gaussian(k_index, 2)

        return np.array(k_index)

    def tt(self):
        # Calculate TT-index
        tt_index = (self.tair_850.values - self.dewp_850.values) - (2 * self.tair_500.values)

        # Smooth the data
        tt_index = mpcalc.smooth_gaussian(tt_index, 2)

        return np.array(tt_index)

    def li(self):
        # Extract LI-index
        li_index = self.data["Best_4_layer_Lifted_Index_surface"][self.time_step]

        # Smooth the data
        li_index = mpcalc.smooth_gaussian(li_index, 2)

        return np.array(li_index)

    def sweat(self):
        """
        SWEAT = 12 Td850hPa + 20 (TTS -49) + 2V850hPa + V500hPa + 125 (cis + 0,2)
        onde V850hPa e V500hPa são a velocidade do vento em nós em 850 hPa e 500 hPa, respectivamente;
        cis = sen [direção (graus) V500hPa -V850hPa]
        """
        # Calculate SWEAT-index
        cis = np.sin((self.wind_dir_500.magnitude -
                      self.wind_dir_850.magnitude))
        sweat_index = ((12 * self.tair_850.values) + 
                       (20 * (self.tt() - 49)) + 
                       (2 * self.wind_spd_850.magnitude) + 
                       self.wind_spd_500.magnitude + 
                       (125 * (cis * 0.2)))

        # Smooth the data
        sweat_index = mpcalc.smooth_gaussian(sweat_index, 2)

        return np.array(sweat_index)