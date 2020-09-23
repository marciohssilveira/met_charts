import metpy.calc as mpcalc
import numpy as np
from metpy.units import units
import xarray as xr

from variables import ExtractVariables


class CalculateIndices:
    def __init__(self, data):
        self.data = data
        self.time_step = 0

        # Use the functions in Class ExtractVariables
        # To create the useful data for indices calculations
        variables = ExtractVariables(self.data)

        self.tair_850 = variables.temperature(850)
        self.tair_700 = variables.temperature(700)
        self.tair_500 = variables.temperature(500)
        self.dewp_850 = variables.temperature(850)
        self.dewp_700 = variables.temperature(700)
        self.wind_spd_850, self.wind_dir_850 = variables.wind(850)
        self.wind_spd_500, self.wind_dir_500 = variables.wind(500)

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
        tt_index = (self.tair_850.values + self.dewp_850.values) - (2 * self.tair_500.values)

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
