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
        variables = ExtractVariables(self.data)

        # Assign variables to the data creation functions on ExtractVariables class
        self.tair = variables.temperature
        self.dewp = variables.dew_point
        self.wind_dir = variables.wind_direction
        self.wind_spd = variables.wind_speed

    def k(self):
        # Calculate K-index
        k_index = ((self.tair(850) - self.tair(500)) +
                   self.dewp(850) - (self.tair(700) - self.dewp(700)))

        # Smooth the data
        k_index = mpcalc.smooth_gaussian(k_index, 2)

        return k_index

    def tt(self):
        # Calculate TT-index
        tt_index = (self.tair(850) + self.dewp(850)) - (2 * self.tair(500))

        # Smooth the data
        tt_index = mpcalc.smooth_gaussian(tt_index, 2)

        return tt_index

    def vt(self):
        # Calculate VT-index
        vt_index = self.tair(850) - self.tair(500)

        return vt_index

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
        cis = np.sin((self.wind_dir(500) - self.wind_dir(850)))
        sweat_index = ((12 * self.tair(850)) + (20 * (self.tt() - 49)) +
                       (2 * self.wind_spd(850)) + self.wind_spd(500) + (125 * (cis * 0.2)))

        # Smooth the data
        sweat_index = mpcalc.smooth_gaussian(sweat_index, 2)

        return sweat_index
