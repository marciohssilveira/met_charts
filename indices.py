import metpy.calc as mpcalc
import numpy as np
from metpy.units import units


class CalculateIndices:

    def __init__(self, data):
        self.data = data

        # Obtain pressure levels (Temperature data has more vertical levels, so lets filter)
        pressure_var = data.variables['isobaric']
        pressure = units.Quantity(pressure_var[:].squeeze(), 'Pa')
        index_pressure = [int(index) for index, value in enumerate(pressure)]

        pressure6_var = data.variables['isobaric6']
        pressure6 = units.Quantity(pressure6_var[:].squeeze(), 'Pa')
        index_pressure6 = [int(index) for index, value in enumerate(pressure6) if value in pressure]

        # Pull out variables
        temp_var = self.data.variables['Temperature_isobaric'][0][index_pressure6]
        rhum_var = self.data.variables['Relative_humidity_isobaric'][0][index_pressure]
        uwnd_var = self.data.variables['u-component_of_wind_isobaric'][0][index_pressure]
        vwnd_var = self.data.variables['v-component_of_wind_isobaric'][0][index_pressure]

        # Get actual data values and assigning units
        tair = units.Quantity(temp_var[:].squeeze() - 273.15, 'degC')
        rhum = units.Quantity(rhum_var[:].squeeze(), 'percent')
        uwnd = units.Quantity(uwnd_var[:].squeeze(), 'knots')
        vwnd = units.Quantity(vwnd_var[:].squeeze(), 'knots')
        wind_spd = mpcalc.wind_speed(uwnd, vwnd)
        wind_dir = mpcalc.wind_direction(uwnd, vwnd)

        # Obtaining useful levels indexes
        lev_850 = np.where(self.data.variables['isobaric'][:] == 850 * 100)[0][0]
        lev_700 = np.where(self.data.variables['isobaric'][:] == 700 * 100)[0][0]
        lev_500 = np.where(self.data.variables['isobaric'][:] == 500 * 100)[0][0]

        # Filter the data indexing per level and make them available throughout the Class
        self.tair_850 = tair[lev_850]
        self.rhum_850 = rhum[lev_850]
        self.wind_spd_850 = wind_spd[lev_850]
        self.wind_dir_850 = wind_dir[lev_850]

        self.tair_700 = tair[lev_700]
        self.rhum_700 = rhum[lev_700]

        self.tair_500 = tair[lev_500]
        self.wind_spd_500 = wind_spd[lev_500]
        self.wind_dir_500 = wind_dir[lev_500]

        # Calculate dewpoints using a metpy function
        self.dewp_850 = mpcalc.dewpoint_from_relative_humidity(self.tair_850, self.rhum_850)
        self.dewp_700 = mpcalc.dewpoint_from_relative_humidity(self.tair_700, self.rhum_700)

    def k(self):
        # Calculate K-index
        k_index = (self.tair_850 - self.tair_500) + self.dewp_850 - (self.tair_700 - self.dewp_700)

        # Smooth the data
        k_index = mpcalc.smooth_gaussian(k_index, 2)

        return k_index

    def tt(self):
        # Calculate TT-index
        tt_index = (self.tair_850 - self.dewp_850) - (2 * self.tair_500)

        # Smooth the data
        tt_index = mpcalc.smooth_gaussian(tt_index, 2)

        return tt_index

    def li(self):
        # Extract LI-index
        li_index = self.data.variables['Best_4_layer_Lifted_Index_surface'][:].squeeze()

        # Smooth the data
        li_index = mpcalc.smooth_gaussian(li_index, 2)

        return li_index

    def sweat(self):
        """
        SWEAT = 12 Td850hPa + 20 (TTS -49) + 2V850hPa + V500hPa + 125 (cis + 0,2)
        onde V850hPa e V500hPa são a velocidade do vento em nós em 850 hPa e 500 hPa, respectivamente;
        cis = sen [direção (graus) V500hPa -V850hPa]
        """
        # Calculate SWEAT-index
        cis = np.sin((self.wind_dir_500.magnitude - self.wind_dir_850.magnitude))
        sweat_index = (12 * self.tair_850.magnitude) + (20 * (self.tt().magnitude - 49)) + (
                2 * self.wind_spd_850.magnitude) + self.wind_spd_500.magnitude + (125 * (cis * 0.2))

        # Smooth the data
        sweat_index = mpcalc.smooth_gaussian(sweat_index, 2)

        return sweat_index
