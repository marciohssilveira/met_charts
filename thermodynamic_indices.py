import metpy.calc as mpcalc
from metpy.calc import dewpoint_from_relative_humidity, cape_cin
import numpy as np
from metpy.units import units


class Indices:

    def __init__(self, data):
        self.data = data

        # Pull out variables
        relh_var = self.data.variables['Relative_humidity_isobaric']
        temp_var = self.data.variables['Temperature_isobaric']

        # Get actual data values and assigning units
        self.relh = units.Quantity(relh_var[:].squeeze(), 'percent')
        self.temp = units.Quantity(temp_var[:].squeeze() - 273.15, 'degC')

    def k_index(self):
        # Convert number of hours since the reference time into an actual date
        # time = num2pydate(time_var[:].squeeze(), time_var.units)

        lev_850 = np.where(self.data.variables['isobaric'][:] == 850 * 100)[0][0]
        relh_850 = self.relh[lev_850]
        temp_850 = self.temp[lev_850]
        dewp_850 = dewpoint_from_relative_humidity(temp_850, relh_850)

        lev_700 = np.where(self.data.variables['isobaric'][:] == 700 * 100)[0][0]
        relh_700 = self.relh[lev_700]
        temp_700 = self.temp[lev_700]
        dewp_700 = dewpoint_from_relative_humidity(temp_700, relh_700)

        lev_500 = np.where(self.data.variables['isobaric'][:] == 500 * 100)[0][0]
        temp_500 = self.temp[lev_500]

        # Calculate K-index
        kindx = (temp_850 - temp_500) + dewp_850 - (temp_700 - dewp_700)

        # Smooth the data
        kindx = mpcalc.smooth_gaussian(kindx, 2)

        return kindx

    # def cape(self):
