import metpy.calc as mpcalc
import numpy as np
from metpy.units import units
import xarray as xr


class ExtractVariables:
    def __init__(self, data):
        self.data = data
        self.time_step = 1

    def coordinates(self):
        """
        Extracts the lon-lat 2d meshgrid of the data for future plotting usage
        """
        # Extract lon and lat grid data values
        lon = self.data['lon'].values
        lat = self.data['lat'].values
        # create a meshgrid for future plotting
        lon_2d, lat_2d = np.meshgrid(lon, lat)
        return lon_2d, lat_2d

    def temperature(self, level):
        """
        Receives the integer value of the desired vertical pressure level
        to extract from data the temperature values in degree Celsius.
        In return you will have a np.array for the desired time and level.
        """
        # Obtaining the index for the given pressure level
        index_level = np.where(np.array(self.data["isobaric6"]) == level * 100)[0][0]
        # Extracting temperature data and converting to celsius
        tair = self.data["Temperature_isobaric"][self.time_step][index_level] - 273.15
        return np.array(tair)

    def relative_humidity(self, level):
        """
        Receives the integer value of the desired vertical pressure level
        to extract from data the relative humidity values in %.
        In return you will have a np.array for the desired time and level.
        """
        # Obtaining the index for the given pressure level
        index_level = np.where(np.array(self.data["isobaric"]) == level * 100)[0][0]
        # Extracting relative humidity data
        rhum = self.data["Relative_humidity_isobaric"][self.time_step][index_level]
        return np.array(rhum)

    def dew_point(self, level):
        """
        Receives the integer value of the desired vertical pressure level
        to extract from datathe dew point temperature values in degree Celsius.
        It will use previously defined self.temperature() function 
        to extract temperature and use it in a metpy function to obtain dewpoint values.
        In return you will have a np.array for the desired time and level.
        """
        # Obtaining the index for the given pressure level
        index_level = np.where(np.array(self.data["isobaric"]) == level * 100)[0][0]
        # Obtain relative humidity data for the given pressure level
        rhum = self.data["Relative_humidity_isobaric"][self.time_step][index_level]
        # Assigning units to make metpy function work
        rhum.attrs["units"] = "percent"
        # Obtain air temperature data for the given pressure level
        index_level = np.where(np.array(self.data["isobaric6"]) == level * 100)[0][0]
        tair = self.data["Temperature_isobaric"][self.time_step][index_level] - 273.15
        # Assigning units to make metpy function work
        tair.attrs["units"] = "degree_Celsius"
        # Calculate dewpoints using a metpy function
        dewp = mpcalc.dewpoint_from_relative_humidity(tair, rhum)
        return np.array(dewp)

    def wind_speed(self, level):
        """
        Receives the integer value of the desired vertical pressure level
        to extract from data a tuple containing values of wind speed.
        In return you will have a pint.Quantity array for the desired time and level.
        """
        # Obtaining the index for the given pressure level
        index_level = np.where(np.array(self.data["isobaric"]) == level * 100)[0][0]
        # Extracting wind components data
        uwnd = self.data["u-component_of_wind_isobaric"][self.time_step][index_level]
        vwnd = self.data["v-component_of_wind_isobaric"][self.time_step][index_level]
        # Calculate wind speed using metpy functions
        wind_spd = mpcalc.wind_speed(uwnd, vwnd)
        return np.array(wind_spd)

    def wind_direction(self, level):
        """
        Receives the integer value of the desired vertical pressure level
        to extract from data a tuple containing values of wind direction.
        In return you will have a pint.Quantity array for the desired time and level.
        """
        # Obtaining the index for the given pressure level
        index_level = np.where(np.array(self.data["isobaric"]) == level * 100)[0][0]
        # Extracting wind components data
        uwnd = self.data["u-component_of_wind_isobaric"][self.time_step][index_level]
        vwnd = self.data["v-component_of_wind_isobaric"][self.time_step][index_level]
        # Calculate wind direction using metpy function
        wind_dir = mpcalc.wind_direction(uwnd, vwnd)
        return np.array(wind_dir)

    def u_wind(self, level):
        """
        Receives the integer value of the desired vertical pressure level
        to extract from data a tuple containing u component of wind.
        In return you will have a np.array for the desired time and level.
        """
        # Obtaining the index for the given pressure level
        index_level = np.where(np.array(self.data["isobaric"]) == level * 100)[0][0]
        # Extracting u-wind component data
        uwnd = self.data["u-component_of_wind_isobaric"][self.time_step][index_level]
        return np.array(uwnd)

    def v_wind(self, level):
        """
        Receives the integer value of the desired vertical pressure level
        to extract from data a tuple containing v component of wind.
        In return you will have a np.array for the desired time and level.
        """
        # Obtaining the index for the given pressure level
        index_level = np.where(np.array(self.data["isobaric"]) == level * 100)[0][0]
        # Extracting v-wind component data
        vwnd = self.data["v-component_of_wind_isobaric"][self.time_step][index_level]
        return np.array(vwnd)

    def geopotential_height(self, level):
        """
        Receives the integer value of the desired vertical pressure level
        to extract from data the geopotential values in gpm.
        In return you will have a np.array for the desired time and level.
        """
        # Obtaining the index for the given pressure level
        index_level = np.where(np.array(self.data["isobaric6"]) == level * 100)[0][0]
        # Extracting geopotential height data
        hgpt = self.data["Geopotential_height_isobaric"][self.time_step][index_level]
        return np.array(hgpt)

    def mean_sea_level_pressure(self):
        """
        Extract values of mean sea level pressure.
        In return you will have a np.array for the desired time and level.
        """
        # Extract mean sea level pressure data and convert it into hPa
        mslp = self.data["Pressure_reduced_to_MSL_msl"][self.time_step] / 100
        return np.array(mslp)

    def omega(self, level):
        """
        Receives the integer value of the desired vertical pressure level
        to extract from data the omega values in Pa/s.
        In return you will have a np.array for the desired time and level.
        """
        # Obtaining the index for the given pressure level
        index_level = np.where(np.array(self.data["isobaric4"]) == level * 100)[0][0]
        # Extracting omega data
        omega = self.data["Vertical_velocity_pressure_isobaric"][self.time_step][index_level]
        return np.array(omega)

    def precipitable_water(self):
        """
        Extract values of precipitable water.
        In return you will have a np.array for the desired time and level.
        """
        # Extracting precipitable water data
        prwt = self.data["Precipitable_water_entire_atmosphere_single_layer"][self.time_step]
        return np.array(prwt)

    def divergence(self, level):
        """
        Uses a metpy function to calculate wind divergence in a given level.
        Uses predefined fuctions to obtain wind and grid data.
        Returns a np.array.
        """
        # Grab lat/lon values
        lat = self.data['lat'].values
        lon = self.data['lon'].values

        # Compute dx and dy spacing for use in divergence calculation
        dx, dy = mpcalc.lat_lon_grid_deltas(lon, lat)

        # Extract wind components
        uwnd = self.u_wind(level)
        vwnd = self.v_wind(level)

        # Use MetPy to compute the divergence for the given wind level
        div = mpcalc.divergence(uwnd, vwnd, dx, dy)
        return np.array(div * 100000)
