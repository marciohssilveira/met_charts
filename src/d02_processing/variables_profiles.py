from types import LambdaType
import metpy.calc as mpcalc
import numpy as np
from metpy.units import units
import xarray as xr


class ExtractVariablesProfile:
    def __init__(self, raw_data, time_step):
        self.raw_data = raw_data
        self.time_step = time_step

    def point(self, data, lat, lon):
        index_lat = np.where(self.raw_data['lat'] == round(lat * 4)/4)[0][0]
        index_lon = np.where(self.raw_data['lon'] - 360 == round(lon * 4)/4)[0][0]
        return np.array(data)[index_lat][index_lon]

if __name__ == "__main__":
    pass
