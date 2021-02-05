import numpy as np
import pandas as pd
from d02_processing.variables import ExtractVariables
from d02_processing.indices import CalculateIndices
import json


class Tables:
    def __init__(self, data, time_step):
        self.data = data
        self.time_step = time_step
        self.output_path = 'images/tables'
        date = pd.to_datetime(
            np.datetime_as_string(data[list(dict(self.data.dims).keys())[-1]][self.time_step].values, unit='h',
                                  timezone='UTC'))
        self.time_stamp = f'{date.day:02d}-{date.month:02d}-{date.year} {date.hour:02d}UTC'
        self.levels = list(self.data['isobaric'].values / 100)[::-1]

        # Variables
        variables = ExtractVariables(self.data, self.time_step)
        self.tair = variables.temperature
        self.dewp = variables.dew_point
        self.rhum = variables.relative_humidity
        self.hgpt = variables.geopotential_height
        self.mslp = variables.mean_sea_level_pressure
        self.omega = variables.omega
        self.wnd_dir = variables.wind_direction
        self.wnd_spd = variables.wind_speed
        self.u_wnd = variables.u_wind
        self.v_wnd = variables.v_wind
        self.div = variables.divergence
        self.prwt = variables.precipitable_water

        # Indices
        indices = CalculateIndices(self.data, self.time_step)
        self.k_index = indices.k
        self.tt_index = indices.tt
        self.li_index = indices.li
        self.vt_index = indices.vt
        self.sweat_index = indices.sweat

        # Loading airports coordinates
        with open('data/02_external_data/airports.json') as json_file:
            self.airports = json.load(json_file)

    def extract_airport_data(self, index, lat, lon):
        i_lat = np.where(self.data['lat'] == round(lat * 4) / 4)[0][0]
        i_lon = np.where(self.data['lon'] - 360 == round(lon * 4) / 4)[0][0]
        airport_data = np.array(index)[i_lat][i_lon]
        return airport_data

    def create_table(self):
        """
        :param adjusted_data: receives the post processed dataframe
        :param valid:
        :return:
        """
        for area in self.airports:
            for airport in self.airports[area]:
                lon = self.airports[area][airport]['lon']
                lat = self.airports[area][airport]['lat']
                k_index = self.extract_airport_data(self.k_index(), lat, lon)
                print(f'{self.time_stamp} - {airport}: {k_index}')
                # table_output_path = f'{self.output_path}/{area}/'
                # Path(table_output_path).mkdir(parents=True, exist_ok=True)
                # plt.savefig(f'{table_output_path}/sounding_{airport}_{self.time_step}.jpg')

        return k_index

    # def calculate_indices(self, station_data):
    #     p = station_data['pressure'].values * units.hPa
    #     T = station_data['Temperature_isobaric'].values * units.degC
    #     Td = station_data['Dewpoint'].replace(
    #         np.nan, 0.0000001).values * units.degC
    #     prof = mpcalc.parcel_profile(p, T[0], Td[0])
    #     cape = mpcalc.cape_cin(pressure=p,
    #                            temperature=T,
    #                            dewpt=Td,
    #                            parcel_profile=prof)
    #     lcl_pressure, lcl_temperature = mpcalc.lcl(pressure=p[0],
    #                                                temperature=T[0],
    #                                                dewpt=Td[0])
    #     el_pressure, el_temperature = mpcalc.el(pressure=p[0],
    #                                             temperature=T[0],
    #                                             dewpt=Td[0])
    #     lfc_pressure, lfc_temperature = mpcalc.lfc(pressure=p[0],
    #                                                temperature=T[0],
    #                                                dewpt=Td[0])


if __name__ == "__main__":
    pass
