import numpy as np
import pandas as pd
from metpy.plots import SkewT
import metpy.calc as mpcalc
import matplotlib.pyplot as plt
from metpy.units import units
from d02_processing.variables import ExtractVariables
import json
from pathlib import Path


class Sounding:
    def __init__(self, data, time_step):
        self.data = data
        self.time_step = time_step
        self.output_path = '/images/soundings'
        date = pd.to_datetime(
            np.datetime_as_string(data[list(dict(self.data.dims).keys())[-1]][self.time_step].values, unit='h',
                                  timezone='UTC'))
        self.time_stamp = f'{date.day:02d}-{date.month:02d}-{date.year} {date.hour:02d}UTC'
        # Variables
        self.levels = list(self.data['isobaric'].values / 100)[::-1]
        variables = ExtractVariables(self.data, self.time_step)
        self.tair = variables.temperature
        self.dewp = variables.dew_point
        self.u_wnd = variables.u_wind
        self.v_wnd = variables.v_wind

        # Loading airports coordinates
        with open('data/02_external_data/airports.json') as json_file:
            self.airports = json.load(json_file)

    def create_profile_variable(self, variable, lat, lon):
        index_lat = np.where(self.data['lat'] == round(lat * 4) / 4)[0][0]
        index_lon = np.where(self.data['lon'] - 360 == round(lon * 4) / 4)[0][0]
        vertical_profile = {}
        for level in self.levels:
            vertical_profile[level] = np.array(variable(level))[index_lat][index_lon]
        return vertical_profile

    def plot_skewt(self):
        """
        :param adjusted_data: receives the post processed dataframe
        :param valid:
        :return:
        """
        for area in self.airports:
            for airport in self.airports[area]:
                lon = self.airports[area][airport]['lon']
                lat = self.airports[area][airport]['lat']
                pressure_levels = self.levels * units.hPa
                tair = list(self.create_profile_variable(self.tair, lat, lon).values()) * units.degC
                dewp = list(self.create_profile_variable(self.dewp, lat, lon).values()) * units.degC
                u_wnd = list(self.create_profile_variable(self.u_wnd, lat, lon).values()) * \
                        units('meters / second').to('knots')
                v_wnd = list(self.create_profile_variable(self.v_wnd, lat, lon).values()) * \
                        units('meters / second').to('knots')

                # Create a new figure. The dimensions here give a good aspect ratio.
                fig = plt.figure(figsize=(12, 9))
                skew = SkewT(fig, rotation=45)

                # Plot the data using normal plotting functions, in this case using
                # log scaling in Y, as dictated by the typical meteorological plot
                skew.plot(pressure_levels, tair, 'r')
                skew.plot(pressure_levels, dewp, 'g')
                skew.plot_barbs(pressure_levels, u_wnd, v_wnd)
                skew.ax.set_ylim(1020, 100)
                skew.ax.set_xlim(-40, 60)

                # Calculate LCL height and plot as black dot
                lcl_pressure, lcl_temperature = mpcalc.lcl(pressure_levels[0], tair[0], dewp[0])
                skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')

                # Calculate full parcel profile and add to plot as black line
                prof = mpcalc.parcel_profile(pressure_levels, tair[0], dewp[0])
                skew.plot(pressure_levels, prof, 'k', linewidth=2)

                # An example of a slanted line at constant T -- in this case the 0
                # isotherm
                skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)

                # Add the relevant special lines
                skew.plot_dry_adiabats()
                skew.plot_moist_adiabats()
                skew.plot_mixing_lines()

                skew.shade_cape(pressure_levels, tair, prof)
                skew.shade_cin(pressure_levels, tair, prof)
                plt.title(f'Perfil vertical (GFS) de {airport} valido para {self.time_stamp}', fontsize=16, ha='center')
                sounding_output_path = f'{self.output_path}/{area}/{airport}'
                Path(sounding_output_path).mkdir(parents=True, exist_ok=True)
                plt.savefig(f'{sounding_output_path}/sounding_{airport}_{self.time_step:02d}.jpg')

        return skew

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
