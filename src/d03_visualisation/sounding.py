import numpy as np
from metpy.plots import SkewT
import metpy.calc as mpcalc
import matplotlib.pyplot as plt
from metpy.units import units


class PlotSkew:

    def __init__(self):
        pass

    def plot_skewt(self, station_data):
        """
        :param adjusted_data: receives the post processed dataframe
        :param valid:
        :return:
        """

        # We will pull the data out of the example dataset into individual variables
        # and assign units.

        p = station_data['pressure'].values * units.hPa
        T = station_data['Temperature_isobaric'].values * units.degC
        Td = station_data['Dewpoint'].replace(
            np.nan, 0.0000001).values * units.degC
        u = station_data['u-component_of_wind_isobaric'].values * \
            units('meters / second').to('knots')
        v = station_data['v-component_of_wind_isobaric'].values * \
            units('meters / second').to('knots')

        # Create a new figure. The dimensions here give a good aspect ratio.
        fig = plt.figure(figsize=(12, 9))
        skew = SkewT(fig, rotation=45)

        # Plot the data using normal plotting functions, in this case using
        # log scaling in Y, as dictated by the typical meteorological plot
        skew.plot(p, T, 'r')
        skew.plot(p, Td, 'g')
        skew.plot_barbs(p, u, v)
        skew.ax.set_ylim(1020, 100)
        skew.ax.set_xlim(-40, 60)

        # Calculate LCL height and plot as black dot
        lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
        skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')

        # Calculate full parcel profile and add to plot as black line
        prof = mpcalc.parcel_profile(p, T[0], Td[0])
        skew.plot(p, prof, 'k', linewidth=2)

        # An example of a slanted line at constant T -- in this case the 0
        # isotherm
        skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)

        # Add the relevant special lines
        skew.plot_dry_adiabats()
        skew.plot_moist_adiabats()
        skew.plot_mixing_lines()

        skew.shade_cape(p, T, prof)
        skew.shade_cin(p, T, prof)

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
