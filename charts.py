import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
import numpy as np
from metpy.units import units

import indices


class CalculateCharts:
    def __init__(self, data):
        self.data = data
        calculate_indices = indices.CalculateIndices(self.data)
        self.k = calculate_indices.k()
        self.tt = calculate_indices.tt()
        self.li = calculate_indices.li()
        self.sweat = calculate_indices.sweat()

    def clouds_humidity(self):
        """
        UR > 70% average(850, 700, 500) (blue contourf)
        UR > 70% average (1000, 850) (green countour)
        Wind > 5m/s (arrows)
        1000-500mb thickness (red contour)
        """

    def showers_heat_humidity(self):
        """
        K > 30 + TTS > 45 (green countourf)
        K > 30 + TTS > 45 + LI < -1 (red countourf)
        LIFT (blue contour)
        300hPa geopotential height (black contour)
        """

    def rain(self):
        """
        OMEGA -0.001 (green contourf)
        OMEGA -0.01 and UR > 40% average(1000/850) (orange contourf)
        OMEGA -0.5 and UR > 70% average(1000/850/700/500) (red contourf)
        500hPa geopotential height (black contour)
        500hPa Streamlines (gray streamlines)
        """

    def thunderstorm_showers(self):
        """
        OMEGA -0.001 and UR > 40% average(1000/500) and K >30 TTS>45 LIF < -1 (red contourf)
        OMEGA -0.001 and UR > 40% average(1000/500) and K >30 TTS>45 (orange contourf)
        250hPa divergence (blue contourf)
        250hPa Streamlines (gray streamlines)
        """

    def storms(self):
        """
        OMEGA -0.001 and UR > 40% average(1000/500) and K >35 TTS>50 LIF < -4 (purple contourf)
        850hPa wind (green contourf)
        850hPa >15m/s wind (vector)
        250hPa jetstream (wind > 35m/s) (yellow contourf)
        Precipitable water 40-60mm (red contour)
        850hPa streamlines (gray streamlines)
        """

    def hail(self):
        """
        OMEGA -0.001 and UR > 70% average(1000/500) and TTS>52 VT > 25 SWEAT>220 LIF < -2 (blue contourf)
        500 hPa temperature (black contour)
        850 hPa temperature (gray contour)
        500hPa OMEGA > -2 (orange contour)
        """

    def instability(self):
        """
        index = ((K > 30) + (TT > 45) + (SWEAT > 220)) / 3
        if index > 1, then unstable
        """

    def temperature_advection(self):
        # Pull out variables you want to use
        hght_var = self.data.variables['Geopotential_height_isobaric']
        temp_var = self.data.variables['Temperature_isobaric']
        u_wind_var = self.data.variables['u-component_of_wind_isobaric']
        v_wind_var = self.data.variables['v-component_of_wind_isobaric']
        lat_var = self.data.variables['lat']
        lon_var = self.data.variables['lon']

        # Get actual data values and remove any size 1 dimensions
        lat = lat_var[:].squeeze()
        lon = lon_var[:].squeeze()
        hght = hght_var[:].squeeze()
        temp = units.Quantity(temp_var[:].squeeze(), temp_var.units)
        u_wind = units.Quantity(u_wind_var[:].squeeze(), u_wind_var.units)
        v_wind = units.Quantity(v_wind_var[:].squeeze(), v_wind_var.units)

        lev_850 = np.where(self.data.variables['isobaric'][:] == 850 * 100)[0][0]
        hght_850 = hght[lev_850]
        temp_850 = temp[lev_850]
        u_wind_850 = u_wind[lev_850]
        v_wind_850 = v_wind[lev_850]

        # Combine 1D latitude and longitudes into a 2D grid of locations
        lon_2d, lat_2d = np.meshgrid(lon, lat)

        # Gridshift for barbs
        lon_2d[lon_2d > 180] = lon_2d[lon_2d > 180] - 360

        #########################################

        # Use helper function defined above to calculate distance
        # between lat/lon grid points
        dx, dy = mpcalc.lat_lon_grid_deltas(lon_var, lat_var)

        # Calculate temperature advection using metpy function
        adv = mpcalc.advection(temp_850, [u_wind_850, v_wind_850],
                               (dx, dy), dim_order='yx')

        # Smooth heights and advection a little
        # Be sure to only put in a 2D lat/lon or Y/X array for smoothing
        Z_850 = mpcalc.smooth_gaussian(hght_850, 2)
        adv = mpcalc.smooth_gaussian(adv, 2)

        # Set Projection of Data
        datacrs = ccrs.PlateCarree()

        # Set Projection of Plot
        plotcrs = ccrs.cartopy.crs.Mercator(central_longitude=0.0,
                                            min_latitude=-80.0,
                                            max_latitude=84.0,
                                            globe=None,
                                            latitude_true_scale=None,
                                            false_easting=0.0,
                                            false_northing=0.0,
                                            scale_factor=None)

        # Create new figure
        fig = plt.figure(figsize=(10, 12))
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,
                               hspace=0.01, wspace=0.01)

        # Add the map and set the extent
        ax = plt.subplot(gs[0], projection=plotcrs)
        # plt.title(f'850mb Temperature Advection for {time:%d %B %Y %H:%MZ}', fontsize=16)
        ax.set_extent([-30., -80., 0., -50.])

        # Plot Height Contours
        clev850 = np.arange(900, 3000, 30)
        cs = ax.contour(lon_2d, lat_2d, Z_850, clev850, colors='black', linewidths=1.5,
                        linestyles='solid', transform=datacrs)
        plt.clabel(cs, fontsize=10, inline=1, inline_spacing=10, fmt='%i',
                   rightside_up=True, use_clabeltext=True)

        # Plot Temperature Contours
        clevtemp850 = np.arange(-20, 20, 2)
        cs2 = ax.contour(lon_2d, lat_2d, temp_850.to(units('degC')), clevtemp850,
                         colors='green', linewidths=1.25, linestyles='dashed',
                         transform=datacrs)
        plt.clabel(cs2, fontsize=10, inline=1, inline_spacing=10, fmt='%i',
                   rightside_up=True, use_clabeltext=True)

        # Plot Colorfill of Temperature Advection
        cint = np.arange(-8, 9)
        cf = ax.contourf(lon_2d, lat_2d, 3 * adv.to(units('delta_degC/hour')), cint[cint != 0],
                         extend='both', cmap='bwr', transform=datacrs)
        cax = plt.subplot(gs[1])
        cb = plt.colorbar(cf, cax=cax, orientation='horizontal', extendrect=True, ticks=cint)
        cb.set_label(r'$^{o}C/3h$', size='large')

        # Add state/country boundaries to plot
        states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none')

        SOURCE = 'Natural Earth'
        LICENSE = 'public domain'

        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(states_provinces, edgecolor='gray')
        ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

        # Plot Wind Barbs
        ax.barbs(lon_2d, lat_2d, u_wind_850.magnitude, v_wind_850.magnitude,
                 length=6, regrid_shape=20, pivot='middle', transform=datacrs)

        plt.savefig('test.png')
        plt.show()
