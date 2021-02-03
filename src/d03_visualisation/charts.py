import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


from d02_processing.indices import CalculateIndices
from d02_processing.variables import ExtractVariables


class CalculateCharts:
    def __init__(self, data, time_step):
        self.data = data
        self.time_step = time_step
        self.output_path = 'images/charts'

        # Getting a time stamp to be used in the images
        date = pd.to_datetime(np.datetime_as_string(
            data[list(dict(self.data.dims).keys())[-1]][self.time_step].values, unit='h', timezone='UTC'))
        self.time_stamp = f'{date.day:02d}-{date.month:02d}-{date.year} {date.hour:02d}UTC'

        # Use the functions in Class ExtractVariables and CalculateIndices
        indices = CalculateIndices(self.data, self.time_step)
        variables = ExtractVariables(self.data, self.time_step)

        # Assign variables to the data creation functions on ExtractVariables class
        # Indices
        self.k_index = indices.k
        self.tt_index = indices.tt
        self.li_index = indices.li
        self.vt_index = indices.vt
        self.sweat_index = indices.sweat

        # Variables
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

        self.area_1 = {'SBSR': (-20.816111, -49.404722),
                       'SBAU': (-21.144167, -50.426389),
                       'SBAQ': (-21.804444, -48.140278),
                       'SBRP': (-21.136389, -47.776667),
                       'SBAX': (-19.560556, -46.965556),
                       'SBUL': (-18.883611, -48.225278),
                       'SBUR': (-19.764722, -47.966111),
                       'SBVG': (-21.588889, -45.473333),
                       'SBZM': (-21.513056, -43.173056),
                       'SBBH': (-19.851944, -43.950556),
                       'SBCF': (-19.624444, -43.971944),
                       'SBMK': (-16.706111, -43.821944),
                       'SBIP': (-19.470556, -42.488056),
                       'SBGV': (-18.896944, -41.986111),
                       'SBBR': (-15.871111, -47.918611),
                       'SBGO': (-16.632500, -49.221111),
                       'SBCN': (-17.724722, -48.610000),
                       'SWLC': (-17.834722, -50.956111),
                       'SBBW': (-15.860833, -52.389444),
                       'SNBR': (-12.079167, -45.009444),
                       'SBLE': (-12.482222, -41.276944),
                       'SBVC': (-14.907778, -40.914722),
                       'SNVB': (-13.296389, -38.992500),
                       'SDIY': (-12.200556, -38.900556),
                       'SBSV': (-12.908611, -38.322500),
                       'SBIL': (-14.815000, -39.033333),
                       'SNTF': (-17.524444, -39.668333),
                       'SBCR': (-19.011944, -57.671389),
                       'SBCG': (-20.469444, -54.670278)}

        self.area_2 = {'SBGR': (-23.435556, -46.473056),
                       'SBMT': (-23.506667, -46.633889),
                       'SBSP': (-23.626111, -46.656389),
                       'SBSJ': (-23.228889, -45.871111),
                       'SBTA': (-23.038889, -45.515833),
                       'SBST': (-23.928056, -46.299722),
                       'SBKP': (-23.006944, -47.134444),
                       'SBJD': (-23.181667, -46.943611),
                       'SBBP': (-22.979167, -46.537500),
                       'SBJH': (-23.426944, -47.165833),
                       'SDCO': (-23.483056, -47.486389),
                       'SBBU': (-22.343611, -49.053889),
                       'SBAE': (-22.157778, -49.068333),
                       'SBML': (-22.195556, -49.926944),
                       'SBDN': (-22.178333, -51.418889),
                       'SBTG': (-20.751389, -51.680278),
                       'SBDB': (-21.247222, -56.452500),
                       'SBDO': (-22.200556, -54.925556),
                       'SBPP': (-22.549722, -55.703056),
                       'SBLO': (-23.330278, -51.136667),
                       'SBMG': (-23.479444, -52.012222),
                       'SBTD': (-24.685278, -53.696389),
                       'SBCA': (-25.002222, -53.501944),
                       'SSGG': (-25.388333, -51.523611),
                       'SBPO': (-26.217222, -52.694444),
                       'SDAM': (-22.859167, -47.108056)}
        # Obtain coordinates
        self.lon_2d, self.lat_2d = variables.coordinates()

    def create_map(self):
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
        plt.figure(figsize=(10, 14))

        # Control the figures with GridSpec
        gs = gridspec.GridSpec(
            nrows=5, ncols=1, height_ratios=[4, .1, .1, .1, .1],
            top=0.99, bottom=0.01, left=0.05, right=0.95,
            hspace=0.4)

        # Add the map to the top row of the gs grid
        ax = plt.subplot(gs[0], projection=plotcrs)

        # Set the extent
        ax.set_extent([-35., -60., -10., -30.])
        # Add state/country boundaries to plot
        states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none')

        ax.add_feature(cfeature.COASTLINE, edgecolor='gray')
        ax.add_feature(cfeature.BORDERS, edgecolor='gray')
        ax.add_feature(states_provinces, edgecolor='gray')
        ax.gridlines(draw_labels=True, dms=True,
                     x_inline=False, y_inline=False)

        # Plot area_1 airports:
        lons_1 = [x[-1] for x in self.area_1.values()]
        lats_1 = [y[0] for y in self.area_1.values()]
        ax.scatter(lons_1, lats_1, transform=ccrs.PlateCarree(),
                   marker="*", c='mediumblue', zorder=99)

        # Plot area_2 airports:
        lons_2 = [x[-1] for x in self.area_2.values()]
        lats_2 = [y[0] for y in self.area_2.values()]
        ax.scatter(lons_2, lats_2, transform=ccrs.PlateCarree(),
                   marker="*", c='black', zorder=99)

        return ax, gs

    def clouds_humidity(self):
        """
        UR > 70% average(850, 700, 500) (blue contourf)
        UR > 70% average (1000, 850) (green countour)
        Wind > 5 m/s (arrows)
        1000-500mb thickness (red contour)
        Sea level pressure (black contour)
        """
        # Create figure/map/grids
        ax, gs = self.create_map()

        # UR > 70% average (1000, 850) (green countour)
        rhum_1 = np.mean([self.rhum(1000), self.rhum(850)], axis=0)
        # Use levels to only plot humidity above 70%
        levels = np.arange(70, 110, 10)
        # Put plot on ax = gs[0] (row 0)
        gs0 = ax.contourf(self.lon_2d, self.lat_2d, rhum_1, cmap='Greens',
                          transform=ccrs.PlateCarree(), alpha=0.5, levels=levels)
        # Put colorbar on gs[1] (row 1)
        gs1 = plt.subplot(gs[1])
        colorbar = plt.colorbar(gs0,
                                cax=gs1,   # cax means where the colorbar is gonna be put on
                                orientation='horizontal',
                                extendrect=True,
                                ticks=levels)
        colorbar.set_label('UR > 70% nos níveis 1000/850 hPa', size=8)
        colorbar.ax.tick_params(labelsize=8)

        # UR > 70% average(850, 700, 500) (blue contourf)
        rhum_2 = np.mean(
            [self.rhum(850), self.rhum(700), self.rhum(500)], axis=0)
        # Put another plot on ax = gs[0] (row 0)
        gs0 = ax.contourf(self.lon_2d, self.lat_2d, rhum_2, cmap='Blues',
                          transform=ccrs.PlateCarree(), alpha=0.5, levels=levels)
        # Put colorbar on gs[2] (row 2)
        gs2 = plt.subplot(gs[2])
        colorbar = plt.colorbar(gs0,
                                cax=gs2,  # cax means where the colorbar is gonna be put on
                                orientation='horizontal',
                                extendrect=True,
                                ticks=levels)
        colorbar.set_label('UR > 70% nos níveis 850/700/500 hPa', size=8)
        colorbar.ax.tick_params(labelsize=8)

        # Wind > 5m/s (barbs)
        # Conditions
        condition = self.wnd_spd(850) > 5
        u_wnd = self.u_wnd(850) * condition
        v_wnd = self.v_wnd(850) * condition
        u_wnd[u_wnd == 0] = np.nan
        v_wnd[v_wnd == 0] = np.nan
        # Create plot
        # As the windspeed data are in m/s, the barb increments thresholds are adapted to look like knots
        ax.barbs(self.lon_2d, self.lat_2d, u_wnd, v_wnd, length=5,
                 regrid_shape=25, transform=ccrs.PlateCarree(),
                 barbcolor='dimgray', alpha=0.5,
                 barb_increments={'half': 2.57222, 'full': 5.14444, 'flag': 25.7222})

        # 1000-500mb thickness (red contour)
        hgpt = self.hgpt(500) - self.hgpt(1000)
        # Create plot
        levels = np.arange(hgpt.min(), hgpt.max(),
                           (hgpt.max() - hgpt.min()) / 10)
        # Put another plot on ax = gs[0] (row 0)
        cs = ax.contour(self.lon_2d, self.lat_2d, hgpt, colors='red',
                        transform=ccrs.PlateCarree(), levels=levels, alpha=0.6)
        ax.clabel(cs, inline=True, fontsize=8, fmt='%0.0f')

        # Sea level pressure (black contour)
        # Create plot
        levels = np.arange(self.mslp().min(), self.mslp().max(),
                           (self.mslp().max() - self.mslp().min()) / 10)
        cs = ax.contour(self.lon_2d, self.lat_2d, self.mslp(), colors='black',
                        transform=ccrs.PlateCarree(), alpha=0.6, levels=levels)
        ax.clabel(cs, inline=True, fontsize=8, fmt='%0.0f')

        # # Insert info on other elements of the plot on row 4
        # label = 'teste\n teste'
        # gs3 = plt.subplot(gs[3])
        # ax.annotate(label,
        #             xy=(0.5, 0), xytext=(0, 10),
        #             xycoords=('axes fraction', 'figure fraction'),
        #             textcoords='offset points',
        #             size=14, ha='center', va='bottom')

        # Create title
        ax.set_title(
            f'Umidade e/ou nebulosidade {self.time_stamp}', fontsize=16, ha='center')
        plt.savefig(f'{self.output_path}/umidade_{self.time_step:02d}.jpg')

    def showers_heat_humidity(self):
        """
        K > 30 + TTS > 45 (green countourf)
        K > 30 + TTS > 45 + LI < -1 (red countourf)
        LIFT (blue contour)
        300hPa geopotential height (black contour)
        """
        # Create figure/map
        ax, gs = self.create_map()

        # K > 30 + TTS > 45 (green countourf)
        # Define conditions
        condition_1 = (self.k_index() > 30) & (self.tt_index() > 45)
        combination_1 = (self.k_index() + self.tt_index()) * condition_1
        combination_1[combination_1 == 0] = np.nan
        # Put plot on ax = gs[0] (row 0)
        levels = np.arange(np.nanmin(combination_1), np.nanmax(combination_1) + 10,
                           (np.nanmax(combination_1) - np.nanmin(combination_1)) / 10)
        gs0 = ax.contourf(self.lon_2d, self.lat_2d, combination_1,
                          cmap='Greens', transform=ccrs.PlateCarree(),
                          alpha=0.3, levels=levels, extend='max')
        # Put colorbar on gs[1] (row 1)
        gs1 = plt.subplot(gs[1])
        colorbar = plt.colorbar(gs0,
                                cax=gs1,   # cax means where the colorbar is gonna be put on
                                orientation='horizontal',
                                extendrect=True,
                                ticks=[int(x) for x in levels])
        colorbar.set_label('Combinação de índices: K > 30 e TTS > 45', size=8)
        colorbar.ax.tick_params(labelsize=8)

        # # K > 30 + TTS > 45 + LI < -1 (red countourf)
        # Define conditions
        condition_2 = (self.k_index() > 30) & (
            self.tt_index() > 45) & (self.li_index() < -1)
        combination_2 = (self.k_index() + self.tt_index() +
                         self.li_index()) * condition_2
        combination_2[combination_2 == 0] = np.nan
        # Put another plot on ax = gs[0] (row 0)
        levels = np.arange(np.nanmin(combination_2), np.nanmax(combination_2) + 10,
                           (np.nanmax(combination_2) - np.nanmin(combination_2)) / 10)
        gs0 = ax.contourf(self.lon_2d, self.lat_2d, combination_2,
                          cmap='Reds', transform=ccrs.PlateCarree(),
                          alpha=0.3, levels=levels, extend='max')
        # Put colorbar on gs[2] (row 2)
        gs2 = plt.subplot(gs[2])
        colorbar = plt.colorbar(gs0,
                                cax=gs2,  # cax means where the colorbar is gonna be put on
                                orientation='horizontal',
                                extendrect=True,
                                ticks=[int(x) for x in levels])
        colorbar.set_label(
            'Combinação de índices: K > 30, TTS > 45 e LI < -1', size=8)
        colorbar.ax.tick_params(labelsize=8)

        # LIFT (blue contour)
        levels = np.arange(self.li_index().min(), self.li_index().max(), 2)
        # Create plot
        cs = ax.contour(self.lon_2d, self.lat_2d, self.li_index(),
                        colors='royalblue', transform=ccrs.PlateCarree(),
                        levels=levels[levels < -1], linestyles='solid', alpha=0.6)
        ax.clabel(cs, inline=True, fontsize=8, fmt='%0.0f')

        # 300hPa geopotential height (black contour)
        cint = np.arange(self.hgpt(300).min(), self.hgpt(300).max(), 50)
        # Create plot
        cs = ax.contour(self.lon_2d, self.lat_2d, self.hgpt(300),
                        colors='black', transform=ccrs.PlateCarree(),
                        alpha=0.5, levels=cint)
        ax.clabel(cs, inline=True, fontsize=8, fmt='%0.0f')
        ax.set_title(
            f'Pancadas por calor e umidade {self.time_stamp}', fontsize=16, ha='center')
        plt.savefig(f'{self.output_path}/pancadas_{self.time_step:02d}.jpg')

    def rain(self):
        """
        OMEGA (500hPa) < -0.001 (green contourf)
        OMEGA (500hPa) < -0.01 and UR > 40% average(1000/850) (orange contourf)
        OMEGA (500hPa) < -0.5 and UR > 70% average(1000/850/700/500) (red contourf)
        500hPa geopotential height (black contour)
        500hPa Streamlines (gray streamlines)
        """
        # Create figure/map
        ax, gs = self.create_map()

        # OMEGA < -0.001 (green contourf)
        # Define conditions
        condition_1 = self.omega(500) < -0.001
        combination_1 = self.omega(500) * condition_1
        combination_1[combination_1 == 0] = np.nan
        # Create plot
        gs0 = ax.contourf(self.lon_2d, self.lat_2d, combination_1, colors='palegreen',
                          transform=ccrs.PlateCarree(), level=1)
        # Put colorbar on gs[1] (row 1)
        gs1 = plt.subplot(gs[1])
        colorbar = plt.colorbar(gs0,
                                cax=gs1,   # cax means where the colorbar is gonna be put on
                                orientation='horizontal',
                                extendrect=False,
                                ticks=[])
        colorbar.set_label('Omega < -0.001 Pa/s em 500 hPa', size=8)
        colorbar.ax.tick_params(labelsize=8)

        # OMEGA < -0.001 and UR > 40% average(850/700/500) (orange contourf)
        # Define conditions
        rhum = np.mean(
            [self.rhum(850), self.rhum(700), self.rhum(500)], axis=0)
        condition_2 = (self.omega(500) < -0.001) & (rhum > 40)
        combination_2 = (self.omega(500) + rhum) * condition_2
        combination_2[combination_2 == 0] = np.nan
        # Create plot
        gs0 = ax.contourf(self.lon_2d, self.lat_2d, combination_2, colors='gold',
                          transform=ccrs.PlateCarree(), levels=1)
        # Put colorbar on gs[2] (row 2)
        gs2 = plt.subplot(gs[2])
        colorbar = plt.colorbar(gs0,
                                cax=gs2,  # cax means where the colorbar is gonna be put on
                                orientation='horizontal',
                                extendrect=False,
                                ticks=[])
        colorbar.set_label(
            'Omega < -0.001 Pa/s em 500 hPa + UR > 40% nos níveis 850/700/500 hPa', size=8)
        colorbar.ax.tick_params(labelsize=8)

        # OMEGA < -0.3 and UR > 70% average(850/700/500) (red contourf)
        # Define conditions
        condition_3 = (self.omega(500) < -0.3) & (rhum > 70)
        combination_3 = (self.omega(500) + rhum) * condition_3
        combination_3[combination_3 == 0] = np.nan
        # Create plot
        gs0 = ax.contourf(self.lon_2d, self.lat_2d, combination_3, colors='red',
                          transform=ccrs.PlateCarree(), levels=1)
        # Put colorbar on gs[3] (row 3)
        gs3 = plt.subplot(gs[3])
        colorbar = plt.colorbar(gs0,
                                cax=gs3,  # cax means where the colorbar is gonna be put on
                                orientation='horizontal',
                                extendrect=False,
                                ticks=[])
        colorbar.set_label(
            'Omega < -0.3 Pa/s em 500 hPa + UR > 70% nos níveis 850/700/500 hPa', size=8)
        colorbar.ax.tick_params(labelsize=8)
        # 500hPa geopotential height (black contour)
        cs = ax.contour(self.lon_2d, self.lat_2d, self.hgpt(500),
                        colors='black', transform=ccrs.PlateCarree(),
                        alpha=0.5)
        ax.clabel(cs, inline=True, fontsize=8, fmt='%0.0f')

        # 500hPa Streamlines (gray streamlines)
        ax.streamplot(self.lon_2d, self.lat_2d, self.u_wnd(500), self.v_wnd(500),
                      transform=ccrs.PlateCarree(), color='slategray')

        ax.set_title(f'Chuva {self.time_stamp}', fontsize=16, ha='center')
        plt.savefig(f'{self.output_path}/chuva_{self.time_step:02d}.jpg')

    def thunderstorm_showers(self):
        """
        OMEGA -0.001 and UR > 40% average(850/700/500) and K >30 TTS>45 LIF < -1 (red contourf)
        OMEGA -0.001 and UR > 40% average(850/700/500) and K >30 TTS>45 (orange contourf)
        250hPa divergence (blue contourf)
        250hPa Streamlines (gray streamlines)
        """
        # Create figure
        ax, gs = self.create_map()

        # 250hPa divergence (blue contourf)
        # Put plot on ax = gs[0] (row 0)
        levels = np.arange(0, self.div(250).max(), 5)
        gs0 = ax.contourf(self.lon_2d, self.lat_2d, self.div(250), cmap='Blues',
                          transform=ccrs.PlateCarree(), alpha=0.7, levels=levels)
        # Put colorbar on gs[1] (row 1)
        gs1 = plt.subplot(gs[1])
        colorbar = plt.colorbar(gs0,
                                cax=gs1,   # cax means where the colorbar is gonna be put on
                                orientation='horizontal',
                                extendrect=True,
                                ticks=levels)
        colorbar.set_label('Divergência (10^6) em 250 hPa', size=8)
        colorbar.ax.tick_params(labelsize=8)

        # OMEGA (500hPa) -0.001 and UR > 40% average(850/700/500) and K >30 TTS>45 (orange contourf)
        rhum = np.mean(
            [self.rhum(850), self.rhum(700), self.rhum(500)], axis=0)
        # Define conditions
        condition_1 = (self.omega(500) < -0.001) & (rhum >
                                                    40) & (self.k_index() > 30) & (self.tt_index() > 45)
        combination_1 = (self.omega(500) + rhum +
                         self.k_index() + self.tt_index()) * condition_1
        combination_1[combination_1 == 0] = np.nan
        # Create plot
        gs0 = ax.contourf(self.lon_2d, self.lat_2d, combination_1, colors='orange',
                          transform=ccrs.PlateCarree(), level=1)
        # Put colorbar on gs[2] (row 2)
        gs2 = plt.subplot(gs[2])
        colorbar = plt.colorbar(gs0,
                                cax=gs2,   # cax means where the colorbar is gonna be put on
                                orientation='horizontal',
                                extendrect=False,
                                ticks=[])
        colorbar.set_label(
            'Omega < -0.001 Pa/s em 500 hPa, UR > 40% nos níveis 850/700/500 e combinação de índices K > 30 e TTS > 45', size=8)
        colorbar.ax.tick_params(labelsize=8)

        # OMEGA (500hPa) -0.001 and UR > 40% average(850/700/500) and K > 30 TTS > 45 LIF < -1 (red contourf)
        # Define conditions
        condition_2 = (self.omega(500) < -0.001) & (rhum > 40) & (self.k_index() >
                                                                  30) & (self.tt_index() > 45) & (self.li_index() < -1)
        combination_2 = (
            self.omega(500) + rhum + self.k_index() + self.tt_index() + self.li_index()) * condition_2
        combination_2[combination_2 == 0] = np.nan
        # Create plot
        gs0 = ax.contourf(self.lon_2d, self.lat_2d, combination_2, colors='red',
                          transform=ccrs.PlateCarree(), level=1)
        # Put colorbar on gs[3] (row 3)
        gs3 = plt.subplot(gs[3])
        colorbar = plt.colorbar(gs0,
                                cax=gs3,   # cax means where the colorbar is gonna be put on
                                orientation='horizontal',
                                extendrect=False,
                                ticks=[])
        colorbar.set_label(
            'Omega < -0.001 Pa/s em 500 hPa, UR > 40% nos níveis 850/700/500 e combinação de índices K > 30, TTS > 45 e LIF < -1', size=8)
        colorbar.ax.tick_params(labelsize=8)

        # 250hPa Streamlines (gray streamlines)
        ax.streamplot(self.lon_2d, self.lat_2d, self.u_wnd(250), self.v_wnd(250),
                      transform=ccrs.PlateCarree(), color='slategray')
        ax.set_title(f'Pancadas de chuva com trovoada {self.time_stamp}',
                     fontsize=16, ha='center')
        plt.savefig(f'{self.output_path}/trovoadas_{self.time_step:02d}.jpg')

    def storms(self):
        """
        OMEGA -0.001 and UR > 40% average(1000/500) and K >35 TTS>50 LIF < -4 (purple contourf)
        850hPa wind (green contourf)
        850hPa >15m/s wind (vector)
        250hPa jetstream (wind > 35m/s) (yellow contourf)
        Precipitable water 40-60mm (red contour)
        850hPa streamlines (gray streamlines)
        """
        # Create figure
        ax, gs = self.create_map()

        # OMEGA -0.001 and UR > 40% average(1000/500) and K > 35 TTS > 50 LIF < -4 (purple contourf)
        # Define conditions
        rhum = np.mean([self.rhum(1000), self.rhum(500)], axis=0)
        condition_1 = (self.omega(500) < -0.001) & (rhum > 40) & (self.k_index()
                                                                  > 35) & (self.tt_index() > 50) & (self.li_index() < -4)
        combination_1 = (self.omega(500) + rhum + self.k_index() +
                         self.tt_index() + self.li_index()) * condition_1
        combination_1[combination_1 == 0] = np.nan
        # Create plot
        # Put colorbar on gs[1] (row 1)
        gs0 = ax.contourf(self.lon_2d, self.lat_2d, combination_1, colors='darkviolet',
                          transform=ccrs.PlateCarree(), level=1)
        gs1 = plt.subplot(gs[1])
        colorbar = plt.colorbar(gs0,
                                cax=gs1,   # cax means where the colorbar is gonna be put on
                                orientation='horizontal',
                                extendrect=True,
                                ticks=[])
        colorbar.set_label(
            'Omega < -0.001 Pa/s em 500 hPa, UR > 40% nos níveis 1000/500 e combinação de índices K > 35, TTS > 50 e LIF < -4', size=8)
        colorbar.ax.tick_params(labelsize=8)

        # 850hPa wind (green contourf) did not plot it yet because i fount it unecessary. wait for team feedback

        # 850hPa >15kt wind (vector)
        # Conditions
        condition = self.wnd_spd(850) > 15
        u_wnd = self.u_wnd(850) * condition
        v_wnd = self.v_wnd(850) * condition
        u_wnd[u_wnd == 0] = np.nan
        v_wnd[v_wnd == 0] = np.nan
        # Create plot
        # As the windspeed data are in m/s, the barb increments thresholds are adapted to look like knots
        ax.barbs(self.lon_2d, self.lat_2d, u_wnd, v_wnd, length=5,
                 regrid_shape=25, transform=ccrs.PlateCarree(),
                 barbcolor='dimgray', alpha=0.5,
                 barb_increments={'half': 2.57222, 'full': 5.14444, 'flag': 25.7222})

        # Precipitable water 40-60mm (red contour)
        # Create plot
        levels = np.arange(40, 60, 5)
        # Put another plot on ax = gs[0] (row 0)
        gs0 = ax.contourf(self.lon_2d, self.lat_2d, self.prwt(),
                          cmap='Blues', transform=ccrs.PlateCarree(),
                          alpha=0.5, levels=levels, extend='max')
        # Put colorbar on gs[2] (row 2)
        gs2 = plt.subplot(gs[2])
        colorbar = plt.colorbar(gs0,
                                cax=gs2,  # cax means where the colorbar is gonna be put on
                                orientation='horizontal',
                                extendrect=True,
                                ticks=levels)
        colorbar.set_label('Água precipitável > 60 kg.m-2', size=8)
        colorbar.ax.tick_params(labelsize=8)

        # 250hPa jetstream (wind > 60kt) (yellow contourf)
        # Create plot
        # As the wind speed values are in m/s, we first need to convert to knots
        wnd_spd = self.wnd_spd(250) * 1.94384
        levels = np.arange(60, np.nanmax(wnd_spd) + 10, 10)
        # Put another plot on ax = gs[0] (row 0)
        gs0 = ax.contourf(self.lon_2d, self.lat_2d, wnd_spd, cmap='YlOrRd',
                          transform=ccrs.PlateCarree(), alpha=0.2, levels=levels)
        # Put colorbar on gs[3] (row 3)
        gs3 = plt.subplot(gs[3])
        colorbar = plt.colorbar(gs0,
                                cax=gs3,  # cax means where the colorbar is gonna be put on
                                orientation='horizontal',
                                extendrect=True,
                                ticks=levels)
        # The values in ticks are converted to knots
        colorbar.set_label(
            'Corrente de jato (Vento em 250 hPa > 60 nós)', size=8)
        colorbar.ax.tick_params(labelsize=8)
        # 850hPa streamlines (gray streamlines)
        ax.streamplot(self.lon_2d, self.lat_2d, self.u_wnd(850), self.v_wnd(850),
                      transform=ccrs.PlateCarree(), color='slategray')

        ax.set_title(f'Tempestades {self.time_stamp}',
                     fontsize=16, ha='center')
        plt.savefig(f'{self.output_path}/tempestades_{self.time_step:02d}.jpg')

    def hail(self):
        """
        OMEGA -0.001 and UR > 70% average(1000/500) and TTS > 52 VT > 25 SWEAT > 220 LIF < -2 (blue contourf)
        500 hPa temperature (black contour)
        850 hPa temperature (gray contour)
        500hPa OMEGA > -2 (orange contour)
        """
        # Create figure
        ax, gs = self.create_map()

        # OMEGA -0.001 and UR > 70% average(1000/500) and TTS > 50 VT > 25 SWEAT > 220 LIF < -2 (blue contourf)
        # Define conditions
        rhum = np.mean([self.rhum(1000), self.rhum(500)], axis=0)
        condition_1 = (self.omega(500) < -0.001) & (rhum > 70) & (self.tt_index() >
                                                                  52) & (self.vt_index() > 25) & (self.sweat_index() > 220) & (self.li_index() < -2)
        combination_1 = (
            self.omega(500) + rhum + self.tt_index() + self.vt_index() + self.sweat_index() + self.li_index()) * condition_1
        combination_1[combination_1 == 0] = np.nan
        # Create plot
        # Put colorbar on gs[1] (row 1)
        gs0 = ax.contourf(self.lon_2d, self.lat_2d, combination_1, colors='red',
                          transform=ccrs.PlateCarree(), level=1)
        gs1 = plt.subplot(gs[1])
        colorbar = plt.colorbar(gs0,
                                cax=gs1,   # cax means where the colorbar is gonna be put on
                                orientation='horizontal',
                                extendrect=False,
                                ticks=[])
        colorbar.set_label(
            'Omega < -0.001 Pa/s em 500 hPa, UR > 40% nos níveis 1000/500 e combinação de índices K > 30 e TTS > 45', size=8)
        colorbar.ax.tick_params(labelsize=8)

        # 500 hPa temperature (black contour)
        cs = ax.contour(self.lon_2d, self.lat_2d, self.tair(500),
                        colors='black', transform=ccrs.PlateCarree(),
                        linestyles='solid', level=0)
        ax.clabel(cs, inline=True, fontsize=8, fmt='%0.0f')

        # 850 hPa temperature (gray contour)
        levels = np.arange(self.tair(850).min(), self.tair(850).max() + 3, 3)
        gs0 = ax.contourf(self.lon_2d, self.lat_2d, self.tair(850), cmap='coolwarm',
                          transform=ccrs.PlateCarree(), alpha=0.3, levels=levels, extend='both')
        gs2 = plt.subplot(gs[2])
        colorbar = plt.colorbar(gs0,
                                cax=gs2,   # cax means where the colorbar is gonna be put on
                                orientation='horizontal',
                                extendrect=False,
                                ticks=[int(x) for x in levels])
        colorbar.set_label('Temperatura em 850 hPa', size=8)
        colorbar.ax.tick_params(labelsize=8)

        # 500hPa OMEGA < -0.3 (orange contourf)
        levels = np.arange(self.omega(500).min(), self.omega(500).max())
        gs0 = ax.contourf(self.lon_2d, self.lat_2d, self.omega(500),
                          cmap='autumn', transform=ccrs.PlateCarree(),
                          levels=levels[levels < -0.3], alpha=0.4)
        gs3 = plt.subplot(gs[3])
        colorbar = plt.colorbar(gs0,
                                cax=gs3,   # cax means where the colorbar is gonna be put on
                                orientation='horizontal',
                                extendrect=False,
                                ticks=[int(x) for x in levels[levels < -0.3]])
        colorbar.set_label('Omega em 500 hPa < 0.3 Pa/s', size=8)
        colorbar.ax.tick_params(labelsize=8)

        ax.set_title(f'Granizo {self.time_stamp}', fontsize=16, ha='center')
        plt.savefig(f'{self.output_path}/granizo_{self.time_step:02d}.jpg')

    def instability(self):
        """
        index() = ((K > 30) + (TT > 45) + (SWEAT > 220)) / 3
        if index() > 1, then unstable
        """
        # Create figure/map
        ax, gs = self.create_map()

        # Define conditions
        condition_1 = (self.k_index() > 30) & (
            self.tt_index() > 45) & (self.sweat_index() > 220)
        combination_1 = (self.k_index() + self.tt_index() +
                         self.sweat_index()) * condition_1
        combination_1[combination_1 == 0] = np.nan
        # Put plot on ax = gs[0] (row 0)
        levels = np.arange(np.nanmin(combination_1), np.nanmax(combination_1) + 10,
                           (np.nanmax(combination_1) - np.nanmin(combination_1)) / 10)
        gs0 = ax.contourf(self.lon_2d, self.lat_2d, combination_1,
                          cmap='Reds', transform=ccrs.PlateCarree(),
                          alpha=0.3, levels=levels, extend='max')
        # Put colorbar on gs[1] (row 1)
        gs1 = plt.subplot(gs[1])
        colorbar = plt.colorbar(gs0,
                                cax=gs1,   # cax means where the colorbar is gonna be put on
                                orientation='horizontal',
                                extendrect=True,
                                ticks=[int(x) for x in levels])
        colorbar.set_label(
            'Combinação de índices: K > 30, TTS > 45 e SWEAT > 220. Valores > 1 indicam instabilidade.', size=8)
        colorbar.ax.tick_params(labelsize=8)

        ax.set_title(
            f'Instabilidade {self.time_stamp}', fontsize=16, ha='center')
        plt.savefig(
            f'{self.output_path}/instabilidade_{self.time_step:02d}.jpg')

if __name__ == "__main__":
    pass
