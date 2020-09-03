from datetime import datetime
from metpy.units import units
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from cftime import num2pydate
from siphon.catalog import TDSCatalog
import indices
import get_data
from metpy.calc import dewpoint_from_relative_humidity, parcel_profile, wind_speed, wind_direction
import metpy.calc as mpcalc


URL = 'http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/catalog.xml'
dataset = 'Latest Collection for GFS Quarter Degree Forecast'

variables = ['Relative_humidity_isobaric',
             'Temperature_isobaric',
             'u-component_of_wind_isobaric',
             'v-component_of_wind_isobaric',
             'Best_4_layer_Lifted_Index_surface',
             'Geopotential_height_isobaric',
             'Precipitable_water_entire_atmosphere_single_layer']

gfs_data = get_data.GetGFSData(variables)
data = gfs_data.get_data()

# Instantiate with the data the indices class containing
# all the functions to calculate the thermodynamic indexes
calculate_indices = indices.CalculateIndices(data)

k_index = calculate_indices.k()
tt_index = calculate_indices.tt()
li_index = calculate_indices.li()
sweat_index = calculate_indices.sweat()

# def find_time_var(var, time_basename='time'):
#     # Helper function for finding proper time variable
#     for coord_name in var.coordinates.split():
#         if coord_name.startswith(time_basename):
#             return coord_name
#     raise ValueError('No time variable found for ' + var.name)
#
#
# temp_var = data.variables['Temperature_isobaric']
# time_var = data.variables[find_time_var(temp_var)]
# lat_var = data.variables['lat']
# lon_var = data.variables['lon']
#
# # Get actual data values and remove any size 1 dimensions
# lat = lat_var[:].squeeze()
# lon = lon_var[:].squeeze()
# lon_2d, lat_2d = np.meshgrid(lon, lat)


# Convert number of hours since the reference time into an actual date
# time = num2pydate(time_var[:].squeeze(), time_var.units)


# cape = data.variables['Convective_available_potential_energy_surface'][:].squeeze()


# #
# def create_map():
#     # Set the projection information
#     proj = ccrs.PlateCarree(central_longitude=0.0, globe=None)
#     # Create a figure with an axes object on which we will plot. Pass the projection to that axes.
#     fig, ax = plt.subplots(subplot_kw=dict(projection=proj), figsize=(10, 12))
#
#     # Zoom in
#
#     ax.set_extent([-35., -60., -5., -30.])
#
#     # Add state/country boundaries to plot
#     states_provinces = cfeature.NaturalEarthFeature(
#         category='cultural',
#         name='admin_1_states_provinces_lines',
#         scale='50m',
#         facecolor='none')
#
#     ax.add_feature(cfeature.LAND)
#     ax.add_feature(cfeature.COASTLINE)
#     ax.add_feature(cfeature.BORDERS)
#     ax.add_feature(states_provinces, edgecolor='gray')
#     ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
#
#     return fig, ax
#
#
# def create_plot(variable):
#     # Get a new background map figure using the predefined function
#     fig, ax = create_map()
#     gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99, hspace=0.01, wspace=0.01)
#
#     # Plot Colorfill of k-index
#     variable_min = int(variable.magnitude.min())
#     variable_max = int(variable.magnitude.max())
#     cint = np.arange(variable_min, variable_max)
#     plot = ax.contourf(lon_2d, lat_2d, tt_index, levels=cint[cint > 60],
#                      extend='max', cmap='Blues', transform=ccrs.PlateCarree())
#     cax = plt.subplot(gs[1])
#     cb = plt.colorbar(k, cax=cax, orientation='horizontal', extendrect=True, ticks=cint)
#     cb.set_label(r'$^{o}C$', size='large')
#
#
# plt.savefig('kindx.png')
# plt.show()
