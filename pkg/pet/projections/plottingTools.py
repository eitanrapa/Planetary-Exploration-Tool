#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import numpy as np
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource


def grid_plot(fig, ax, globe, geodetic_coordinates, north_extent=90, south_extent=-90, east_extent=180,
              west_extent=-180):

    ax.set_extent([west_extent, east_extent, south_extent, north_extent], crs=ccrs.PlateCarree(globe=globe))

    trimmed_coordinates = [(long, lat, height * 1e3) for lat, long, height in
                           geodetic_coordinates if south_extent < lat < north_extent and
                           west_extent < long < east_extent]
    longitudes, latitudes, heights = zip(*trimmed_coordinates)

    # Make a grid to down-sample the topographical map for visualization
    n = 1000
    extent = [-180, 180, -90, 90]
    grid_x, grid_y = np.mgrid[extent[0]:extent[1]:n * 1j, extent[2]:extent[3]:n * 1j]

    # Interpolation grid
    grid = griddata((longitudes, latitudes), heights, (grid_x, grid_y), method='cubic', fill_value=0)

    # Make a light-source for topographic hillshading
    cmap = plt.get_cmap('terrain')
    ls = LightSource(azdeg=315, altdeg=45)

    # Shade the grid with some arbitrary vertical exaggeration
    rgb = ls.shade(grid.T, cmap=cmap, blend_mode='soft', vert_exag=100)

    # Plot the grid and colorbar
    try:
        im = ax.imshow(rgb, extent=[extent[0], extent[1], extent[3], extent[2]],
                       transform=ccrs.PlateCarree(globe=globe), vmin=min(heights), vmax=max(heights), cmap=cmap)
    except ValueError:
        print("Exception: Try a smaller size plot")
        quit()

    # Add colorbar
    plt.colorbar(im, label="Heights [m]")

    # Add latitude and longitude lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(globe=globe), linewidth=0.5, color='black', alpha=0.5,
                      linestyle='--', draw_labels=True)
    gl.top_labels = True
    gl.left_labels = True
    gl.right_labels = True
    gl.xlines = True
    gl.xlocator = mticker.FixedLocator([-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180])
    gl.ylocator = mticker.FixedLocator(
        [-90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 8, 'color': 'gray'}
    gl.ylabel_style = {'size': 8, 'color': 'grey'}

    return fig, ax


def scatter_plot(fig, ax, globe, geodetic_coordinates, north_extent=90, south_extent=-90, east_extent=180,
                 west_extent=-180):

    ax.set_extent([west_extent, east_extent, south_extent, north_extent], crs=ccrs.PlateCarree(globe=globe))

    coordinates = [(long, lat, height * 1e3) for lat, long, height in geodetic_coordinates]
    longitudes, latitudes, heights = zip(*coordinates)

    # Plot points on the map
    ax.scatter(longitudes, latitudes, transform=ccrs.PlateCarree(globe=globe),
               color='red', marker='o', s=0.1, alpha=0.25)

    # Add latitude and longitude lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(globe=globe), linewidth=0.5, color='black', alpha=0.5,
                      linestyle='--', draw_labels=True)
    gl.top_labels = True
    gl.left_labels = True
    gl.right_labels = True
    gl.xlines = True
    gl.xlocator = mticker.FixedLocator([-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180])
    gl.ylocator = mticker.FixedLocator(
        [-90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 8, 'color': 'gray'}
    gl.ylabel_style = {'size': 8, 'color': 'grey'}

    return fig, ax

# end of file
