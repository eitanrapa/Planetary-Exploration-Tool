#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import numpy as np
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource


class CartopyViz(pet.component):
    """
    Class that encapsulates the methods use to plot with Cartopy
    """

    folder_path = pet.properties.str()
    folder_path.doc = "Path to file save"

    north_extent = pet.properties.float()
    north_extent.default = 90
    north_extent.doc = "North extent of plot"

    south_extent = pet.properties.float()
    south_extent.default = -90
    south_extent.doc = "South extent of plot"

    east_extent = pet.properties.float()
    east_extent.default = 180
    east_extent.doc = "East extent of plot"

    west_extent = pet.properties.float()
    west_extent.default = -180
    west_extent.doc = "West extent of plot"

    def grid_plot(self, fig, ax, globe, geodetic_coordinates):
        """
        Plot a grid of values such as topography
        :param fig: Matplotlib fig
        :param ax: Matplotlib ax
        :param globe: Cartopy globe
        :param geodetic_coordinates: geodetic coordinates to plot
        """

        # Set the extent
        ax.set_extent([self.west_extent, self.east_extent,
                       self.south_extent, self.north_extent], crs=ccrs.PlateCarree(globe=globe))

        # Get reduced set of coordinates
        trimmed_coordinates = [(long, lat, height) for lat, long, height in
                               geodetic_coordinates if self.south_extent < lat < self.north_extent and
                               self.west_extent < long < self.east_extent]
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

        # Return with formatting
        return fig, ax, globe

    def scatter_plot(self, fig, ax, globe, geodetic_coordinates, color='black'):
        """
        Plot a scatter of values such as a swath
        :param fig: Matplotlib fig
        :param ax: Matplotlib ax
        :param globe: Cartopy globe
        :param geodetic_coordinates: geodetic coordinates to plot
        :param color: color of scatterplot
        """

        # Set the extent
        ax.set_extent([self.west_extent, self.east_extent,
                       self.south_extent, self.north_extent], crs=ccrs.PlateCarree(globe=globe))

        # Get the coordinates
        coordinates = [(long, lat, height) for lat, long, height in geodetic_coordinates]
        longitudes, latitudes, heights = zip(*coordinates)

        # Plot points on the map
        ax.scatter(longitudes, latitudes, transform=ccrs.PlateCarree(globe=globe),
                   color=color, marker='o', s=0.1, alpha=0.25)

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

        # Return with formatting
        return fig, ax, globe

# end of file
