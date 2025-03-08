#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


class BiaxialConic(pet.component, family="pet.projections.biaxialProjections.biaxialConic",
                   implements=pet.protocols.projections.biaxialProjection):
    """
    Class that represents a biaxial conic projection of geodetic coordinates
    """

    central_latitude = pet.properties.float()
    central_latitude.default = 39
    central_latitude.doc = "Latitude of natural origin (lat_0)"

    central_longitude = pet.properties.float()
    central_longitude.default = -96
    central_longitude.doc = "Defines the latitude where scale is not distorted (lat_ts)"

    false_easting = pet.properties.float()
    false_easting.default = 0
    false_easting.doc = "False easting (x_0)"

    false_northing = pet.properties.float()
    false_northing.default = 0
    false_northing.doc = "False northing (y_0)"

    first_standard_parallel = pet.properties.float()
    first_standard_parallel.default = 33
    first_standard_parallel.doc = "First standard parallel latitudes"

    second_standard_parallel = pet.properties.float()
    second_standard_parallel.default = 45
    second_standard_parallel.doc = "Second standard parallel latitudes"

    cutoff = pet.properties.float()
    cutoff.default = -30
    cutoff.doc = "Latitude of map cutoff"

    north_extent = pet.properties.float()
    north_extent.default = 90
    north_extent.doc = "North extent of map"

    south_extent = pet.properties.float()
    south_extent.default = -90
    south_extent.doc = "South extent of map"

    west_extent = pet.properties.float()
    west_extent.default = -180
    west_extent.doc = "West extent of map"

    east_extent = pet.properties.float()
    east_extent.default = 180
    east_extent.doc = "East extent of map"

    folder_path = pet.properties.str()
    folder_path.doc = "Path to save file"

    @pet.export
    def proj(self, planet):
        """"
        Provides the cartopy proj
        :param planet: Planet to do projection on
        :return: fig, ax, globe to plot projection on
        """

        # Get the planet axes
        planet_axes = planet.get_axes()

        # Define Enceladus globe
        img_globe = ccrs.Globe(semimajor_axis=planet_axes[0].value, semiminor_axis=planet_axes[2].value,
                               ellipse=None)

        # Create a circular map using Cylindrical projection
        fig, ax = plt.subplots(subplot_kw={
            'projection': ccrs.LambertConformal(central_latitude=self.central_latitude,
                                                central_longitude=self.central_longitude,
                                                false_easting=self.false_easting,
                                                false_northing=self.false_northing,
                                                standard_parallels=[self.first_standard_parallel,
                                                                    self.second_standard_parallel],
                                                cutoff=self.cutoff,
                                                globe=img_globe)})

        # Set the extent
        ax.set_extent([self.west_extent, self.east_extent,
                       self.south_extent, self.north_extent], crs=ccrs.PlateCarree(globe=img_globe))

        # Add latitude and longitude lines
        gl = ax.gridlines(crs=ccrs.PlateCarree(globe=img_globe), linewidth=0.5, color='black', alpha=0.5,
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

        return fig, ax, img_globe

# end of file
