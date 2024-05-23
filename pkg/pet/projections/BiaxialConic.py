#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import pyproj


class BiaxialConic(pet.component, family="pet.projections.biaxialconic", implements=pet.protocols.projection):
    """

    """

    def __init__(self, name, locator, implicit, planet):
        super().__init__(name, locator, implicit)
        planet_axes = planet.get_axes()
        self.axes = planet.get_axes[0], planet_axes[2]

    def project(self, central_longitude=0.0, central_latitude=0.0, standard_parallels=(20.0, 50.0)):

        # Define Enceladus globe
        img_globe = ccrs.Globe(semimajor_axis=self.axes[0], semiminor_axis=self.axes[1], ellipse=None)

        # Create a circular map using a Conic projection
        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.AlbersEqualArea(central_longitude=central_longitude,
                                                                              central_latitude=central_latitude,
                                                                              standard_parallels=standard_parallels,
                                                                              globe=img_globe)})

        # Plot the circular map
        ax.set_global()

        # Zoom in on South Pole
        ax.set_extent([-180, 180, -90, -30], crs=ccrs.PlateCarree(globe=img_globe))

        # Plot points on the map
        ax.scatter(positions[:, 1], positions[:, 2], transform=ccrs.PlateCarree(globe=img_globe),
                   color='black', marker='o', s=0.1)

        # Add latitude and longitude lines
        gl = ax.gridlines(crs=ccrs.PlateCarree(globe=img_globe),
                          linewidth=0.5, color='black', alpha=0.5, linestyle='--', draw_labels=True)
        gl.top_labels = True
        gl.left_labels = True
        gl.right_labels = True
        gl.xlines = True
        gl.xlocator = mticker.FixedLocator([-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180])
        gl.ylocator = mticker.FixedLocator([-90, -80, -70, -60, -50, -40, -30, -20, -10, 0])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 8, 'color': 'gray'}
        gl.ylabel_style = {'size': 8, 'color': 'grey'}

        # Show the plot
        plt.savefig('/home/erapapor/kraken-bak/Enceladus-Exploration-Twin-files/coverage_maps/Orbit_map',
                    format='png', dpi=500)

# end of file
