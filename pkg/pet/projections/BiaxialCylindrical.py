#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


class BiaxialCylindrical(pet.component, family="pet.projections.biaxialcylindrical",
                         implements=pet.protocols.projection):
    """
    Class that represents a biaxial cylindrical projection of geodetic coordinates
    """

    central_longitude = pet.properties.float()
    central_longitude.default = 0
    central_longitude.doc = "Defines the latitude where scale is not distorted (lat_ts)"

    @pet.export
    def proj(self, planet):

        # Get the planet axes
        planet_axes = planet.get_axes()

        # Define Enceladus globe
        img_globe = ccrs.Globe(semimajor_axis=planet_axes[0] * 1e3, semiminor_axis=planet_axes[2] * 1e3, ellipse=None)

        # Create a circular map using Cylindrical projection
        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.LambertCylindrical(central_longitude=
                                                                                 self.central_longitude,
                                                                                 globe=img_globe)})
        return fig, ax, img_globe

# end of file
