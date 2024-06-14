#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


class BiaxialPlanar(pet.component, family="pet.projections.biaxialplanar", implements=pet.protocols.projection):
    """
    Class that represents a biaxial planar projection of geodetic coordinates
    """

    central_latitude = pet.properties.float()
    central_latitude.default = 0
    central_latitude.doc = "Latitude of natural origin (lat_0)"

    central_longitude = pet.properties.float()
    central_longitude.default = 0
    central_longitude.doc = "Defines the latitude where scale is not distorted (lat_ts)"

    false_easting = pet.properties.float()
    false_easting.default = 0
    false_easting.doc = "False easting (x_0)"

    false_northing = pet.properties.float()
    false_northing.default = 0
    false_northing.doc = "False northing (y_0)"

    @pet.export
    def proj(self, planet):

        # Get the planet axes
        planet_axes = planet.get_axes()

        # Define Enceladus globe
        img_globe = ccrs.Globe(semimajor_axis=planet_axes[0] * 1e3, semiminor_axis=planet_axes[2] * 1e3, ellipse=None)

        # Create a circular map using Cylindrical projection
        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Stereographic(central_latitude=self.central_latitude,
                                                                            central_longitude=self.central_longitude,
                                                                            false_easting=self.false_easting,
                                                                            false_northing=self.false_northing,
                                                                            globe=img_globe)})

        return fig, ax, img_globe

# end of file
