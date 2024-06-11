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


class BiaxialConic(pet.component, family="pet.projections.biaxialconic", implements=pet.protocols.projection):
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

    first_standard_paralell = pet.properties.float()
    first_standard_paralell.default = 33
    first_standard_paralell.doc = "First standard parallel latitudes"

    second_standard_paralell = pet.properties.float()
    second_standard_paralell.default = 45
    second_standard_paralell.doc = "Second standard parallel latitudes"

    cutoff = pet.properties.float()
    cutoff.default = -30
    cutoff.doc = "Latitude of map cutoff"

    @pet.export
    def proj(self, planet):

        planet_axes = planet.get_axes()

        # Define Enceladus globe
        img_globe = ccrs.Globe(semimajor_axis=planet_axes[0] * 1e3, semiminor_axis=planet_axes[2] * 1e3, ellipse=None)

        # Create a circular map using Cylindrical projection
        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.LambertConformal(central_latitude=self.central_latitude,
                                                                               central_longitude=self.central_longitude,
                                                                               false_easting=self.false_easting,
                                                                               false_northing=self.false_northing,
                                                                               standard_parallels=
                                                                               [self.first_standard_paralell,
                                                                                self.second_standard_paralell],
                                                                               cutoff=self.cutoff,
                                                                               globe=img_globe)})

        return fig, ax, img_globe

# end of file
