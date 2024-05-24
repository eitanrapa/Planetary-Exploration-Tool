#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
from pyproj.crs.datum import CustomDatum, CustomEllipsoid
from pyproj import Transformer
from pyproj import CRS

class BiaxialConic(pet.component, family="pet.projections.biaxialconic", implements=pet.protocols.projection):
    """

    """

    latitude_first_parallel = pet.properties.float()
    latitude_first_parallel.default = 0
    latitude_first_parallel.doc = "First standard parallel (lat_1)"

    latitude_second_parallel = pet.properties.float()
    latitude_second_parallel.default = 0
    latitude_second_parallel.doc = "Second standard parallel (lat_2)"

    latitude_false_origin = pet.properties.float()
    latitude_false_origin.default = 0
    latitude_false_origin.doc = "Latitude of projection center (lat_0)"

    longitude_false_origin = pet.properties.float()
    longitude_false_origin.default = 0
    longitude_false_origin.doc = "Longitude of projection center (lon_0)"

    easting_false_origin = pet.properties.float()
    easting_false_origin.default = 0
    easting_false_origin = "False easting (x_0)"

    northing_false_origin = pet.properties.float()
    northing_false_origin.default = 0
    northing_false_origin = "False northing (y_0)"

    def __init__(self, name, locator, implicit, planet):
        super().__init__(name, locator, implicit)
        planet_axes = planet.get_axes()
        self.axes = planet.get_axes[0], planet_axes[2]

    # def transform(self, geodetic_coordinates):
    #     ell = CustomEllipsoid(semi_major_axis=self.axes[0], semi_minor_axis=self.axes[1])
    #     cd = CustomDatum(ellipsoid=ell)
    #     geodetic_crs = CRS.from_string("proj="
    #     transformer = Transformer.fr
    #     geodetic_crs = Proj(proj="latlong",datum=cd)
    #     albers_equal_area_proj = Proj(proj="aea", datum=cd)
    #
    #     x, y = transform(p1=geodetic_proj, p2=albers_equal_area_proj, x=geodetic_coordinates[:, 0], y=)
    #



# end of file
