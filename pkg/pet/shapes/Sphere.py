#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet


class Sphere(pet.component, family="pet.shapes.sphere", implements=pet.protocols.shape):
    """

    """

    r = pet.properties.float()
    r.default = 1
    r.doc = "radius of the sphere"

    @pet.export
    def project(self, cartesian):
        """
        Calculates the projection of a point on a sphere
        :param cartesian: x, y, z points to use as a vector from the origin
        :return: x, y, z at the point the vector intersects the ellipsoid
        """
        # Unpack the coordinate system
        x, y, z = cartesian

        return []

    @pet.export
    def geodetic(self, cartesian):
        """
        A function that converts between cartesian and geodetic coordinates in a triaxial ellipsoid
        """
        # Unpack the coordinate system
        x, y, z = cartesian

        return []

    @pet.export
    def cartesian(self, geodetic):
        """
        A function that converts between geodetic and cartesian coordinates in a triaxial ellipsoid
        """
        # Unpack the coordinate system
        long, lat, height = geodetic

        return []

# end of file
