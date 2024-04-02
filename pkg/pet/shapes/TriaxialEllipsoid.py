#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import scipy.optimize


class TriaxialEllipsoid(pet.component, family="pet.shapes.triaxialellipsoid", implements=pet.protocols.shape):
    """

    """

    a = pet.properties.float()
    a.default = 1
    a.doc = "first of the semiaxes"

    b = pet.properties.float()
    b.default = 1
    b.doc = "second of the semiaxes"

    c = pet.properties.float()
    c.default = 1
    c.doc = "third of the semiaxes"

    @pet.export
    def intersect(self, position):
        """
        Calculates the geocentric radius of a point on Enceladus.
        :param position: x, y, z points to use as a vector from the origin
        :return: x, y, z at the point the vector intersects the ellipsoid
        """
        # Unpack the position
        x, y, z = position

        # Define an ellipsoid function using the parametric versions of the position vector and set it equal to 0
        def f(t): return (x * t) ** 2 / self.a ** 2 + (y * t) ** 2 / self.b ** 2 + (z * t) ** 2 / self.c ** 2 - 1

        # Find the root
        root = scipy.optimize.fsolve(f, 1)

        # return the intersecting point by plugging in the root
        return x * root, y * root, z * root

    @pet.export
    def cartesian_to_geodetic(self, x, y, z):
        """
        """

        return []

    @pet.export
    def geodetic_to_cartesian(self, long, lat, height):
        """

        """
        return []

# end of file
