#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import scipy.optimize


class TriaxialEllipsoid(pet.component, family="pet.shapes.triaxialEllipsoid", implements=pet.protocols.shape):
    """
    A function that represents the triaxial ellipsoid shape and its geometric functions
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
    def project(self, points):
        """
        Calculates the projection of a point on the triaxial ellipsoid
        :param points: set of x, y, z points to use as a vector from the origin
        :return: set of x, y, z at the point the vector intersects the ellipsoid
        """

        # Define an ellipsoid function using the parametric versions of the position vector and set it equal to 0
        def f(t, x0, y0, z0): return ((x0 * t) ** 2 / self.a ** 2 +
                                      (y0 * t) ** 2 / self.b ** 2 + (z0 * t) ** 2 / self.c ** 2 - 1)

        # Create a generator
        for point in points:
            # Unpack the coordinate system
            x, y, z = point

            # Find the root, with arbitrary starting estimate
            starting_estimate = 1
            root = scipy.optimize.fsolve(f, starting_estimate, args=point)

            # return the intersecting point by plugging in the root
            yield x * root, y * root, z * root

        return

    @pet.export
    def geodetic(self, points):
        """
        A function that converts Cartesian points to geodetic in a triaxial ellipsoid
        :param points: set of x, y, z points
        :return: long, lat, height of the points
        """

        # Create a generator
        for point in points:
            # Unpack the coordinate system
            x, y, z = point

        return

    @pet.export
    def cartesian(self, points):
        """
        A function that converts geodetic points to Cartesian in a triaxial ellipsoid
        :param points: set of long, lat, height points
        :return: x, y, z of the points
        """

        # Create a generator
        for point in points:
            # Unpack the coordinate system
            long, lat, height = point

        return

# end of file
