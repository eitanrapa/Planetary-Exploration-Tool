#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import scipy.optimize
import numpy as np


class Sphere(pet.component, family="pet.shapes.sphere", implements=pet.protocols.shape):
    """
    A function that represents the sphere shape and its geometric functions
    """

    r = pet.properties.float()
    r.default = 1
    r.doc = "radius of the sphere"

    @pet.export
    def project(self, points):
        """
        Calculates the projection of a point on a sphere
        :param points: set of x, y, z points to use as a vector from the origin
        :return: set of x, y, z at the point the vector intersects the ellipsoid
        """

        # Define an ellipsoid function using the parametric versions of the position vector and set it equal to 0
        def f(t, x0, y0, z0): return (x0 ** 2 + y0 ** 2 + z0 ** 2) * t ** 2 + (
                (2 * x0) + (2 * y0) + (2 * z0)) * t - self.r ** 2

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

            # Calculate the height above the radius
            height = np.sqrt(x ** 2 + y ** 2 + z ** 2) - self.r

            # Raise exception if undefined
            if x == y == z == 0:
                raise Exception("Undefined")

            # Calculate the latitude
            if z > 0:
                lat = np.pi / 2 - np.arctan(np.sqrt(x ** 2 + y ** 2) / z)
            elif z < 0:
                lat = -np.pi / 2 - np.arctan(np.sqrt(x ** 2 + y ** 2) / z)
            else:
                lat = 0

            # Calculate the longitude
            if x < 0 <= y:
                long = np.arctan(y / x) + np.pi
            elif x < 0 and y < 0:
                long = np.arctan(y / x) - np.pi
            elif x == 0 and y > 0:
                long = np.pi / 2
            elif x == 0 and y < 0:
                long = - np.pi / 2
            else:
                long = np.arctan(y / x)

            yield long, lat, height

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

            # Calculate the x value
            x = (self.r + height) * np.cos(lat) * np.cos(long)

            # Calculate the y value
            y = (self.r + height) * np.cos(lat) * np.sin(long)

            # Calculate the z values
            z = (self.r + height) * np.sin(lat)

            yield x, y, z

        return

# end of file
