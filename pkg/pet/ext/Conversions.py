#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved
import numpy as np
import pet


class Conversions(pet.component):
    """
    Class that encapsulates conversions between Cartesian and Geodetic coordinates
    """

    a = pet.properties.float()
    a.doc = "first semi-axis [m]"

    b = pet.properties.float()
    b.doc = "second semi-axis [m]"

    c = pet.properties.float()
    c.doc = "third semi-axis [m]"

    def _cartesian(self, geodetic_coordinates):
        """
        Python binding for conversion between geodetic and Cartesian coordinates in a triaxial ellipsoid
        :param geodetic_coordinates: latitude, longitude, and height of the point to be converted in [deg, deg, m]
        :return: Cartesian coordinates corresponding to the geodetic coordinates
        """

        # Define a TriaxialEllipsoid struct
        te = pet.ext.pet.TriaxialEllipsoid(a=self.a, b=self.b, c=self.c)

        # Create empty list of coordinates
        cartesian_coordinates = []

        # Iterate over the geodetic coordinates
        for i in range(len(geodetic_coordinates)):

            # Define a GeodeticPoint struct
            gp = pet.ext.pet.GeodeticPoint(latitude=geodetic_coordinates[i, 0],
                                           longitude=geodetic_coordinates[i, 1], height=geodetic_coordinates[i, 2])

            # Define an empty CartesianPoint struct to be modified
            cp = pet.ext.pet.CartesianPoint(x=0, y=0, z=0)

            # Call the C++ code
            pet.ext.pet.cartesian(te=te, gp=gp, cp=cp)

            # Append the new Cartesian coordinates to list
            cartesian_coordinates.append([cp.x, cp.y, cp.z])

        return np.asanyarray(cartesian_coordinates)

    def cartesian(self, geodetic_coordinates):
        """
        Helper function for _cartesian
        """

        # Make sure geodetic_coordinates is a numpy array
        geodetic_coordinates = np.asanyarray(geodetic_coordinates)

        # Make sure dimensions is 2
        if geodetic_coordinates.ndim == 1:

            geodetic_coordinates = np.asanyarray([geodetic_coordinates])

        # Call _Cartesian function for results
        return self._cartesian(geodetic_coordinates=geodetic_coordinates)

    def _geodetic(self, cartesian_coordinates, tol):
        """
        Python binding for conversion between Cartesian and geodetic coordinates in a triaxial ellipsoid
        :param cartesian_coordinates: x, y, z of the point to be converted in [m, m, m]
        :param tol: Error tolerance for convergence
        :return: Geodetic coordinates corresponding to the Cartesian coordinates plus number of iterations to result
        """

        # Define a TriaxialEllipsoid struct
        te = pet.ext.pet.TriaxialEllipsoid(a=self.a, b=self.b, c=self.c)

        # Create an empty list of coordinates
        geodetic_coordinates = []

        for i in range(len(cartesian_coordinates)):

            # Define an empty GeodeticPoint struct to be modified
            gp = pet.ext.pet.GeodeticPoint(latitude=0, longitude=0, height=0)

            # Define a CartesianPoint struct
            cp = pet.ext.pet.CartesianPoint(x=cartesian_coordinates[i, 0], y=cartesian_coordinates[i, 1],
                                            z=cartesian_coordinates[i, 2])

            # Call the C++ code and get the number of iterations
            iterations = pet.ext.pet.geodetic(te=te, cp=cp, gp=gp, tol=tol)

            # Append the new geodetic coordinates to list
            geodetic_coordinates.append([gp.latitude, gp.longitude, gp.height, iterations])

        return np.asanyarray(geodetic_coordinates)

    def geodetic(self, cartesian_coordinates, tol=1.0e-9):
        """
        Helper function for _geodetic
        """

        # Make sure cartesian_coordinates is a numpy array
        cartesian_coordinates = np.asanyarray(cartesian_coordinates)

        # Make sure dimensions is 2
        if cartesian_coordinates.ndim == 1:

            cartesian_coordinates = np.asanyarray([cartesian_coordinates])

        # Call _Geodetic function for results
        return self._geodetic(cartesian_coordinates=cartesian_coordinates, tol=tol)

# end of file
