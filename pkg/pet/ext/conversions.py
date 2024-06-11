#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved
import numpy as np
import pet


def cartesian(planet, geodetic_coordinates):
    """
    Python binding for conversion between geodetic and Cartesian coordinates in a triaxial ellipsoid
    :param planet: Target planet
    :param geodetic_coordinates: latitude, longitude, and height of the point to be converted
    :return: Cartesian coordinates corresponding to the geodetic coordinates
    """

    # Unpack the axes
    a, b, c = planet.get_axes()

    # Define a TriaxialEllipsoid struct
    te = pet.ext.pet.TriaxialEllipsoid(a=a, b=b, c=c)

    cartesian_coordinates = []
    geodetic_coordinates = np.asarray(geodetic_coordinates)

    if geodetic_coordinates.ndim == 1:

        # Define a GeodeticPoint struct
        gp = pet.ext.pet.GeodeticPoint(latitude=geodetic_coordinates[0],
                                       longitude=geodetic_coordinates[1], height=geodetic_coordinates[2])

        # Define an empty CartesianPoint struct to be modified
        cp = pet.ext.pet.CartesianPoint(x=0, y=0, z=0)

        # Call the C++ code
        pet.ext.pet.cartesian(te=te, gp=gp, cp=cp)

        return [cp.x, cp.y, cp.z]

    else:

        for i in range(len(geodetic_coordinates)):

            # Define a GeodeticPoint struct
            gp = pet.ext.pet.GeodeticPoint(latitude=geodetic_coordinates[i, 0],
                                           longitude=geodetic_coordinates[i, 1], height=geodetic_coordinates[i, 2])

            # Define an empty CartesianPoint struct to be modified
            cp = pet.ext.pet.CartesianPoint(x=0, y=0, z=0)

            # Call the C++ code
            pet.ext.pet.cartesian(te=te, gp=gp, cp=cp)

            cartesian_coordinates.append([cp.x, cp.y, cp.z])

    return cartesian_coordinates


def geodetic(planet, cartesian_coordinates, tol=1.0e-9):
    """
    Python binding for conversion between Cartesian and geodetic coordinates in a triaxial ellipsoid
    :param planet: Target planet
    :param cartesian_coordinates: x, y, z of the point to be converted
    :param tol: Error tolerance for convergence
    :return: Geodetic coordinates corresponding to the Cartesian coordinates plus number of iterations to result
    """

    # Unpack the axes
    a, b, c = planet.get_axes()

    # Define a TriaxialEllipsoid struct
    te = pet.ext.pet.TriaxialEllipsoid(a=a, b=b, c=c)

    geodetic_coordinates = []
    cartesian_coordinates = np.asarray(cartesian_coordinates)

    if cartesian_coordinates.ndim == 1:

        # Define an empty GeodeticPoint struct to be modified
        gp = pet.ext.pet.GeodeticPoint(latitude=0, longitude=0, height=0)

        # Define a CartesianPoint struct
        cp = pet.ext.pet.CartesianPoint(x=cartesian_coordinates[0], y=cartesian_coordinates[1],
                                        z=cartesian_coordinates[2])

        # Call the C++ code and get the number of iterations
        iterations = pet.ext.pet.geodetic(te=te, cp=cp, gp=gp, tol=tol)

        return [gp.latitude, gp.longitude, gp.height, iterations]

    else:

        for i in range(len(cartesian_coordinates)):

            # Define an empty GeodeticPoint struct to be modified
            gp = pet.ext.pet.GeodeticPoint(latitude=0, longitude=0, height=0)

            # Define a CartesianPoint struct
            cp = pet.ext.pet.CartesianPoint(x=cartesian_coordinates[i, 0], y=cartesian_coordinates[i, 1],
                                            z=cartesian_coordinates[i, 2])

            # Call the C++ code and get the number of iterations
            iterations = pet.ext.pet.geodetic(te=te, cp=cp, gp=gp, tol=tol)

            geodetic_coordinates.append([gp.latitude, gp.longitude, gp.height, iterations])

        return geodetic_coordinates

# end of file
