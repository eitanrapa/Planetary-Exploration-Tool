#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet


class Conversions(pet.component):
    """

    """

    def cartesian(self, planet_axes, geodetic_coordinates):
        a, b, c = planet_axes
        te = pet.ext.pet.TriaxialEllipsoid(a=a, b=b, c=c)
        gp = pet.ext.pet.GeodeticPoint(latitude=geodetic_coordinates[0],
                                       longitude=geodetic_coordinates[1], height=geodetic_coordinates[2])
        cp = pet.ext.pet.CartesianPoint(x=0, y=0, z=0)
        pet.ext.pet.cartesian(te=te, gp=gp, cp=cp)

        return cp.x, cp.y, cp.z

    def geodetic(self, planet_axes, cartesian_coordinates, tol=1.0e-9):
        a, b, c = planet_axes
        te = pet.ext.pet.TriaxialEllipsoid(a=a, b=b, c=c)
        gp = pet.ext.pet.GeodeticPoint(latitude=0, longitude=0, height=0)
        cp = pet.ext.pet.CartesianPoint(x=cartesian_coordinates[0], y=cartesian_coordinates[1],
                                        z=cartesian_coordinates[2])

        iterations = pet.ext.pet.geodetic(te=te, cp=cp, gp=gp, tol=tol)
        return gp.latitude, gp.longitude, gp.height, iterations
# end of file
