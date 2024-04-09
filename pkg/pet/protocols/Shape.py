#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet


class Shape(pet.protocol, family="pet.shapes"):
    """
    The abstract specification for shapes
    """

    @pet.provides
    def project(self, points):
        """
        A function to calculate the intersect of ray created by a position and the shape
        """

    @pet.provides
    def geodetic(self, points):
        """
        Function to convert Cartesian coordinates to geodetic
        """

    @pet.provides
    def cartesian(self, points):
        """
        Function to convert geodetic coordinates to Cartesian
        """

# end of file
