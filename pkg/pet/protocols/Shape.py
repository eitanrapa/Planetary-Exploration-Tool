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
    def intersect(self, position):
        """
        A function to calculate the intersect of ray created by a position and the shape
        """

    @pet.provides
    def cartesian_to_geodetic(self, x, y, z):
        """
        Function to convert between cartesian and geodetic coordinates
        """

    @pet.provides
    def geodetic_to_cartesian(self):
        """
        Function to convert between geodetic and cartesian coordinates
        """

# end of file
