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
    def intersect(self):
        """

        """
        return []

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
