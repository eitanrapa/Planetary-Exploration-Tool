#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet


class Projection(pet.protocol, family="pet.projections"):
    """
    The abstract specification for projections
    """

    @pet.provides
    def proj(self, geodetic_coordinates):
        """
        Provides the cartopy proj
        """

# end of file
