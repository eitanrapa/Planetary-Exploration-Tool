#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import pet


class BiaxialProjection(pet.protocol, family="pet.projections"):
    """
    The abstract specification for biaxial projections
    """

    @pet.provides
    def proj(self, planet, west_extent, east_extent, south_extent, north_extent):
        """
        Provides the cartopy proj
        """

# end of file
