#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import pet


class InSAR(pet.protocol, family="pet.instruments"):
    """
    The abstract specification for InSAR instruments
    """

    @pet.provides
    def get_instrument_noise(self, planet, baseline, satellite_velocity, look_angles, incidence_angles, distances):
        """
        Get the instrument noise for InSAR
        """

# end of file
