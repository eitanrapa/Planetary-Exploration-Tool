#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import cspyce as spice
import numpy as np


class CreateOrbit(pet.component):
    """

    """

    def make_SPK_from_states(self):

        handle = spice.spkopn(fname="orbit_file.bsp", ifname="orbit", ncomch=10)

    def make_circular_orbit_SPK(self, altitude):

        handle = spice.spkopn(fname="orbit_file.bsp", ifname="orbit", ncomch=10)


# end of file
