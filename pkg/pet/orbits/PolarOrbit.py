#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet


class PolarOrbit(pet.component, family="pet.orbits.fivetoone", implements=pet.protocols.orbit):
    """
    Defines an orbit which is synchronous to a pole
    """

# end of file
