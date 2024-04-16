#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import spiceypy

class PolarOrbit(pet.component, family="pet.orbits.fivetoone", implements=pet.protocols.orbit):
    """
    Defines an orbit which is synchronous to a pole
    """

    orbit = pet.properties.path()
    orbit.doc = "path to spk file"


# end of file
