#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
from pyre.units.SI import meter, second
from pyre.units.angle import degree


class Nightingale(pet.component, family="pet.instruments.nightingale", implements=pet.protocols.instrument):
    """
    Defines the Nightingale mission instrument
    """

# end of file
