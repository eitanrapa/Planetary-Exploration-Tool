# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import pyre
from pyre.units.SI import meter, second
from pyre.units.angle import degree


class Enceladus(pet.component, family="pet.planets.enceladus", implements=pet.protocols.planet):
    """
    An encapsulation of the data and parameters used to represent the planetary body of Enceladus
    """

    # geometry
    major_equatorial_semiaxis = 256.6 * 1e3 * meter
    minor_equatorial_semiaxis = 251.4 * 1e3 * meter
    polar_semiaxis = 248.3 * 1e3 * meter
    orbital_time = 118386.8352 * second
    inclination = 0.009 * degree

    @pet.export
    def topography(self):
        """
        Provides a list of control points that defines the surface of Enceladus
        """
        return []

# end of file
