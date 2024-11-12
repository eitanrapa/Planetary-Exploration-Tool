#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import numpy as np


class Surface(pet.component):
    """

    """

    def sigma_wye(self, incidence_angle):
        """

        """

        # Convert to radians
        incidence_angle_radians = np.radians(incidence_angle)

        # Modify the amplitude
        modified_amplitude = 3.52 * np.power(np.cos(incidence_angle_radians), 1.23)
        return modified_amplitude

# end of file
