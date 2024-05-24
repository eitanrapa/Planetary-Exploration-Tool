#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet


class GroundTarget(pet.component):
    """
    Class that represents a target point on the ground
    """

    x = pet.properties.float()
    x.doc = "x coordinate"

    y = pet.properties.float()
    y.doc = "y coordinate"

    z = pet.properties.float()
    z.doc = "z coordinate"

    def get_position(self):
        """
        Return the x, y, z positions of the GroundTarget
        :return: x, y, z coordinates
        """

        # Return x, y, z
        return self.x, self.y, self.z

    def attach_geodetic(self):
        """
        Attach geodetic coordinates to the GroundTarget
        """

# end of file
