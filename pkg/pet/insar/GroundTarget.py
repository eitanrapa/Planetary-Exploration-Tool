#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet


class GroundTarget:
    """
    Class that represents a target point on the ground
    """

    def __init__(self, x, y, z, seen):
        self.x = x
        self.y = y
        self.z = z
        self.seen = seen

    def get_position(self):
        """
        Return the x, y, z positions of the GroundTarget
        :return: x, y, z coordinates
        """

        # Return x, y, z
        return self.x, self.y, self.z

# end of file
