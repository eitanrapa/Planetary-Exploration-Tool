# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# california institute of technology
# (c) 2023-2024 all rights reserved
#

import pet


class Planet(pet.protocol, family="pet.planets"):
    """
    The abstract specification for planetary bodies
    """

    @pet.provides
    def topography(self):
        """
        Generate an array of points that defines surface of planet
        """

    @pet.provides
    def visualize_topography(self):
        """
        Visualizes the topography that defines the surface of a planet
        """

# end of file
