# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# california institute of technology
# (c) 2023-2024 all rights reserved

import pet


class Planet(pet.protocol, family="pet.planets"):
    """
    The abstract specification for planetary bodies
    """

    @pet.provides
    def get_axes(self):
        """

        """

    @pet.provides
    def get_surface_intersects(self, vectors):
        """

        """

    @pet.provides
    def get_sub_obs_points(self, times, instrument_id):
        """

        """

    @pet.provides
    def visualize_topography(self):
        """
        Visualizes the topography that defines the surface of a planet
        """

# end of file
