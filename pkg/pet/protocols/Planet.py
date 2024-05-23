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
        Return the axes of the planet
        """

    @pet.provides
    def get_surface_intersects(self, vectors):
        """
        Find the intersects of a set of vectors with the planet DSK
        """

    @pet.provides
    def get_sub_obs_points(self, times, instrument_id):
        """
        Get the nadir sub-observation points of the instrument with the DSK
        """

    @pet.provides
    def visualize_topography(self):
        """
        Creates a visualization of the surface of Enceladus
        """

# end of file
