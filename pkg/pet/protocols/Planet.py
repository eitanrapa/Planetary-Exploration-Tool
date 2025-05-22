# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# california institute of technology
# (c) 2023-2025 all rights reserved

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
    def get_surface_intersects_local_angles(self, satellite_position, raydirs):
        """
        Find the intersects of a set of vectors with the planet DSK
        """

    @pet.provides
    def get_closest_point_to_surface(self, points):
        """
        Get the closest intersect of a point with the planet and the distance to the surface
        """

    @pet.provides
    def get_sub_obs_points(self, times, campaign):
        """
        Get the nadir sub-observation points of the instrument with the DSK
        """

    @pet.provides
    def visualize_topography(self, projection, fig, globe, ax, return_fig):
        """
        Creates a visualization of the surface of Enceladus
        """

# end of file
