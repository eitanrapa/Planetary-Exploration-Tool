# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# california institute of technology
# (c) 2023-2024 all rights reserved

import pet


class Orbiter(pet.protocol, family="pet.campaigns"):
    """
    The abstract specification for orbiting campaigns
    """

    @pet.provides
    def get_states(self, target_body_id, times, reference_body):
        """
        Get the states of an instrument
        """

    @pet.provides
    def plot_orbit(self, visualization, planet, start_time, end_time, temporal_resolution, return_fig):
        """
        Plot the orbit of the instrument
        """


# end of file
