#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet


class Instrument(pet.protocol, family="pet.instruments"):
    """
    The abstract specification for radar instruments
    """

    @pet.provides
    def convert_time(self, time):
        """

        """

    @pet.provides
    def get_states(self, target_body_id, times, reference_body):
        """

        """

    @pet.provides
    def plot_orbit(self, target_body_id, start_time, end_time, reference_body):
        """

        """

# end of file
