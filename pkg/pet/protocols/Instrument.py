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
    def convert_times(self, times):
        """
        Converts time strings to Ephemeris Time using a loaded leap seconds file
        """

    @pet.provides
    def convert_utcs(self, time):
        """
        Converts ETs to time strings using the loaded leap seconds file
        """

    @pet.provides
    def get_states(self, target_body_id, times, reference_body):
        """
        Get the states of an instrument
        """

    @pet.provides
    def plot_orbit(self, visualization, planet, start_time, end_time, time_interval, return_fig):
        """
        Plot the orbit of the instrument
        """

# end of file
