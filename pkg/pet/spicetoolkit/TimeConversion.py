#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import cspyce as spice
import numpy as np


class TimeConversion(pet.component):
    """

    """

    def _convert_times(self, times):
        """
        Converts a time string to Ephemeris Time using a loaded leap seconds file
        :param times: String to be converted
        :return: Ephemeris time float
        """

        # Convert time from string to ephemeris time
        ets = [spice.str2et(str=time) for time in times]

        # Return ephemeris time
        return np.asanyarray(ets)

    def convert_times(self, times):
        """
        Helper function for _convert_times
        """

        # Make sure times is a numpy array
        times = np.asanyarray(times)

        # Make sure dimensions is 1
        if times.ndim == 0:
            times = np.asanyarray([times])

        # Call _Cartesian function for results
        return self._convert_times(times=times)

    def _convert_utcs(self, ets):
        """
        Converts ETs to time strings using the loaded leap seconds file
        :param ets: Ephemeris times to be converted
        :return: Time string
        """

        # Convert time to UTC
        utcs = [spice.et2utc(et=et, format="c", prec=14) for et in ets]

        # Return the strings
        return utcs

    def convert_utcs(self, ets):
        """
        Helper function for _convert_utcs
        """

        # Make sure times is a numpy array
        ets = np.asanyarray(ets)

        # Make sure dimensions is 1
        if ets.ndim == 0:
            ets = np.asanyarray([ets])

        # Call _Cartesian function for results
        return self._convert_utcs(ets=ets)

# end of file
