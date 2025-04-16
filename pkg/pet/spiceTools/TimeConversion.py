#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import pet
import cspyce as spice
import numpy as np


class TimeConversion(pet.component):
    """
    Class that encapsulates time conversions
    """

    def _convert_ets(self, utcs):
        """
        Converts a time string to Ephemeris Time using a loaded leap seconds file
        :param utcs: String to be converted
        :return: Ephemeris time float
        """

        # Convert time from string to ephemeris time
        ets = [spice.str2et(str=utc) for utc in utcs]

        # Return ephemeris time
        return np.asarray(ets)

    def convert_ets(self, utcs):
        """
        Helper function for _convert_ets
        """

        # Make sure times is a numpy array
        utcs = np.asarray(utcs)

        # Make sure dimensions is 1
        single_return = False
        if utcs.ndim == 0:
            single_return = True
            utcs = np.asarray([utcs])

        # Call _Cartesian function for results
        ets = self._convert_ets(utcs=utcs)

        # Strip array if ndim is 0
        if single_return:
            return ets[0]
        else:
            return ets

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
        ets = np.asarray(ets)

        # Make sure dimensions is 1
        single_return = False
        if ets.ndim == 0:
            single_return = True
            ets = np.asarray([ets])

        # Call _Cartesian function for results
        utcs = self._convert_utcs(ets=ets)

        # Strip array if ndim is 1
        if single_return:
            return utcs[0]
        else:
            return utcs

# end of file
