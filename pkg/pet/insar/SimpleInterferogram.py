#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved
import pet


class SimpleInterferogram(pet.component):
    """

    """
    acquisition_strategy = pet.properties.dict()
    acquisition_strategy.doc = ""

    pairing_strategy = pet.properties.dict()
    pairing_strategy.doc = ""

    swath_times = pet.properties.list()
    swath_times.doc = ""

    def get_flattened_angle(self):
        """

        """

    def get_flattened_phase(self, ):
        """

        """

    def get_interferogram(self, track_number, pairing_number):
        """

        """

    def plot_interferogram(self):
        """

        """

# end of file
