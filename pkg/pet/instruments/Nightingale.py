#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import spiceypy as spice


class Nightingale(pet.component, family="pet.instruments.nightingale", implements=pet.protocols.instrument):
    """
    Defines the Nightingale mission instrument
    """

    body_id = pet.properties.int()
    body_id.doc = "spk file body id"

    start_look_angle = pet.properties.float()
    end_look_angle = pet.properties.float()

    @pet.export
    def convert_time(self, time):
        """

        """

        # Convert time from string to Ephemeris Time
        et = spice.str2et(time=time)
        return et

    @pet.export
    def get_state(self, target_body_id, time, reference_body):
        if isinstance(time, str):
            et = self.convert_time(time)
        else:
            et = time
        current = spice.spkez(targ=int(target_body_id), et=et, ref=reference_body, abcorr="None", obs=self.body_id)[0]
        position = current[:3]  # in km
        velocity = current[3:6]  # in km

        return position, velocity

    @pet.export
    def plot_orbit(self, target_body_id, start_time, end_time, reference_body):

        return

# end of file
