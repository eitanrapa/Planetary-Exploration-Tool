#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import cspyce as spice


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
        et = spice.str2et(str=time)
        return et

    @pet.export
    def get_states(self, target_body_id, times, reference_body):
        if isinstance(times[0], str):
            ets = [self.convert_time(time=time) for time in times]
        else:
            ets = times
        states = spice.spkez_vector(targ=int(target_body_id), et=ets, ref=reference_body, abcorr="None",
                                    obs=self.body_id)

        positions = states[0][:, :3]  # in km
        velocities = states[0][:, 3:6]  # in km
        return positions, velocities

    @pet.export
    def plot_orbit(self, target_body_id, start_time, end_time, reference_body):

        return

# end of file
