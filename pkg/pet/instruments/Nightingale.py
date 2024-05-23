#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import cspyce as spice
from ..ext import conversions


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
        Converts a time string to Ephemeris Time using a loaded leap seconds file
        :param time: String to be converted
        :return: Ephemeris time float
        """

        # Convert time from string to ephemeris time
        et = spice.str2et(str=time)

        # Return ephemeris time
        return et

    @pet.export
    def get_states(self, target_body_id, times, reference_body):
        """
        Get the states of an instrument
        :param target_body_id: ID of the target body
        :param times: The times to get the position and velocity for
        :param reference_body: IAU reference frame
        :return: List of positions, velocities of the instrument
        """

        # If the time is a string, convert first. If not, use as ephemeris time
        if isinstance(times[0], str):
            ets = [self.convert_time(time=time) for time in times]
        else:
            ets = times

        # Get the states using SPICE toolkit
        states = spice.spkez_vector(targ=int(target_body_id), et=ets, ref=reference_body, abcorr="None",
                                    obs=self.body_id)

        # Separate positions and velocities
        positions = states[0][:, :3]  # in km
        velocities = states[0][:, 3:6]  # in km

        # Return positions and velocities
        return positions, velocities

    def get_five_tracks(self, planet, latitude_cutoff=0, start_time="2046 DEC 20 15:10:40.134"):
        """
        Get the Nightingale first five tracks split at a given latitude cutoff
        :param planet: Target planet
        :param latitude_cutoff: Latitude to split tracks at
        :param start_time: Start time of the tracks
        :return: Times at which the cutoffs occur
        """

        # Generate the ephemeris times to look through for 1 second of temporal spacing
        ets = [self.convert_time(time=start_time) + i for i in range(360000)]

        # Get the states using the SPICE toolkit
        states = spice.spkez_vector(targ=int(planet.body_id), et=ets, ref=planet.reference_id, abcorr="None",
                                    obs=self.body_id)

        # Get the positions from the states
        positions = states[0][:, :3]  # in km

        # Get the planet axes
        planet_axes = planet.get_axes()

        # Iterate through the positions
        times = []
        for i in range(len(positions) - 1):
            # If the latitude cutoff is found, attach the time of the cutoff
            if ((conversions.geodetic(planet_axes=planet_axes,
                                      cartesian_coordinates=positions[i])[0] < latitude_cutoff) and
                    (conversions.geodetic(planet_axes=planet_axes, cartesian_coordinates=positions[i + 1])[0] >
                     latitude_cutoff)):
                times.append([ets[i]])

        # Return the times
        return times

    @pet.export
    def plot_orbit(self, target_body_id, start_time, end_time, reference_body):
        """
        Plot the orbit of the instrument
        """
        return

# end of file
