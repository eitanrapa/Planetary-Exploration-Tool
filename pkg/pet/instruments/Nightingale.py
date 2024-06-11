#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import cspyce as spice
import numpy as np
import os
import matplotlib.pyplot as plt
from ..ext import conversions
from ..projections import plottingTools


class Nightingale(pet.component, family="pet.instruments.nightingale", implements=pet.protocols.instrument):
    """
    Defines the Nightingale mission instrument
    """

    body_id = pet.properties.int()
    body_id.doc = "spk file body id"

    start_look_angle = pet.properties.float()
    start_look_angle.doc = ""

    end_look_angle = pet.properties.float()
    end_look_angle.doc = ""

    orbit_cycle = pet.properties.float()
    orbit_cycle.doc = ""

    wavelength = pet.properties.float()
    wavelength.doc = ""

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

    def convert_utc(self, et):
        """

        """

        utc = spice.et2utc(et=et, format="c", prec=14)

        return utc

    @pet.export
    def get_states(self, planet, times):
        """
        Get the states of an instrument
        :param planet: Target planet
        :param times: The times to get the position and velocity for
        :return: List of positions, velocities of the instrument
        """

        times = np.array(times, ndmin=1)
        if times.ndim == 1:

            # If the time is a string, convert first. If not, use as ephemeris time
            if isinstance(times, str):
                ets = self.convert_time(time=times)
            else:
                ets = times
        else:

            if isinstance(times[0], str):
                ets = self.convert_time(time=times)
            else:
                ets = times

        # Get the states using SPICE toolkit
        states = spice.spkez_vector(targ=int(planet.body_id), et=ets, ref=planet.reference_id, abcorr="None",
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
        ets = [self.convert_time(time=start_time) + i for i in range(162000)]  # 45 hour search

        # Get the positions from the states
        positions, velocities = self.get_states(planet=planet, times=ets)

        # Iterate through the positions
        times = [start_time]
        geodetic_coordinates = np.asarray(conversions.geodetic(planet=planet, cartesian_coordinates=positions))

        for i in range(len(geodetic_coordinates) - 1):
            # If the latitude cutoff is found, attach the time of the cutoff
            if (geodetic_coordinates[i, 0] < latitude_cutoff) and (geodetic_coordinates[i + 1, 0] > latitude_cutoff):
                times.append(self.convert_utc(ets[i]))

        # Return the times
        return times

    @pet.export
    def plot_orbit(self, planet, projection, start_time, end_time, north_extent=90, south_extent=-90, east_extent=180,
                   west_extent=-180, time_interval=10, return_fig=False):
        """
        Plot the orbit of the instrument
        """

        # Create list of times to observe at
        times = np.arange(self.convert_time(time=start_time), self.convert_time(time=end_time) + time_interval,
                          time_interval)

        positions, velocities = self.get_states(planet, times)

        geodetic_coordinates = np.asarray(
            conversions.geodetic(planet=planet, cartesian_coordinates=positions))

        fig, ax, globe = planet.visualize_topography(projection=projection, north_extent=north_extent,
                                                     south_extent=south_extent,
                                                     east_extent=east_extent,
                                                     west_extent=west_extent, return_fig=True)

        fig, ax = plottingTools.scatter_plot(fig=fig, ax=ax, globe=globe,
                                             geodetic_coordinates=geodetic_coordinates[:, :3],
                                             north_extent=north_extent, south_extent=south_extent,
                                             east_extent=east_extent, west_extent=west_extent)

        if return_fig:
            return fig, ax, globe

        # Add labels and legend
        ax.set_title('Orbit')

        # Get the current working directory
        current_dir = os.getcwd()

        # Construct the path three directories up
        path_three_dirs_up = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir, os.pardir))

        # Define the target directory within the three-up directory
        target_dir = os.path.join(path_three_dirs_up, 'figs')

        # Define the full path to save the plot
        full_path = os.path.join(target_dir, 'orbit.png')

        # Show the plot
        plt.savefig(fname=full_path, format='png', dpi=500)

# end of file
