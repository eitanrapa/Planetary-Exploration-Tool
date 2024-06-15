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


class Nightingale(pet.component, family="pet.instruments.nightingale", implements=pet.protocols.instrument):
    """
    Defines the Nightingale mission instrument
    """

    body_id = pet.properties.int()
    body_id.doc = "spk file body id"

    start_look_angle = pet.properties.float()
    start_look_angle.doc = "Starting look angle"

    end_look_angle = pet.properties.float()
    end_look_angle.doc = "Ending look angle"

    orbit_cycle = pet.properties.float()
    orbit_cycle.doc = "Time to make one orbit"

    wavelength = pet.properties.float()
    wavelength.doc = "Radar wavelength"

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

    @pet.export
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

    @pet.export
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

    def _get_states(self, planet, times):
        """
         Get the states of an instrument
         :param planet: Target planet
         :param times: The times to get the position and velocity for
         :return: List of positions, velocities of the instrument
         """

        # Get the states using SPICE toolkit
        states = spice.spkez_vector(targ=int(planet.body_id), et=times, ref=planet.reference_id, abcorr="None",
                                    obs=self.body_id)

        # Separate positions and velocities, convert to meters
        positions = np.asanyarray(states[0])[:, :3] * 1e3
        velocities = np.asanyarray(states[0])[:, 3:6] * 1e3

        # Return positions and velocities
        return positions, velocities

    @pet.export
    def get_states(self, planet, times):
        """
        Helper function for _get_states
        """

        # Make sure times is a numpy array
        times = np.asanyarray(times)

        # Make sure dimensions is 1
        if times.ndim == 0:
            times = np.asanyarray([times])

        # Return positions and velocities
        return self._get_states(planet=planet, times=times)

    def get_five_tracks(self, planet, latitude_cutoff=0, start_time="2046 DEC 20 15:10:40.134"):
        """
        Get the Nightingale first five tracks split at a given latitude cutoff
        :param planet: Target planet
        :param latitude_cutoff: Latitude to split tracks at
        :param start_time: Start time of the tracks
        :return: Times at which the cutoffs occur
        """

        # Generate the ephemeris times to look through for 1 second of temporal spacing
        ets = [self.convert_times(times=start_time)[0] + i for i in range(162000)]  # 45 hour search

        # Get the positions from the states
        positions, velocities = self.get_states(planet=planet, times=ets)

        # Iterate through the positions
        times = [start_time]

        # Get planet axes
        a, b, c = planet.get_axes()

        # Create a coordinate conversion object
        convert = pet.ext.conversions(name="conversions", a=a, b=b, c=c)

        # Get the geodetic coordinates
        geodetic_coordinates = np.asarray(convert.geodetic(cartesian_coordinates=positions))

        for i in range(len(geodetic_coordinates) - 1):
            # If the latitude cutoff is found, attach the time of the cutoff
            if (geodetic_coordinates[i, 0] < latitude_cutoff) and (geodetic_coordinates[i + 1, 0] > latitude_cutoff):
                times.append(self.convert_utcs(ets[i])[0])

        # Return the times
        return times

    @pet.export
    def plot_orbit(self, visualization, projection, planet, start_time, end_time, time_interval=10, return_fig=False):
        """
        Plot the orbit of the instrument
        :param visualization: Visualization tool to use
        :param projection: Cartopy projection to use
        :param planet: The planet to visualize as well
        :param start_time: Start time of orbit plotting
        :param end_time: End time of orbit plotting
        :param time_interval: Time interval for plot points
        :param return_fig: Whether to return the fig, ax, globe objects
        """

        # Create list of times to observe at
        times = np.arange(self.convert_times(times=start_time), self.convert_times(times=end_time) + time_interval,
                          time_interval)

        # Get the positions and velocities
        positions, velocities = self.get_states(planet, times)

        # Get planet axes
        a, b, c = planet.get_axes()

        # Create a coordinate conversion object
        convert = pet.ext.conversions(name="conversions", a=a, b=b, c=c)

        # Get geodetic coordinates
        geodetic_coordinates = convert.geodetic(cartesian_coordinates=positions)

        # Get fig, ax, globe from topography
        fig, ax, globe = planet.visualize_topography(visualization=visualization,
                                                     projection=projection, return_fig=True)

        # Use visualization tool to plot
        fig, ax, globe = visualization.scatter_plot(fig=fig, ax=ax, globe=globe,
                                                    geodetic_coordinates=geodetic_coordinates[:, :3], color='black')

        # Return fig, ax, globe if necessary
        if return_fig:
            return fig, ax, globe

        # Add labels and legend
        ax.set_title('Orbit')

        # Save the plot
        plt.savefig(fname=visualization.folder_path + '/orbit.png', format='png', dpi=500)

# end of file
