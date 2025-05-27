#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import pet
import cspyce as spice
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


class NightingaleFoundational(pet.component, family="pet.campaigns.orbiter.nightingaleFoundational",
                      implements=pet.protocols.campaigns.orbiter):
    """
    Defines the Nightingale Foundational campaign
    """

    body_id = pet.properties.int()
    body_id.doc = "spk file body id"

    start_time = pet.properties.str()
    start_time.doc = "start time of coverage"

    planet = pet.protocols.planet()
    planet.doc = "target planet"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Get the start time in ephemeris time
        self.start_time_et = self.get_et_start_time()

        return

    def _get_states(self, times):
        """
         Get the states of an instrument
         :param times: The times to get the position and velocity for
         :return: List of positions, velocities of the instrument
         """

        # Get the states using SPICE toolkit
        states = spice.spkez_vector(targ=int(self.planet.body_id), et=times, ref=self.planet.reference_id,
                                    abcorr="None", obs=self.body_id)

        # Separate positions and velocities, convert to meters
        positions = np.asanyarray(states[0])[:, :3] * 1e3
        velocities = np.asanyarray(states[0])[:, 3:6] * 1e3

        # Return positions and velocities
        return positions, velocities

    @pet.export
    def get_states(self, times):
        """
        Helper function for _get_states
        """

        # Make sure times is a numpy array
        times = np.asarray(times)

        # Make sure dimensions is 1
        single_return = False
        if times.ndim == 0:
            single_return = True
            times = np.asarray([times])

        # Get positions and velocities
        positions, velocities = self._get_states(times=times)

        # Strip array if ndim is 0
        if single_return:
            return positions[0], velocities[0]
        else:
            return positions, velocities

    def get_et_start_time(self):
        """
        Get the Nightingale first five tracks split at a given latitude cutoff
        :param latitude_cutoff: Latitude to split tracks at
        :return: Times at which the cutoffs occur
        """

        # Create a time conversion instance
        time_conversion = pet.spiceTools.timeConversion()

        # Generate the ephemeris times to look through for 1 second of temporal spacing
        et = time_conversion.convert_ets(utcs=self.start_time)

        return et

    @pet.export
    def plot_orbit(self, projection, start_time, end_time, temporal_resolution=10, fig=None, globe=None, ax=None,
                   return_fig=False):
        """
        Plot the orbit of the instrument
        :param projection: Cartopy projection to use
        :param start_time: Start time of orbit plotting
        :param end_time: End time of orbit plotting
        :param temporal_resolution: Time interval for plot points
        :param fig: matplotlib figure
        :param globe: cartopy globe
        :param ax: matplotlib ax
        :param return_fig: Whether to return the fig, ax, globe objects
        :return: fig, ax, globe if return_fig is True
        """

        # Create a time conversion instance
        time_conversion = pet.spiceTools.timeConversion()

        # Create list of times to observe at
        times = np.arange(time_conversion.convert_ets(utcs=start_time),
                          time_conversion.convert_ets(utcs=end_time) + temporal_resolution, temporal_resolution)

        # Get the positions and velocities
        positions, velocities = self.get_states(times)

        # Get planet axes
        a, b, c = self.planet.get_axes()

        # Create a coordinate conversion object
        coordinate_conversions = pet.conversions.coordinateConversions(name="conversions", a=a, b=b, c=c)

        # Get geodetic coordinates
        geodetic_coordinates = coordinate_conversions.geodetic(cartesian_coordinates=positions)[:, :3]

        if fig is None:
            # Get the projection
            fig, ax, globe = projection.proj(planet=self.planet)

        # Get the coordinates
        coordinates = [(long, lat, height) for lat, long, height in geodetic_coordinates]
        longitudes, latitudes, heights = zip(*coordinates)

        # Plot points on the map
        ax.scatter(longitudes, latitudes, transform=ccrs.PlateCarree(globe=globe),
                   color='black', marker='o', s=0.1, alpha=0.25)

        # Return fig, ax, globe if necessary
        if return_fig:
            return fig, ax, globe

        # Add labels and legend
        plt.title('Orbit', pad=20)

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/nightingale_orbit_' + str(start_time) + '_' +
                    str(end_time) + '.png', format='png', dpi=500)

# end of file
