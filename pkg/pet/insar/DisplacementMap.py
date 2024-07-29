#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import h5py
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d
import numpy as np


class DisplacementMap(pet.component):
    """
    Class that represents an instance of displacement values created by a Finite Element Model.
    """

    displacement_data_path = pet.properties.str()
    displacement_data_path.doc = "path to hdf5 file containing surface displacement"

    def __init__(self, name, locator, implicit, planet):
        super().__init__(name, locator, implicit)
        self.planet = planet

    def read_displacements(self):
        """
        Returns the latitudes, longitudes, times, and cycle time of the data cubes representing the surface displacement
        over a tidal cycle, plus the cubes themselves.
        :return: Latitudes, longitudes, times, cycle time, corresponding to the data cubes representing the east, north,
        and up vectors at each point
        """

        # Read the file with h5py
        displacement_file = h5py.File(name=self.displacement_data_path, mode='r')

        # Get the cubes in the displacement_values folder
        east_cube = displacement_file["displacement_values/east_cube"][:]
        north_cube = displacement_file["displacement_values/north_cube"][:]
        up_cube = displacement_file["displacement_values/up_cube"][:]

        # Get the latitudes, longitudes corresponding to the cubes
        latitudes = displacement_file["latitudes"][:]
        longitudes = displacement_file["longitudes"][:]

        # Get the cycle time and the fractional times
        cycle_time = displacement_file["times"].attrs["Cycle Time"]
        times = displacement_file["times"][:]

        # Close the file
        displacement_file.close()

        # Return the retrieved values
        return latitudes, longitudes, times, cycle_time, east_cube, north_cube, up_cube

    def attach(self, swath1, swath2, time_displacement=0, use_mid_point=False):
        """
        Attach to each GroundTarget of a GroundSwath object the displacement at that point at the time it was
        observed using a regular grid interpolation.
        :param swath1: GroundSwath object containing GroundTarget objects
        :param swath2: second GroundSwath object containing GroundTarget objects
        :param time_displacement: How much time to advance the model by [s]
        :param use_mid_point: Use the middle point time of the swath as the measuring point of the displacement field
        :return: Nothing returned
        """

        # Get list of times to observe at
        time_space1 = swath1.time_space

        time_space2 = swath2.time_space

        # Read the necessary data
        latitudes, longitudes, times, cycle_time, east_cube, north_cube, up_cube = self.read_displacements()

        # Create interpolating functions for each dimension
        east_interp = RegularGridInterpolator(points=(longitudes, latitudes, times),
                                              values=east_cube, method='slinear', bounds_error=False,
                                              fill_value=None)
        north_interp = RegularGridInterpolator(points=(longitudes, latitudes, times),
                                               values=north_cube, method='slinear', bounds_error=False,
                                               fill_value=None)
        up_interp = RegularGridInterpolator(points=(longitudes, latitudes, times),
                                            values=up_cube, method='slinear', bounds_error=False, fill_value=None)

        # Get planet axes
        a, b, c = self.planet.get_axes()

        # Create a coordinate conversion object
        convert = pet.ext.conversions(name="conversions", a=a, b=b, c=c)

        # Populate the GroundTarget displacement values at each time
        for i in range(len(swath1.swath_beams)):
            for point in swath1.swath_beams[i]:

                # Pack the coordinates
                cartesian_coordinates = point.x, point.y, point.z

                # Get the corresponding latitude and longitude with a triaxial ellipsoid conversion
                lat, long = convert.geodetic(cartesian_coordinates=cartesian_coordinates)[0, :2]

                # Use mid_point if True
                if use_mid_point:
                    mid_index = len(time_space1) // 2

                    # Calculate the fractional time corresponding to the cycle
                    modified_time = ((time_space1[mid_index]) % cycle_time) / cycle_time

                else:
                    # Calculate the fractional time corresponding to the cycle
                    modified_time = ((time_space1[i]) % cycle_time) / cycle_time

                # Attach the displacement to the point
                point.displacement = [east_interp((long, lat, modified_time)).item(),
                                      north_interp((long, lat, modified_time)).item(),
                                      up_interp((long, lat, modified_time)).item(),
                                      lat, long]

        # Populate the GroundTarget displacement values at each time
        for i in range(len(swath2.swath_beams)):
            for point in swath2.swath_beams[i]:

                # Pack the coordinates
                cartesian_coordinates = point.x, point.y, point.z

                # Get the corresponding latitude and longitude with a triaxial ellipsoid conversion
                lat, long = convert.geodetic(cartesian_coordinates=cartesian_coordinates)[0, :2]

                # Use mid_point if True
                if use_mid_point:
                    mid_index = len(time_space2) // 2

                    # Calculate the fractional time corresponding to the cycle
                    modified_time = ((time_space2[mid_index] + time_displacement) % cycle_time) / cycle_time

                else:
                    # Calculate the fractional time corresponding to the cycle
                    modified_time = ((time_space2[i] + time_displacement) % cycle_time) / cycle_time

                # Attach the displacement to the point
                point.displacement = [east_interp((long, lat, modified_time)).item(),
                                      north_interp((long, lat, modified_time)).item(),
                                      up_interp((long, lat, modified_time)).item(),
                                      lat, long]

    def visualize(self, time_point, projection, direction):
        """
        Visualize the displacement map at a specific time
        :param time_point: Point in cycle at which to view displacements [s]
        :param projection: Cartopy projection
        :param direction: Vector direction to view, east, north, or up
        """

        # Get the projection
        fig, ax, globe = projection.proj(planet=self.planet)

        # Read the necessary data
        latitudes, longitudes, times, cycle_time, east_cube, north_cube, up_cube = self.read_displacements()

        modified_time = (time_point % cycle_time) / cycle_time

        # Visualize correct direction
        if direction == "east":

            # Initialize an array to store the interpolated values
            interpolated_values = np.zeros((len(east_cube[:, 0, 0]), len(east_cube[0, :, 0])))

            # Loop through each (i, j) point in the 2D grid
            for i in range(len(east_cube[:, 0, 0])):
                for j in range(len(east_cube[0, :, 0])):

                    # Extract the time series for the (i, j) point
                    time_series = east_cube[i, j, :]

                    # Create the interpolator function
                    interpolator = interp1d(times, time_series, kind='linear', bounds_error=True)

                    # Interpolate the value at time t
                    interpolated_values[i, j] = interpolator(modified_time)

            im = ax.imshow(interpolated_values.T, extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(globe=globe))

        elif direction == "north":

            # Initialize an array to store the interpolated values
            interpolated_values = np.zeros((len(north_cube[:, 0, 0]), len(north_cube[0, :, 0])))

            # Loop through each (i, j) point in the 2D grid
            for i in range(len(north_cube[:, 0, 0])):
                for j in range(len(north_cube[0, :, 0])):

                    # Extract the time series for the (i, j) point
                    time_series = north_cube[i, j, :]

                    # Create the interpolator function
                    interpolator = interp1d(times, time_series, kind='linear', bounds_error=True)

                    # Interpolate the value at time t
                    interpolated_values[i, j] = interpolator(modified_time)

            im = ax.imshow(interpolated_values.T, extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(globe=globe))

        elif direction == "up":

            # Initialize an array to store the interpolated values
            interpolated_values = np.zeros((len(up_cube[:, 0, 0]), len(up_cube[0, :, 0])))

            # Loop through each (i, j) point in the 2D grid
            for i in range(len(up_cube[:, 0, 0])):
                for j in range(len(up_cube[0, :, 0])):

                    # Extract the time series for the (i, j) point
                    time_series = up_cube[i, j, :]

                    # Create the interpolator function
                    interpolator = interp1d(times, time_series, kind='linear', bounds_error=True)

                    # Interpolate the value at time t
                    interpolated_values[i, j] = interpolator(modified_time)

            im = ax.imshow(interpolated_values.T, extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(globe=globe))

        else:
            raise Exception("direction must be east, north, or up")

        # Add a colorbar
        plt.colorbar(im, ax=ax, orientation='vertical', shrink=0.25, label="[m]")

        # Add labels and legend
        ax.set_title('Displacements ' + direction)

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/displacements_' + direction + '.png', format='png', dpi=500)

# end of file
