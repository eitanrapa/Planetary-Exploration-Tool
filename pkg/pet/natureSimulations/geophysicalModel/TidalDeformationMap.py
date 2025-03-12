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
from scipy.optimize import curve_fit
import numpy as np
from tqdm import tqdm


class TidalDeformationMap(pet.component, family="pet.natureSimulations.geophysicalModel.tidalDeformationMap",
                          implements=pet.protocols.natureSimulations.geophysicalModel):

    """
    Class that represents an instance of displacement values created by a Finite Element Model.
    """

    displacement_data_path = pet.properties.path()
    displacement_data_path.doc = "path to the displacement data file"

    planet = pet.protocols.planet()
    planet.doc = "target planet"

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
        colatitudes = displacement_file["colatitudes"][:]
        longitudes = displacement_file["longitudes"][:]

        # Get the cycle time and the fractional times
        cycle_time = displacement_file["times"].attrs["Cycle Time"]
        times = displacement_file["times"][:]

        # Close the file
        displacement_file.close()

        # Return the retrieved values
        return colatitudes, longitudes, times, cycle_time, east_cube, north_cube, up_cube

    def get_displacements(self, track):
        """
        Get surface displacements for the positions defined by a track
        :param track: object for which to find displacements for longitudes and latitudes
        :return: displacements in local cartesian vectors
        """

        # Read the necessary data
        colatitudes, longitudes, times, cycle_time, east_cube, north_cube, up_cube = self.read_displacements()

        # Convert colatitudes to latitudes
        latitudes = 90 - colatitudes

        vector_conversions = pet.conversions.vectorConversions()

        # Convert cubes of east, north, up vectors to u, v, w (local x, y, z coordinates)
        u_cube, v_cube, w_cube = vector_conversions.enu_to_cartesian_cubes(east_cube=east_cube, north_cube=north_cube,
                                                                           up_cube=up_cube,
                                                                           latitudes=latitudes, longitudes=longitudes)

        # Stack the u, v, w arrays along a new axis to create a 4D array
        total_cube = np.stack((u_cube, v_cube, w_cube), axis=-1)

        # Create interpolating function
        vector_interp = RegularGridInterpolator(points=(longitudes, latitudes, times),
                                                values=total_cube, method='linear', bounds_error=False,
                                                fill_value=None)

        # Access longitude and latitude coordinates
        longitudes = track.data["longitude"].values
        latitudes = track.data["latitude"].values
        time_space = track.data["time"].values

        # Convert longitudes to 0 - 360 format
        modified_longs = np.asarray(longitudes) % 360

        # Calculate the fractional time corresponding to the cycle
        modified_times = (np.array(time_space) % cycle_time) / cycle_time

        # Attach the displacement to the point
        u_displacements, v_displacements, w_displacements = vector_interp((modified_longs, latitudes, modified_times)).T

        return u_displacements, v_displacements, w_displacements

    def time_series(self, spatial_points, direction):
        """
        Get the time series of displacements for a set of spatial points
        :param spatial_points: List of spatial points in the form (latitude, longitude)
        :param direction: Vector direction to view, east, north, or up
        :return: Time series of displacements for each spatial point
        """

        # Read the necessary data
        colatitudes, longitudes, times, cycle_time, east_cube, north_cube, up_cube = self.read_displacements()

        # Convert colatitudes to latitudes
        latitudes = 90 - colatitudes

        # Define the sine function
        def sine_func(t, a, phi):
            return a * np.sin((2 * np.pi) / self.planet.tidal_cycle * t + phi)

        a_fits = []
        phi_fits = []

        if direction == "east":

            # Interpolate the east cube
            interpolator = RegularGridInterpolator(points=(longitudes, latitudes,
                                                           times * self.planet.tidal_cycle),
                                                   values=east_cube, method='linear', bounds_error=False,
                                                   fill_value=None)
        elif direction == "north":

            # Interpolate the north cube
            interpolator = RegularGridInterpolator(points=(longitudes, latitudes,
                                                           times * self.planet.tidal_cycle),
                                                   values=north_cube, method='linear', bounds_error=False,
                                                   fill_value=None)
        elif direction == "up":

            # Interpolate the up cube
            interpolator = RegularGridInterpolator(points=(longitudes, latitudes,
                                                           times * self.planet.tidal_cycle),
                                                   values=up_cube, method='linear', bounds_error=False,
                                                   fill_value=None)

        else:
            raise Exception("direction must be east, north, or up")

        for point in tqdm(spatial_points, "Interpolating points, fitting curves..."):

            # Create an array of time values from 0 to the tidal cycle
            time_values = np.linspace(0, self.planet.tidal_cycle, len(times))

            # Extract the time series for the spatial point, convert longitude to 0 - 360 forma
            time_series = interpolator((point[1] % 360, point[0], time_values))

            # Fit the sine function to the data
            popt, pcov = curve_fit(f=sine_func, xdata=time_values, ydata=time_series, p0=(1, 1))

            # Extract the fitted parameters
            a_fit, phi_fit = popt

            # If a_fit is negative, take the absolute value and add pi to phi_fit
            if a_fit < 0:
                a_fit = np.abs(a_fit)
                phi_fit += np.pi

            a_fits.append(a_fit)
            phi_fits.append(phi_fit)

        return a_fits, phi_fits

    def visualize(self, time_point, projection, direction, fig=None, globe=None, ax=None, return_fig=False):
        """
        Visualize the displacement map at a specific time
        :param time_point: Point in cycle at which to view displacements [s]
        :param projection: Cartopy projection
        :param direction: Vector direction to view, east, north, or up
        :param fig: matplotlib figure
        :param globe: cartopy globe
        :param ax: matplotlib ax
        :param return_fig: Whether to return fig, globe, ax
        :return: fig, globe, ax if return_fig is True
        """

        if fig is None:
            # Get the projection
            fig, ax, globe = projection.proj(planet=self.planet)

        # Read the necessary data
        colatitudes, longitudes, times, cycle_time, east_cube, north_cube, up_cube = self.read_displacements()

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
                    interpolator = interp1d(x=times, y=time_series, kind='linear', bounds_error=True)

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
                    interpolator = interp1d(x=times, y=time_series, kind='linear', bounds_error=True)

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
                    interpolator = interp1d(x=times, y=time_series, kind='linear', bounds_error=True)

                    # Interpolate the value at time t
                    interpolated_values[i, j] = interpolator(modified_time)

            im = ax.imshow(interpolated_values.T, extent=[-180, 180, -90, 90], transform=ccrs.PlateCarree(globe=globe))

        else:
            raise Exception("direction must be east, north, or up")

        # return fig, ax, globe if necessary
        if return_fig:
            return fig, ax, globe

        # Add a colorbar
        plt.colorbar(im, ax=ax, orientation='vertical', shrink=0.25, label="[m]")

        # Add labels and legend
        plt.title('Displacements ' + direction, pad=20)

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/displacements_' + str(self.pyre_name) + '_' +
                    direction + '_' + str(modified_time) + '.png', format='png', dpi=500)

# end of file
