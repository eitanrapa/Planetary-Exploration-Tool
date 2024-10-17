#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import h5py
import git
import cspyce as spice
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


class Interferogram(pet.component):
    """
    Class that creates a single interferogram between two points of a displacement map given an instrument orbit
    """

    def __init__(self, planet, instrument, deformation_map, swath1, swath2, time_space1,
                 time_space2, baseline, **kwargs):
        super().__init__(**kwargs)
        self.planet = planet
        self.instrument = instrument
        self.deformation_map = deformation_map
        self.swath1 = swath1
        self.swath2 = swath2
        self.time_space1 = time_space1
        self.time_space2 = time_space2
        self.baseline = baseline
        self.igram = []
        self.los_displacements1 = []
        self.los_displacements2 = []
        self.longitudes = []
        self.latitudes = []

    def save(self):
        """
        Save the interferogram to an HDF5 file
        :return: Nothing returned
        """

        # Open HDF5 file
        repo = git.Repo('.', search_parent_directories=True)
        f = h5py.File(name=repo.working_tree_dir + '/files/interferogram.hdf5', mode="w")

        # Get data shape
        igram_shape = [len(row) for row in self.igram]

        # Save data shape
        f.create_dataset(name="data_shape", data=igram_shape, chunks=True)

        # Flatten igram
        flattened_igram = [item for row in self.igram for item in row]

        # Save flat igram data
        f.create_dataset(name="flattened_igram", data=flattened_igram, chunks=True)

        # Flatten LOS displacements 1
        flattened_los_displacements1 = [item for row in self.los_displacements1 for item in row]

        # Save displacements 1
        f.create_dataset(name="flattened_los_displacements1", data=flattened_los_displacements1, chunks=True)

        # Flatten LOS displacements 2
        flattened_los_displacements2 = [item for row in self.los_displacements2 for item in row]

        # Save displacements 2
        f.create_dataset(name="flattened_los_displacements2", data=flattened_los_displacements2, chunks=True)

        # Flatten longitudes
        flattened_longitudes = [item for row in self.longitudes for item in row]

        # Save flat longitude data
        f.create_dataset(name="flattened_longitudes", data=flattened_longitudes, chunks=True)

        # Flatten latitudes
        flattened_latitudes = [item for row in self.latitudes for item in row]

        # Save flat latitude data
        f.create_dataset(name="flattened_latitudes", data=flattened_latitudes, chunks=True)

        # Save time space
        f.create_dataset(name="time_space1", data=self.time_space1, chunks=True)

        # Save flat second time space
        f.create_dataset(name="time_space2", data=self.time_space2, chunks=True)

        # Save baseline
        f.attrs["baseline"] = self.baseline

    def load(self):
        """
        Load a hdf5 file containing a previous SimpleInterferogram object
        :return: Nothing returned
        """

        # Open the HDF5 file in read mode
        repo = git.Repo('.', search_parent_directories=True)
        f = h5py.File(name=repo.working_tree_dir + '/files/interferogram.hdf5', mode='r')

        # Get the data shape
        data_shape = f["data_shape"]

        # Load all the saved data
        flattened_igram = f["flattened_igram"]
        flattened_los_displacements1 = f["flattened_los_displacements1"]
        flattened_los_displacements2 = f["flattened_los_displacements1"]
        flattened_longitudes = f["flattened_longitudes"]
        flattened_latitudes = f["flattened_latitudes"]
        self.time_space1 = f["flattened_time_space1"]
        self.time_space2 = f["flattened_time_space2"]
        self.baseline = f.attrs["baseline"]

        # Use the date shape and known structure to populate the object instance variables
        index = 0
        for length in data_shape:
            self.igram.append(flattened_igram[index:index + length])
            self.los_displacements1.append(flattened_los_displacements1[index:index + length])
            self.los_displacements2.append(flattened_los_displacements2[index:index + length])
            self.longitudes.append(flattened_longitudes[index:index + length])
            self.latitudes.append(flattened_latitudes[index:index + length])
            index += length

    def get_angles_cartesian(self, time, satellite_position, flat_positions):
        """
        Get the look angles from the satellite to the flattened positions
        :param time: Times of satellite
        :param satellite_position: Positions of satellite at times
        :param flat_positions: Flattened ellipsoid points
        :return: Look angles in radians
        """

        # Calculate the satellite intersect and height with the shape
        intersect, satellite_height = self.planet.get_sub_obs_points(times=time, instrument=self.instrument)

        # Calculate the distance between center and satellite intersects
        satellite_radius = np.asarray(np.linalg.norm(intersect - [0, 0, 0]))

        # Get the distance between the satellite and the surface points
        distances = np.asarray([np.linalg.norm(flat_position - satellite_position)
                                for flat_position in flat_positions])

        # Calculate cosines of the angles
        z_plus_re = satellite_height + satellite_radius
        altitudes = np.asarray([np.linalg.norm(flat_position - [0, 0, 0])
                                for flat_position in flat_positions])
        cosine_look_angle = (distances ** 2 + z_plus_re ** 2 - altitudes ** 2) / (2 * distances * z_plus_re)

        # Use arc cosine to get the angles in radians
        look_angles_radians = np.arccos(cosine_look_angle)

        # Return the look angles 2-D array
        return look_angles_radians

    def get_flattened_angles(self, positions, time, satellite_position):
        """
        Calculate the flattened angle from a given set of positions at a time in the satellite orbit
        :param positions: Ground positions
        :param time: Time at which the satellite is in satellite_position
        :param satellite_position: Position of the satellite at a specific time
        :return: Flattened angles of satellite to ground positions
        """

        # Get planet axes
        a, b, c = self.planet.get_axes()

        # Calculate vectors from the ground positions to the satellite
        positions = np.asarray(positions)
        ground_satellite_vectors = positions - satellite_position

        # Create plates normal to the vector from the ground to the satellite at the ground position
        planes = spice.nvp2pl_vector(normal=ground_satellite_vectors, point=positions)

        # Generate ellipses from these planes that intersect the planet ellipsoid
        ellipses = spice.inedpl_vector(a=a, b=b, c=c, plane=planes)[0]

        # Get the "flattened" points from the intersection of ellipses and ground positions
        flat_points = spice.npelpt_vector(point=positions, ellips=ellipses)[0]

        # Calculate the flat angles at the time of the satellite position given the flattened points
        flat_angles = self.get_angles_cartesian(time=time,
                                                satellite_position=satellite_position, flat_positions=flat_points)

        return flat_angles

    def calculate_igram(self):
        """
        Calculate the flattened phases between two swaths given a baseline
        :return: Nothing returned 
        """

        total_longitudes = []
        total_latitudes = []
        total_phases = []
        total_los_displacements1 = []
        total_los_displacements2 = []

        # Get planet axes
        a, b, c = self.planet.get_axes()

        # For speed, calculate all points at once
        u_displacements, v_displacements, w_displacements, lats, longs = (
            self.deformation_map.get_displacements(time_space=self.time_space1, swath=self.swath1))
        transformed_displacements = np.array([u_displacements, v_displacements, w_displacements]).T
        displacements_1 = []
        latitudes = []
        longitudes = []
        index = 0

        for sublist in self.swath1:
            sublist_length = len(sublist)  # Number of elements in this sublist
            displacements_1.append(transformed_displacements[index:index + sublist_length])
            latitudes.append(lats[index:index + sublist_length])
            longitudes.append(longs[index:index + sublist_length])
            index += sublist_length

        u_displacements, v_displacements, w_displacements, lats, longs = (
            self.deformation_map.get_displacements(time_space=self.time_space2, swath=self.swath2))
        transformed_displacements = np.array([u_displacements, v_displacements, w_displacements]).T
        displacements_2 = []
        index = 0

        for sublist in self.swath2:
            sublist_length = len(sublist)  # Number of elements in this sublist
            displacements_2.append(transformed_displacements[index:index + sublist_length])
            index += sublist_length

        for i in range(len(self.swath1)):

            # Get the positions of the groundTargets
            positions = np.asarray([(point[0], point[1], point[2]) for point in self.swath1[i]])

            # Get satellite positions and velocities
            satellite_position, sat_velocity = self.instrument.get_states(times=self.time_space1[i])

            # Get flattened angles for the satellite positions and ground points
            angles = self.get_flattened_angles(positions=positions, time=self.time_space1[i],
                                               satellite_position=satellite_position)

            # Calculate geodetic heights of each of the ground points
            geodetic_heights = spice.nearpt_vector(positn=positions, a=a, b=b, c=c)[1]

            # Get the distances from the ground points to the satellites
            distances = np.asarray([np.linalg.norm(position - satellite_position) for position in positions])

            # Calculate the vectors from the ground points to the satellite positions
            ground_satellite_vectors = satellite_position - positions

            # Calculate the line of sight displacements given the satellite positions for the first swath
            los_displacements1 = np.asarray([(np.dot(displacement, ground_satellite_vector) /
                                              np.linalg.norm(ground_satellite_vector))
                                             for displacement, ground_satellite_vector
                                             in zip(displacements_1[i], ground_satellite_vectors)])

            # Calculate the line of sight displacements given the satellite positions for the second swath
            los_displacements2 = np.asanyarray([(np.dot(displacement, ground_satellite_vector) /
                                                 np.linalg.norm(ground_satellite_vector))
                                                for displacement, ground_satellite_vector
                                                in zip(displacements_2[i], ground_satellite_vectors)])

            # Calculate the phases using the LOS displacements, baseline, geodetic heights, distances, angles
            phases = 4 * np.pi * (1 / self.instrument.wavelength) * (
                    (los_displacements2 - los_displacements1) -
                    ((self.baseline * geodetic_heights) /
                     (distances * np.sin(angles))))

            # Modulate phases to be in the 0 to 2 pi range
            for k in range(len(phases)):
                phases[k] = np.fmod(phases[k], 2 * np.pi)
                if phases[k] < 0:
                    phases[k] = phases[k] + (2 * np.pi)

            # Append the phases, latitudes, and longitudes
            total_phases.append(phases)
            total_latitudes.append([latitude for latitude in latitudes[i]])
            total_longitudes.append([longitude for longitude in longitudes[i]])
            total_los_displacements1.append(los_displacements1)
            total_los_displacements2.append(los_displacements2)

        self.igram = total_phases
        self.los_displacements1 = total_los_displacements1
        self.los_displacements2 = total_los_displacements2
        self.longitudes = total_longitudes
        self.latitudes = total_latitudes
        return self.igram, self.los_displacements1, self.los_displacements2, self.longitudes, self.latitudes

    def visualize_interferogram(self, projection, fig=None, globe=None, ax=None, return_fig=False):
        """
        Visualize the interferogram
        :param projection: Projection object to use for plotting
        :return: Nothing returned
        """

        if fig is None:
            # Get the projection
            fig, ax, globe = projection.proj(planet=self.planet)

        # Make the colormap cyclical
        cm = plt.cm.get_cmap('hsv')

        # Iterate through the interferogram beams
        for i in range(len(self.igram)):
            # Plot points on the map
            im = ax.scatter(self.longitudes[i], self.latitudes[i], vmin=0, vmax=2 * np.pi, cmap=cm,
                            transform=ccrs.PlateCarree(globe=globe), c=self.igram[i], marker='o', s=0.1)

        # return fig, ax, globe if necessary
        if return_fig:
            return fig, ax, globe

        # Add colorbar
        plt.colorbar(im, label="Phase")

        # Add labels and legend
        ax.set_title('Interferogram', pad=20)

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/' + 'interferogram.png', format='png', dpi=500)

    def visualize_displacements(self, projection, fig=None, globe=None, ax=None, return_fig=False):
        """
        Visualize the interferogram
        :param projection: Projection object to use for plotting
        :return: Nothing returned
        """

        if fig is None:
            # Get the projection
            fig, ax, globe = projection.proj(planet=self.planet)

        # Iterate through the interferogram beams
        for i in range(len(self.igram)):
            # Plot points on the map
            im = ax.scatter(self.longitudes[i], self.latitudes[i],
                            transform=ccrs.PlateCarree(globe=globe),
                            c=(self.los_displacements2[i] - self.los_displacements1[i]), marker='o', s=0.1)

        # return fig, ax, globe if necessary
        if return_fig:
            return fig, ax, globe

        # Add colorbar
        plt.colorbar(im, label="Phase")

        # Add labels and legend
        ax.set_title('Interferogram', pad=20)

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/' + 'displacements.png', format='png', dpi=500)

# end of file
