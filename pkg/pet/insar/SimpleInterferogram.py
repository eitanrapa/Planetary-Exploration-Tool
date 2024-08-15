#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import copy
import cspyce as spice
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import h5py
from .GroundTarget import GroundTarget
from .GroundSwath import GroundSwath


class SimpleInterferogram(pet.component):
    """
    Class that creates a single interferogram between two points of a displacement map given an instrument orbit
    """

    pairing_one = pet.properties.int()
    pairing_one.default = 0
    pairing_one.doc = "First instance of track number"

    pairing_two = pet.properties.int()
    pairing_two.default = 1
    pairing_two.doc = "Second instance of track number"

    track_number = pet.properties.int()
    track_number.default = 0
    track_number.doc = "Track number to get first and second index of"

    def __init__(self, name, locator, implicit, planet, instrument, displacements, load_path=None, time_interval=10,
                 ground_resolution=2000, baseline=10):
        super().__init__(name, locator, implicit)
        self.planet = planet
        self.instrument = instrument
        self.displacements = displacements
        self.igram = []
        self.los_displacements1 = []
        self.los_displacements2 = []
        self.satellite_vectors = []
        self.displacements1 = []
        self.displacements2 = []
        self.longitudes = []
        self.latitudes = []
        self.groundSwath = None
        self.base_time_space = None
        if load_path is not None:
            self.load(load_path)
        else:
            self.calculate_igram(time_interval=time_interval, ground_resolution=ground_resolution, baseline=baseline)

    def save(self, path):
        """
        Save the interferogram to an HDF5 file
        :param path: Path where to save interferogram HDF5 file
        :return: Nothing returned
        """

        # Open HDF5 file
        f = h5py.File(name=path + '/' + self.displacements.pyre_name + '_' + str(self.track_number) + '_' +
                      str(self.pairing_one) + '_' + str(self.pairing_two) + '_interferogram.hdf5', mode="w")

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

        # Flatten satellite vectors
        flattened_satellite_vectors = [item for row2d in self.satellite_vectors for row in row2d for item in row]

        # Save satellite vectors
        f.create_dataset(name="flattened_satellite_vectors", data=flattened_satellite_vectors, chunks=True)

        # Flatten ground displacements 1
        flattened_displacements1 = [item for row2d in self.displacements1 for row in row2d for item in row]

        # Save ground displacements 1
        f.create_dataset(name="flattened_displacements1", data=flattened_displacements1, chunks=True)

        # Flatten ground displacements 2
        flattened_displacements2 = [item for row2d in self.displacements2 for row in row2d for item in row]

        # Save ground displacements 2
        f.create_dataset(name="flattened_displacements2", data=flattened_displacements2, chunks=True)

        # Flatten longitudes
        flattened_longitudes = [item for row in self.longitudes for item in row]

        # Save flat longitude data
        f.create_dataset(name="flattened_longitudes", data=flattened_longitudes, chunks=True)

        # Flatten latitudes
        flattened_latitudes = [item for row in self.latitudes for item in row]

        # Save flat latitude data
        f.create_dataset(name="flattened_latitudes", data=flattened_latitudes, chunks=True)

        # Save pairing attributes
        f.attrs["pairing_one"] = self.pairing_one
        f.attrs["pairing_two"] = self.pairing_two

        # Save track number attribute
        f.attrs["track_number"] = self.track_number

        # Save the base time space
        f.create_dataset(name="flattened_base_time_space", data=self.base_time_space, chunks=True)

        # Create a groundSwath group
        grp1 = f.create_group(name="groundswath")

        # Save groundSwath attributes
        grp1.attrs["start_time"] = self.groundSwath.start_time
        grp1.attrs["end_time"] = self.groundSwath.end_time
        grp1.attrs["time_interval"] = self.groundSwath.time_interval
        grp1.attrs["ground_resolution"] = self.groundSwath.ground_resolution

        # Save groundSwath time space
        grp1.create_dataset(name="time_space", data=self.groundSwath.time_space, chunks=True)

        # Flatten the beams
        flattened_beams = [item.get_position() for row in self.groundSwath.swath_beams for item in row]

        # Save the beams
        grp1.create_dataset(name="flattened_beams", data=flattened_beams, chunks=True)

    def load(self, path):
        """
        Load a hdf5 file containing a previous SimpleInterferogram object
        :param path: Path to file to load
        :return: Nothing returned
        """

        # Open the HDF5 file in read mode
        f = h5py.File(name=path, mode='r')

        # Get the data shape
        data_shape = f["data_shape"]

        # Load all the saved data
        flattened_igram = f["flattened_igram"]
        flattened_los_displacements1 = f["flattened_los_displacements1"]
        flattened_los_displacements2 = f["flattened_los_displacements1"]
        flattened_satellite_vectors = f["flattened_satellite_vectors"]
        flattened_displacements1 = f["flattened_displacements1"]
        flattened_displacements2 = f["flattened_displacements2"]
        flattened_longitudes = f["flattened_longitudes"]
        flattened_latitudes = f["flattened_latitudes"]

        # Use the date shape and known structure to populate the object instance variables
        index = 0
        index3d = 0
        for length in data_shape:
            self.igram.append(flattened_igram[index:index + length])
            self.los_displacements1.append(flattened_los_displacements1[index:index + length])
            self.los_displacements2.append(flattened_los_displacements2[index:index + length])
            self.satellite_vectors.append(np.asanyarray((flattened_satellite_vectors[index3d:index3d + length:3],
                                                         flattened_satellite_vectors[index3d +
                                                                                     1:(index3d + length) * 3:3],
                                                         flattened_satellite_vectors[index3d +
                                                                                     2:(index3d + length) * 3:3])).T)
            self.displacements1.append(np.asanyarray((flattened_displacements1[index3d:index3d + length:3],
                                                      flattened_displacements1[index3d + 1:index3d + length:3],
                                                      flattened_displacements1[index3d + 2:index3d + length:3])).T)
            self.displacements2.append(np.asanyarray((flattened_displacements2[index3d:index3d + length:3],
                                                      flattened_displacements2[index3d + 1:index3d + length:3],
                                                      flattened_displacements2[index3d + 2:index3d + length:3])).T)
            self.longitudes.append(flattened_longitudes[index:index + length])
            self.latitudes.append(flattened_latitudes[index:index + length])
            index += length
            index3d = index3d + length * 3

        # Get the pairing numbers, track number, and base time space
        self.pairing_one = f.attrs["pairing_one"]
        self.pairing_two = f.attrs["pairing_two"]
        self.track_number = f.attrs["track_number"]
        self.base_time_space = f["flattened_base_time_space"]

        # Load the groundSwath parameters
        start_time = f["groundswath"].attrs["start_time"]
        end_time = f["groundswath"].attrs["end_time"]
        time_interval = f["groundswath"].attrs["time_interval"]
        ground_resolution = f["groundswath"].attrs["ground_resolution"]

        time_space = f["groundswath/time_space"]

        # Populate the GroundSwath with GroundTarget objects
        flattened_beams = f["groundswath/flattened_beams"]
        swath_beams = []
        index = 0
        for length in data_shape:
            swath_beams.append([GroundTarget(x=element[0], y=element[1], z=element[2]) for
                                element in flattened_beams[index:index + length]])
            index += length

        self.groundSwath = GroundSwath(start_time=start_time, end_time=end_time, time_interval=time_interval,
                                       ground_resolution=ground_resolution, planet=self.planet,
                                       instrument=self.instrument, copy=True)
        self.groundSwath.time_space = time_space
        self.groundSwath.swath_beams = swath_beams

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

    def calculate_flattened_phases(self, ground_swath1, ground_swath2, baseline):
        """
        Calculate the flattened phases between two swaths given a baseline and an attached displacement map
        :param ground_swath1: First GroundSwath object
        :param ground_swath2: Second GroundSwath object
        :param baseline: Baseline of interferogram
        :return: Nothing returned 
        """

        longitudes = []
        latitudes = []
        total_phases = []

        # Get planet axes
        a, b, c = self.planet.get_axes()

        for i in range(len(ground_swath1.swath_beams)):

            # Get the positions of the groundTargets
            positions = np.asarray([point.get_position() for point in ground_swath1.swath_beams[i]])

            # Get satellite positions and velocities
            satellite_position, sat_velocity = self.instrument.get_states(times=ground_swath1.time_space[i])

            # Get flattened angles for the satellite positions and ground points
            angles = self.get_flattened_angles(positions=positions, time=ground_swath1.time_space[i],
                                               satellite_position=satellite_position)

            # Calculate geodetic heights of each of the ground points
            geodetic_heights = spice.nearpt_vector(positn=positions, a=a, b=b, c=c)[1]

            # Get the distances from the ground points to the satellites
            distances = np.asarray([np.linalg.norm(position - satellite_position) for position in positions])

            # Calculate the vectors from the ground points to the satellite positions
            ground_satellite_vectors = satellite_position - positions

            # Get the displacements from the first swath
            displacements_1 = [point.displacement for point in ground_swath1.swath_beams[i]]

            # Get the displacements from the second swath
            displacements_2 = [point.displacement for point in ground_swath2.swath_beams[i]]

            # Calculate the line of sight displacements given the satellite positions for the first swath
            los_displacements1 = np.asarray([(np.dot(displacement[:3], ground_satellite_vector) /
                                              np.linalg.norm(ground_satellite_vector))
                                             for displacement, ground_satellite_vector
                                             in zip(displacements_1, ground_satellite_vectors)])

            # Calculate the line of sight displacements given the satellite positions for the second swath
            los_displacements2 = np.asanyarray([(np.dot(displacement[:3], ground_satellite_vector) /
                                                 np.linalg.norm(ground_satellite_vector))
                                                for displacement, ground_satellite_vector
                                                in zip(displacements_2, ground_satellite_vectors)])

            # Calculate the phases using the LOS displacements, baseline, geodetic heights, distances, angles
            phases = 4 * np.pi * (1 / self.instrument.wavelength) * (
                    (los_displacements2 - los_displacements1) -
                    ((baseline * geodetic_heights) /
                     (distances * np.sin(angles))))

            # Modulate phases to be in the 0 to 2 pi range
            for k in range(len(phases)):
                phases[k] = np.fmod(phases[k], 2 * np.pi)
                if phases[k] < 0:
                    phases[k] = phases[k] + (2 * np.pi)

            # Append the phases, latitudes, and longitudes
            total_phases.append(phases)
            latitudes.append([displacement[3] for displacement in displacements_1])
            longitudes.append([displacement[4] for displacement in displacements_1])

        self.igram = total_phases
        self.longitudes = longitudes
        self.latitudes = latitudes

    def calculate_igram(self, time_interval, ground_resolution, baseline):
        """
        Calculate the interferogram
        :param time_interval: Interval of time sampling of satellite
        :param ground_resolution: Azimuthal ground resolution per beam
        :param baseline: Baselines for interferometry
        :return: Nothing returned
        """

        # Get the times defining the first five tracks
        times = self.instrument.get_five_tracks()

        # Get the orbit cycle time of the instrument
        orbit_cycle_time = self.instrument.orbit_cycle

        # Create a groundSwath object from the track times, time interval, ground resolution, planet, and instrument
        self.groundSwath = GroundSwath(start_time=times[self.track_number],
                                       end_time=times[self.track_number + 1], time_interval=time_interval,
                                       ground_resolution=ground_resolution, planet=self.planet,
                                       instrument=self.instrument)

        # Define the base time space
        self.base_time_space = self.groundSwath.time_space

        # Set the correct time space for the first GroundSwath given the first pairing number
        self.groundSwath.time_space = self.base_time_space + orbit_cycle_time * self.pairing_one

        # Create the second groundSwath object
        gs2 = GroundSwath(start_time=self.groundSwath.start_time,
                          end_time=self.groundSwath.end_time, time_interval=self.groundSwath.time_interval,
                          ground_resolution=self.groundSwath.ground_resolution, planet=self.planet,
                          instrument=self.instrument,
                          copy=True)

        gs2.swath_beams = copy.deepcopy(self.groundSwath.swath_beams)

        # Set the correct time space for the second GroundSwath given the second pairing number
        gs2.time_space = self.base_time_space + orbit_cycle_time * self.pairing_two

        # Attach the displacements of both GroundSwaths
        self.displacements.attach(swath1=self.groundSwath, swath2=gs2, use_mid_point=True)

        # Calculate the flattened phases from the GroundSwaths and baseline
        self.calculate_flattened_phases(ground_swath1=self.groundSwath, ground_swath2=gs2, baseline=baseline)

    def recalculate_igram(self, baseline, pairing_one=None, pairing_two=None):
        """
        Calculate the new interferogram using a saved interferogram
        :param baseline: Baselines for interferometry
        :param pairing_one: Change first pairing if necessary
        :param pairing_two: Change second pairing if necessary
        :return: Nothing returned
        """

        # Get the instrument orbit cycle
        orbit_cycle_time = self.instrument.orbit_cycle

        # Change the pairings if specified
        if pairing_one is not None:
            self.pairing_one = pairing_one
        if pairing_two is not None:
            self.pairing_two = pairing_two

        # Calculate the first GroundSwath time space according to the first pairing number
        self.groundSwath.time_space = self.base_time_space + orbit_cycle_time * self.pairing_one

        # Create the second groundSwath object
        gs2 = GroundSwath(start_time=self.groundSwath.start_time,
                          end_time=self.groundSwath.end_time, time_interval=self.groundSwath.time_interval,
                          ground_resolution=self.groundSwath.ground_resolution, planet=self.planet,
                          instrument=self.instrument,
                          copy=True)

        gs2.swath_beams = copy.deepcopy(self.groundSwath.swath_beams)

        # Calculate the second GroundSwath time space according to the second pairing number
        gs2.time_space = self.base_time_space + orbit_cycle_time * self.pairing_two

        # Attach the displacements of both GroundSwaths
        self.displacements.attach(swath1=self.groundSwath, swath2=gs2, use_mid_point=True)

        # Calculate the flattened phases from the GroundSwaths and baseline
        self.calculate_flattened_phases(ground_swath1=self.groundSwath, ground_swath2=gs2, baseline=baseline)

    def visualize(self, projection):
        """
        Visualize the interferogram
        :param projection: Projection object to use for plotting
        :return: Nothing returned
        """

        # Get the fig, ax, globe from the instrument orbit
        fig, ax, globe = self.instrument.plot_orbit(projection=projection,
                                                    start_time=self.groundSwath.start_time,
                                                    end_time=self.groundSwath.end_time,
                                                    time_interval=self.groundSwath.time_interval, return_fig=True)

        # Make the colormap cyclical
        cm = plt.cm.get_cmap('hsv')

        # Iterate through the interferogram beams
        for i in range(len(self.igram)):

            # Plot points on the map
            im = ax.scatter(self.longitudes[i], self.latitudes[i], vmin=0, vmax=2 * np.pi, cmap=cm,
                            transform=ccrs.PlateCarree(globe=globe), c=self.igram[i], marker='o', s=0.1)

        # Add colorbar
        plt.colorbar(im, label="Phase")

        # Add labels and legend
        ax.set_title('Interferogram')

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/' + self.displacements.pyre_name + '_' + str(self.track_number) +
                    '_' + str(self.pairing_one) + '_' + str(self.pairing_two) + '_interferogram.png', format='png',
                    dpi=500)

# end of file
