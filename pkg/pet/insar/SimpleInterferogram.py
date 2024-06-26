#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import copy
import cspyce as spice
import snaphu
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
    pairing_two.default = 0
    pairing_two.doc = "Second instance of track number"

    track_number = pet.properties.int()
    track_number.default = 0
    track_number.doc = "Track number to get first and second index of"

    def __init__(self, name, locator, implicit, planet, instrument, displacements, load_path=None, time_interval=10,
                 ground_resolution=2000, baseline=10):
        super().__init__(name, locator, implicit)
        self.planet = planet
        self.instrument = instrument
        self.displacement = displacements
        self.igram = []
        self.longitudes = []
        self.latitudes = []
        self.groundSwath = None
        if load_path is not None:
            self.load(load_path)
        else:
            self.calculate_igram(displacements=displacements,
                                 time_interval=time_interval, ground_resolution=ground_resolution, baseline=baseline)

    def save(self, path):
        """
        Save the interferogram to an HDF5 file
        :param path: Path where to save interferogram HDF5 file
        """

        # Open HDF5 file
        f = h5py.File(name=path + '/interferogram.hdf5', mode="w")

        # Get data shape
        igram_shape = [len(row) for row in self.igram]

        # Save data shape
        f.create_dataset(name="data_shape", data=igram_shape, chunks=True)

        # Flatten igram
        flattened_igram = [item for row in self.igram for item in row]

        # Save flat igram data
        f.create_dataset(name="flattened_igram", data=flattened_igram, chunks=True)

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

        """

        f = h5py.File(name=path + '/interferogram.hdf5', mode='r')

        data_shape = f["data_shape"]

        flattened_igram = f["flattened_igram"]
        index = 0
        for length in data_shape:
            self.igram.append(flattened_igram[index:index + length])
            index += length

        flattened_longitudes = f["flattened_longitudes"]
        index = 0
        for length in data_shape:
            self.longitudes.append(flattened_longitudes[index:index + length])
            index += length

        flattened_latitudes = f["flattened_latitudes"]
        index = 0
        for length in data_shape:
            self.latitudes.append(flattened_latitudes[index:index + length])
            index += length

        self.pairing_one = f.attrs["pairing_one"]
        self.pairing_two = f.attrs["pairing_two"]
        self.track_number = f.attrs["track_number"]

        start_time = f["groundswath"].attrs["start_time"]
        end_time = f["groundswath"].attrs["end_time"]
        time_interval = f["groundswath"].attrs["time_interval"]
        ground_resolution = f["groundswath"].attrs["ground_resolution"]

        time_space = f["groundswath/time_space"]

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

        """

        # Get planet axes
        a, b, c = self.planet.get_axes()

        positions = np.asarray(positions)
        ground_satellite_vectors = positions - satellite_position

        planes = spice.nvp2pl_vector(normal=ground_satellite_vectors, point=positions)

        ellipses = spice.inedpl_vector(a=a, b=b, c=c, plane=planes)[0]

        flat_points = spice.npelpt_vector(point=positions, ellips=ellipses)[0]

        flat_angles = self.get_angles_cartesian(time=time,
                                                satellite_position=satellite_position, flat_positions=flat_points)

        return flat_angles

    def calculate_flattened_phases(self, ground_swath1, ground_swath2, baseline):
        """

        """

        longitudes = []
        latitudes = []
        total_phases = []

        # Get planet axes
        a, b, c = self.planet.get_axes()

        # Create a coordinate conversion object
        convert = pet.ext.conversions(name="conversions", a=a, b=b, c=c)

        for i in range(len(ground_swath1.swath_beams)):

            positions = np.asarray([point.get_position() for point in ground_swath1.swath_beams[i]])

            geodetic_coordinates = convert.geodetic(cartesian_coordinates=positions)[:, :3]

            satellite_position, sat_velocity = self.instrument.get_states(times=ground_swath1.time_space[i])

            angles = self.get_flattened_angles(positions=positions, time=ground_swath1.time_space[i],
                                               satellite_position=satellite_position)

            geodetic_heights = spice.nearpt_vector(positn=positions, a=a, b=b, c=c)[1]

            distances = np.asarray([np.linalg.norm(position - [0, 0, 0]) for position in positions])

            ground_satellite_vectors = satellite_position - positions

            displacements_1 = [point.displacement for point in ground_swath1.swath_beams[i]]

            los_displacements1 = (
                np.asarray([np.dot(displacement, ground_satellite_vector) / np.linalg.norm(ground_satellite_vector)
                            for displacement, ground_satellite_vector
                            in zip(displacements_1, ground_satellite_vectors)]))

            displacements_2 = [point.displacement for point in ground_swath2.swath_beams[i]]
            los_displacements2 = (
                np.asarray([np.dot(displacement, ground_satellite_vector) / np.linalg.norm(ground_satellite_vector)
                            for displacement, ground_satellite_vector
                            in zip(displacements_2, ground_satellite_vectors)]))

            phases = 4 * np.pi * (1 / self.instrument.wavelength) * (
                    (los_displacements2 - los_displacements1) -
                    ((baseline * geodetic_heights) /
                     (distances * np.sin(angles))))

            for k in range(len(phases)):
                phases[k] = np.fmod(phases[k], 2 * np.pi)
                if phases[k] < 0:
                    phases[k] = phases[k] + (2 * np.pi)

            total_phases = np.concatenate((total_phases, phases))
            latitudes = np.concatenate((latitudes, geodetic_coordinates[:, 0]))
            longitudes = np.concatenate((longitudes, geodetic_coordinates[:, 1]))

        self.igram = total_phases
        self.longitudes = longitudes
        self.latitudes = latitudes

    def calculate_igram(self, displacements, time_interval, ground_resolution, baseline):
        """

        """

        times = self.instrument.get_five_tracks()

        orbit_cycle_time = self.instrument.orbit_cycle

        gs1 = GroundSwath(start_time=times[self.track_number - 1],
                                    end_time=times[self.track_number], time_interval=time_interval,
                                    ground_resolution=ground_resolution, planet=self.planet,
                                    instrument=self.instrument)

        self.groundSwath = gs1

        gs2 = copy.copy(gs1)

        displacements.attach(swath=gs1, use_mid_point=True)

        displacements.attach(swath=gs2, use_mid_point=True, time_displacement=orbit_cycle_time)

        self.calculate_flattened_phases(ground_swath1=gs1, ground_swath2=gs2, baseline=baseline)

    def recalculate_igram(self, baseline):
        """

        """

        orbit_cycle_time = self.instrument.orbit_cycle

        gs1 = self.groundSwath

        gs2 = copy.deepcopy(gs1)

        self.displacements.attach(swath=gs1, use_mid_point=True)

        self.displacements.attach(swath=gs2, use_mid_point=True, time_displacement=orbit_cycle_time)

        self.calculate_flattened_phases(ground_swath1=gs1, ground_swath2=gs2, baseline=baseline)

    def unwrap(self):
        """

        """

        # Sample coherence for an interferogram with no noise.
        corr = np.ones(self.igram.shape, dtype=np.float32)

        # Unwrap using the 'SMOOTH' cost mode and 'MCF' initialization method.
        unw, conncomp = snaphu.unwrap(self.igram.shape, corr, nlooks=1.0, cost="smooth", init="mcf")

        return unw, conncomp

    def visualize(self, projection):
        """

        """

        # Get the fig, ax, globe from the instrument orbit
        fig, ax, globe = self.instrument.plot_orbit(projection=projection,
                                                    start_time=self.groundSwath.start_time,
                                                    end_time=self.groundSwath.end_time,
                                                    time_interval=self.groundSwath.time_interval, return_fig=True)

        cm = plt.cm.get_cmap('hsv')

        for i in range(len(self.igram)):

            # Plot points on the map
            im = ax.scatter(self.longitudes[i], self.latitudes[i], vmin=0, vmax=2 * np.pi, cmap=cm,
                            transform=ccrs.PlateCarree(globe=globe), c=self.igram[i], marker='o', s=0.1)

        # Add colorbar
        plt.colorbar(im, label="Phase")

        # Add labels and legend
        ax.set_title('Interferogram')

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/interferogram.png', format='png', dpi=500)

# end of file
