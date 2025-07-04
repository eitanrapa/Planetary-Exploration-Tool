#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import sys
import cartopy.crs as ccrs
import cspyce as spice
import matplotlib.pyplot as plt
import numpy as np
import pet
import xarray as xr
from tqdm import tqdm


class Track(pet.component, family="pet.dataAcquisition.track", implements=pet.protocols.dataAcquisition):
    """
    An object containing the swath of a satellite pass. Contains the start, end, and temporal_resolution of the
    satellite pass
    as well as given ground resolution to calculate vectors with. Also contains the satellite time vector and an
    xarray dataset of ground intersects and times
    """

    start_time = pet.properties.float()
    start_time.doc = "start time of track"

    end_time = pet.properties.float()
    end_time.doc = "end time of track"

    planet = pet.protocols.planet()
    planet.doc = "target planet"

    campaign = pet.protocols.campaigns.orbiter()
    campaign.doc = "campaign (must be an orbiter)"

    instrument = pet.protocols.instruments.inSAR()
    instrument.doc = "observation instrument"

    temporal_resolution = pet.properties.float()
    temporal_resolution.default = 10
    temporal_resolution.doc = "temporal resolution of orbit [s]"

    spatial_resolution = pet.properties.float()
    spatial_resolution.default = 200
    spatial_resolution.doc = "spatial resolution of ground track [m]"

    data = None

    @classmethod
    def from_file(cls, planet, campaign, instrument, file_name):
        """
        Load a track from an HDF5 file
        :param planet: Planet object
        :param campaign: Campaign object
        :param instrument: Instrument object
        :param file_name: File name of HDF5 file
        :return: Track object
        """

        # Open the HDF5 file in read mode
        data = xr.open_dataarray(filename_or_obj=file_name)

        # Get the start and end times
        start_time = data.attrs["start_time"]
        end_time = data.attrs["end_time"]

        # Create the object
        obj = cls(name="track" + str(np.random.rand()), start_time=start_time, end_time=end_time,
                  planet=planet, campaign=campaign, instrument=instrument,
                  temporal_resolution=data.attrs["temporal_resolution"],
                  spatial_resolution=data.attrs["spatial_resolution"])

        obj.data = data  # Restore computed result

        return obj

    @classmethod
    def from_data_array(cls, planet, campaign, instrument, data):
        """
        Create a track object from a data array
        :param planet: Planet object
        :param campaign: Campaign object
        :param instrument: Instrument object
        :param data: Data array
        :return: Track object
        """

        # Get the start and end times
        start_time = data.attrs["start_time"]
        end_time = data.attrs["end_time"]

        # Create the object
        obj = cls(name="track" + str(start_time) + str(end_time),
                  start_time=start_time, end_time=end_time,
                  planet=planet, campaign=campaign, instrument=instrument,
                  temporal_resolution=data.attrs["temporal_resolution"],
                  spatial_resolution=data.attrs["spatial_resolution"])

        obj.data = data  # Restore computed result

        return obj

    @classmethod
    def from_files(cls, planet, campaign, instrument, file_list):
        """
        Load a list of tracks from HDF5 files
        :param planet: Planet object
        :param campaign: Campaign object
        :param instrument: Instrument object
        :param file_list: List of file names
        :return: List of track objects
        """

        # Load all the files
        return [cls.from_file(planet=planet, campaign=campaign, instrument=instrument, file_name=file) for
                file in file_list]

    def save(self, file_name):
        """
        Save the track to an HDF5 file
        :return: Nothing returned
        """

        # Open HDF5 file
        self.data.to_netcdf(file_name, engine="netcdf4")

    def create_data_array(self, cartesian_coordinates, geodetic_coordinates, look_angles, incidence_angles,
                          non_projected_incidence_angles, times):
        """
        Create a xarray with the input data
        :param cartesian_coordinates: x, y, z coordinates to be saved
        :param geodetic_coordinates: lat, long, height coordinates to be saved
        :param look_angles: Calculated look angles from satellite
        :param incidence_angles: Calculated incidence angles from satellite
        :param non_projected_incidence_angles: Calculated non-projected incidence angles from satellite
        :param times: Times of observation for points
        :return: Nothing returned
        """

        # Create the xarray Dataset
        da = xr.DataArray(
            data=look_angles,
            dims=["points"],
            coords={
                "time": ("points", times),
                "x": ("points", np.asarray([point[0] for point in cartesian_coordinates])),
                "y": ("points", np.asarray([point[1] for point in cartesian_coordinates])),
                "z": ("points", np.asarray([point[2] for point in cartesian_coordinates])),
                "latitude": ("points", np.asarray([point[0] for point in geodetic_coordinates])),
                "longitude": ("points", np.asarray([point[1] for point in geodetic_coordinates])),
                "height": ("points", np.asarray([point[2] for point in geodetic_coordinates])),
                "incidence_angle": ("points", np.asarray(incidence_angles)),
                "non_projected_incidence_angle": ("points", np.asarray(non_projected_incidence_angles))},
            name="look_angles",
            attrs=dict(
                body_id=self.campaign.body_id,
                start_time=self.start_time,
                end_time=self.end_time,
                temporal_resolution=self.temporal_resolution,
                spatial_resolution=self.spatial_resolution
            ),
        )

        # Save xarray to object
        self.data = da

    def convert_positions(self, cartesian_positions):
        """
        Convert flattened points to geodetic coordinates
        :param cartesian_positions: Points to convert
        :return: Points in geodetic coordinates
        """

        # Get planet axes
        a, b, c = self.planet.get_axes()

        # Create a coordinate conversion object
        coordinate_conversions = pet.conversions.coordinateConversions(name="conversions", a=a, b=b, c=c)

        # Get the corresponding latitude and longitude with a triaxial ellipsoid conversion
        flattened_swath_coordinates = np.asarray(coordinate_conversions.geodetic(
            cartesian_coordinates=cartesian_positions))

        return flattened_swath_coordinates

    @staticmethod
    def create_line_of_sight(position, velocity, thetas, look_direction, force_zero_doppler):
        """ Find the LoS vector corresponding to a certain position and velocity
            vector, taking into consideration the look angle and squint

            :param position: antenna origins (satellite's position) [m].
            :param velocity: velocity of the antenna [m/s].
            :param thetas: vector angles [radians]
            :param look_direction: look direction (left or right)
            :param force_zero_doppler: if True, the LoS is forced to be perpendicular
             to the velocity vector (zero Doppler).
            returns: float 3-D array
            """

        n_la = thetas.size

        if velocity.ndim == 1:
            ax = 0
            n_v = 1

        else:
            ax = 1
            n_v = velocity.shape[0]

        # Calculate velocity and positions versor
        v_ver = velocity / np.linalg.norm(velocity, axis=ax).reshape((n_v, 1))
        r_ver = position / np.linalg.norm(position, axis=ax).reshape((n_v, 1))
        n_ver = np.cross(v_ver, r_ver)  # cross product of versors

        if force_zero_doppler:
            r_ver2 = np.cross(n_ver, v_ver)

        else:
            r_ver2 = r_ver

        if look_direction == "left":
            thetas = -thetas

        line_of_sight = (-np.cos(thetas.reshape(1, n_la, 1)) * r_ver2.reshape(n_v, 1, 3) +
                         np.sin(thetas.reshape(1, n_la, 1)) * n_ver.reshape(n_v, 1, 3))

        return line_of_sight

    def calculate_ground_swath(self, force_zero_doppler=True):
        """
        Vectorized calculation of the ground swath.
        :return: Nothing returned
        """

        print("Starting Swath Computation...", file=sys.stderr)

        # Get the time space
        time_space = np.arange(self.start_time, self.end_time + self.temporal_resolution,
                               self.temporal_resolution)

        # Get planet axes
        a, b, c = self.planet.get_axes()

        print("     Getting satellite positions...", file=sys.stderr)

        # Get positions and velocities of the satellite for the observation times
        satellite_positions, satellite_velocities = self.campaign.get_states(times=time_space)

        print("     Getting antenna beam sampling...", file=sys.stderr)

        # Get the number of rays to intercept with DSK in accordance with spatial resolution (rough approximation)
        ang_span = np.deg2rad(self.instrument.end_look_angle - self.instrument.start_look_angle)
        mean_height = np.mean(np.linalg.norm(satellite_positions, axis=1) - np.mean((a, b, c)))
        mean_range = mean_height / np.cos(np.mean(ang_span))
        num_theta = int((ang_span / self.spatial_resolution) * mean_range)

        # Get look angle array
        thetas = np.linspace(self.instrument.start_look_angle, self.instrument.end_look_angle, num_theta)
        thetas = np.deg2rad(thetas)

        print("     Getting lines of sight...", file=sys.stderr)
        line_of_sight_vectors = self.create_line_of_sight(position=satellite_positions,
                                                          velocity=satellite_velocities,
                                                          thetas=thetas, look_direction=self.instrument.look_direction,
                                                          force_zero_doppler=force_zero_doppler)

        print("     Getting intercept points and local look, incident angles...", file=sys.stderr)
        icp = np.zeros(line_of_sight_vectors.shape)
        incidence_angles = np.zeros((line_of_sight_vectors.shape[0], line_of_sight_vectors.shape[1]))
        non_projected_incidence_angles = np.zeros((line_of_sight_vectors.shape[0], line_of_sight_vectors.shape[1]))
        look_angles = np.zeros((line_of_sight_vectors.shape[0], line_of_sight_vectors.shape[1]))

        for pp in tqdm(range(satellite_positions.shape[0])):
            icp[pp, :, :], incidence_angles[pp, :], non_projected_incidence_angles[pp, :], look_angles[pp, ::] = (
                self.planet.get_surface_intersects_local_angles(
                    satellite_position=satellite_positions[pp, :],
                    raydirs=line_of_sight_vectors[pp, :, :]))

        cartesian_positions = np.reshape(icp, (icp.shape[0] * icp.shape[1], 3))
        flat_incidence_angles = np.reshape(incidence_angles, incidence_angles.size)
        flat_non_projected_incidence_angles =\
            np.reshape(non_projected_incidence_angles, non_projected_incidence_angles.size)
        flat_look_angles = np.reshape(look_angles, look_angles.size)
        flat_times = np.reshape(np.tile(time_space[:, np.newaxis], (1, icp.shape[1])), icp.shape[0] * icp.shape[1])

        print("     Coordinate conversion...", file=sys.stderr)
        swath_coordinates = self.convert_positions(cartesian_positions=cartesian_positions)

        self.create_data_array(cartesian_coordinates=cartesian_positions, geodetic_coordinates=swath_coordinates,
                               incidence_angles=flat_incidence_angles,
                               non_projected_incidence_angles= flat_non_projected_incidence_angles,
                               look_angles=flat_look_angles, times=flat_times)

        print("Swath calculation finished!", file=sys.stderr)

    def visualize_swath(self, projection, fig=None, globe=None, ax=None, return_fig=False):
        """
        Visualize the ground swath
        :param projection: Cartopy projection
        :param fig: matplotlib figure
        :param globe: cartopy globe
        :param ax: matplotlib ax
        :param return_fig: Whether to return fig, globe, ax
        :return: fig, ax, globe if return_fig is True
        """

        if fig is None:
            # Get the projection
            fig, ax, globe = projection.proj(planet=self.planet)

        # Access longitude and latitude coordinates
        longitudes = self.data["longitude"].values
        latitudes = self.data["latitude"].values

        # Plot points on the map
        ax.scatter(longitudes, latitudes, transform=ccrs.PlateCarree(globe=globe),
                   color='red', marker='o', s=0.1, alpha=0.25)

        # return fig, ax, globe if necessary
        if return_fig:
            return fig, ax, globe

        # Add labels and legend
        plt.title('Swath', pad=20)

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/swath_' + str(self.campaign.body_id) + '_' +
                    str(self.start_time) + '_' + str(self.end_time) + '.png', format='png', dpi=500)

# end of file
