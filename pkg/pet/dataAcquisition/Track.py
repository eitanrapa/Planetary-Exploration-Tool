#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

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
    as well as given ground resolution to calculate vectors with. Also contains the satellite time vector and a
    2-D array of the azimuthal beams by row and GroundTargets populating the rows.
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

    def create_data_array(self, cartesian_coordinates, geodetic_coordinates, look_angles, times):
        """
        Create a xarray with the input data
        :param cartesian_coordinates: x, y, z coordinates to be saved
        :param geodetic_coordinates: lat, long, height coordinates to be saved
        :param look_angles: Calculated look angles from satellite
        :param times: Times of observation for points
        :return: Nothing returned
        """

        # Create the xarray Dataset
        da = xr.DataArray(
            data=look_angles,
            dims=["points"],
            coords={
                "sat_pos_time": ("points", times),
                "time": ("points", times),
                "x": ("points", np.asarray([point[0] for point in cartesian_coordinates])),
                "y": ("points", np.asarray([point[1] for point in cartesian_coordinates])),
                "z": ("points", np.asarray([point[2] for point in cartesian_coordinates])),
                "latitude": ("points", np.asarray([point[0] for point in geodetic_coordinates])),
                "longitude": ("points", np.asarray([point[1] for point in geodetic_coordinates])),
                "height": ("points", np.asarray([point[2] for point in geodetic_coordinates]))},
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

    def modify_time(self, time):
        """
        Modify the time of the track
        :param time: Time to add to the track
        :return: Nothing returned
        """

        self.start_time = self.start_time + time
        self.end_time = self.end_time + time
        self.data["time"].values = self.data["time"].values + time
        self.data.attrs["start_time"] = self.data.attrs["start_time"] + time
        self.data.attrs["end_time"] = self.data.attrs["end_time"] + time

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

    def calculate_ground_swath(self):
        """
        Vectorized calculation of the ground swath.
        :return: Nothing returned
        """

        # Get the time space
        time_space = np.arange(self.start_time, self.end_time + self.temporal_resolution,
                               self.temporal_resolution)

        # Get planet axes
        a, b, c = self.planet.get_axes()

        # Get the number of rays to intercept with DSK in accordance with spatial resolution
        num_theta = int((2 * np.pi / self.spatial_resolution) * np.average((a, b, c)))

        print("Getting satellite positions...", file=sys.stderr)

        # Get positions and velocities of the satellite for the observation times
        satellite_positions, satellite_velocities = self.campaign.get_states(times=time_space)

        print("Calculating SPICE planes...", file=sys.stderr)

        # Create the SPICE planes with the normal as the satellite velocities originating at planet center
        planes = spice.nvp2pl_vector(normal=satellite_velocities, point=[0, 0, 0])

        print("Calculating SPICE ellipses...", file=sys.stderr)

        # Calculate SPICE ellipses that intercept the SPICE planes with planet triaxial ellipsoid.
        ellipses = spice.inedpl_vector(a=a, b=b, c=c, plane=planes)[0]

        print("Calculating generating vectors...", file=sys.stderr)

        # Get the generating vectors for the SPICE ellipses
        centers, smajors, sminors = spice.el2cgv_vector(ellipse=ellipses)

        # Get the angles to iterate through for the spanning vectors
        thetas = np.linspace(np.pi, -np.pi, num=num_theta, endpoint=False)[::-1]

        swath_beams = []

        # Iterate through each azimuthal beam
        for i in tqdm(range(len(time_space)), desc="Calculating planet surface intersects..."):

            # Create the vectors per temporal point
            vectors = [np.cos(theta) * smajors[i] + np.sin(theta) * sminors[i] for theta in thetas]

            # Get the intersects with the planet DSK
            intersects = self.planet.get_surface_intersects(vectors=vectors)

            # Make groundTarget objects with the intersects
            swath_beams.append([[intersect[0], intersect[1], intersect[2]] for intersect in intersects])

        # Calculate the relative positions of each GroundTarget for each beam
        relative_positions = self.point_positions_relative_to_satellite(swath_beams=swath_beams,
                                                                        satellite_positions=satellite_positions,
                                                                        satellite_velocities=satellite_velocities)

        # Calculate the look angle of each GroundTarget for each beam
        look_angles = self.get_angles_cartesian(swath_beams=swath_beams,
                                                times=time_space, satellite_positions=satellite_positions)

        # Calculate if the positions of the GroundTargets are past the limb point
        limb_checks = self.check_limbs(swath_beams=swath_beams, satellite_positions=satellite_positions)

        # Check the beams for off-limb, look_angles, and relative positions
        for i in range(len(swath_beams)):

            indices = []
            for j in range(len(swath_beams[i])):

                # Check if the look angle and relative position is correct
                if ((self.instrument.start_look_angle < look_angles[i][j] < self.instrument.end_look_angle)
                        and (relative_positions[i][j] == self.instrument.look_direction) and limb_checks[i][j] == 0):
                    # Append the groundTarget to azimuthal beam
                    indices.append(j)

            # Replace previous beam with all accepted points
            swath_beams[i] = [swath_beams[i][index] for index in indices]
            look_angles[i] = [look_angles[i][index] for index in indices]

        cartesian_positions = np.asanyarray([coord for sublist in swath_beams for coord in sublist])
        flat_angles = np.asarray([angle for sublist in look_angles for angle in sublist])
        times = [[time] * len(beams) for time, beams in zip(time_space, swath_beams)]
        flat_times = np.asarray([time for sublist in times for time in sublist])
        swath_coordinates = self.convert_positions(cartesian_positions=cartesian_positions)

        self.create_data_array(cartesian_coordinates=cartesian_positions, geodetic_coordinates=swath_coordinates,
                               look_angles=flat_angles, times=flat_times)

    def check_limbs(self, swath_beams, satellite_positions):
        """
        Check if all the swath_beam rays from satellite to surface positions intersect the limb ellipse.
        :param swath_beams: Azimuthal beams of the swath
        :param satellite_positions: Array of x, y, z, positions of the satellite [m]
        :return: A list of values whether the points on the ground are before or after the limb point.
        """

        # Get planet axes
        a, b, c = self.planet.get_axes()

        # Get the limb ellipses from satellite positions
        limb_ellipses = spice.edlimb_vector(a=a, b=b, c=c, viewpt=satellite_positions)

        # Get the generating vectors for the SPICE ellipses
        centers, smajors, sminors = spice.el2cgv_vector(ellipse=limb_ellipses)

        # Get planes generated by spanning vectors
        planes = spice.psv2pl_vector(point=centers, span1=smajors, span2=sminors)

        # Iterate through the beams
        intersect_booleans_per_beam = []
        for i in tqdm(range(len(swath_beams)), desc="Checking limb points..."):

            # Get the surface positions per beam
            surface_positions = np.asarray([(point[0], point[1], point[2]) for point in swath_beams[i]])

            # Calculate vectors from satellite to point and the satellite's direction
            vectors_to_sat = satellite_positions[i] - surface_positions

            # Get intersect values
            intersect_values = spice.inrypl_vector(vertex=surface_positions, dir=vectors_to_sat, plane=planes[i])[0]

            # Append the intersects
            intersect_booleans_per_beam.append(intersect_values)

        # Return the intersections in a 2D array
        return intersect_booleans_per_beam

    def point_positions_relative_to_satellite(self, swath_beams, satellite_positions, satellite_velocities):
        """
        Determine whether a set of points on the surface of the Earth are to the left or right of the satellite's
        velocity vectors
        :param swath_beams: Azimuthal beams of the swath
        :param satellite_positions: Array of x, y, z, positions of the satellite [m]
        :param satellite_velocities: Array of x, y, z, velocity vectors of the satellite [m/s]
        :return: A list containing whether the points are to the left or right of their corresponding satellite velocity
        vector
        """

        # Iterate through the beams
        relative_positions_per_beam = []
        for i in tqdm(range(len(swath_beams)), desc="Calculating relative positions..."):

            # Get surface positions per beam
            surface_positions = np.asarray([(point[0], point[1], point[2]) for point in swath_beams[i]])

            # Calculate vectors from satellite to point and the satellite's direction
            vectors_to_point = surface_positions - satellite_positions[i]

            # Calculate the cross products of the vectors
            # noinspection PyUnreachableCode
            cross_products = np.cross(satellite_velocities[i], vectors_to_point)

            # Project cross_products onto corresponding radial vector
            radial_vector = satellite_positions[i] / np.linalg.norm(satellite_positions[i])
            projections = np.dot(cross_products, radial_vector)

            # Determine the position of the point relative to the satellite
            relative_positions = []
            for projection in projections:

                if projection > 0:
                    relative_positions.append("left")

                elif projection < 0:
                    relative_positions.append("right")

            # Append the relative position
            relative_positions_per_beam.append(relative_positions)

        # Return the relative position 2-D array
        return relative_positions_per_beam

    def get_angles_cartesian(self, swath_beams, times, satellite_positions):
        """
        Get the look angles between the GroundTargets and the satellite in degrees.
        :param swath_beams: Azimuthal beams of the swath
        :param times: Times of observation
        :param satellite_positions: Array of x, y, z, positions of the satellite
        :return: Absolute look angle between the surface points and the corresponding satellite position in degrees
        """

        # Calculate the satellite intersects and heights with the shape
        intersects, satellite_heights = self.planet.get_sub_obs_points(times=times, campaign=self.campaign)

        # Calculate the distance between center and satellite intersects
        satellite_radii = np.asarray([np.linalg.norm(intersect - [0, 0, 0]) for intersect in intersects])

        # Iterate through the beams
        look_angles_per_beam = []
        for i in tqdm(range(len(swath_beams)), desc="Calculating look angles..."):

            # Get surface positions per beam
            surface_positions = np.asarray([(point[0], point[1], point[2]) for point in swath_beams[i]])

            # Get the distance between the satellite and the surface points
            distances = np.asarray([np.linalg.norm(surface_position - satellite_positions[i])
                                    for surface_position in surface_positions])

            # Calculate cosines of the angles
            z_plus_re = satellite_heights[i] + satellite_radii[i]
            altitude = np.asarray([np.linalg.norm(surface_position - [0, 0, 0])
                                   for surface_position in surface_positions])
            cosine_look_angle = (distances ** 2 + z_plus_re ** 2 - altitude ** 2) / (2 * distances * z_plus_re)

            # Use arc cosine to get the angles in radians
            look_angles_radians = np.arccos(cosine_look_angle)

            # Convert radians to degrees
            look_angles_degrees = np.degrees(look_angles_radians)

            # Yield output
            look_angles_per_beam.append(look_angles_degrees)

        # Return the look angles 2-D array
        return look_angles_per_beam

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
