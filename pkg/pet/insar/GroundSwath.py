#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import numpy as np
import cspyce as spice
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from .GroundTarget import GroundTarget


class GroundSwath:
    """
    An object containing the swath of a satellite pass. Contains the start, end, and time_interval of the satellite pass
    as well as given ground resolution to calculate vectors with. Also contains the satellite time vector and a
    2-D array of the azimuthal beams by row and GroundTargets populating the rows.
    """

    def __init__(self, start_time, end_time, planet, instrument, time_interval=10, ground_resolution=2000, copy=False):
        self.start_time = start_time
        self.end_time = end_time
        self.time_interval = time_interval
        self.ground_resolution = ground_resolution
        self.time_space = None
        self.swath_beams = []
        self.planet = planet
        self.instrument = instrument
        if copy is False:
            self.calculate_swath()

    def calculate_swath(self):
        """
        Takes an instrument and planet, and uses the defined swath object parameters to calculate the detected ground
        swath.
        """

        # Create list of times to observe at
        self.time_space = np.arange(self.instrument.convert_times(times=self.start_time),
                                    self.instrument.convert_times(times=self.end_time) +
                                    self.time_interval, self.time_interval)

        # Get planet axes
        a, b, c = self.planet.get_axes()

        # Get the number of rays to intercept with DSK in accordance with spatial resolution
        num_theta = int((2 * np.pi / self.ground_resolution) * np.average((a, b, c)))

        # Get positions and velocities of the satellite for the observation times
        satellite_positions, satellite_velocities = self.instrument.get_states(times=self.time_space[:])

        # Create the SPICE planes with the normal as the satellite velocities originating at planet center
        planes = spice.nvp2pl_vector(normal=satellite_velocities, point=[0, 0, 0])

        # Calculate SPICE ellipses that intercept the SPICE planes with planet triaxial ellipsoid.
        ellipses = spice.inedpl_vector(a=a, b=b, c=c, plane=planes)[0]

        # Get the generating vectors for the SPICE ellipses
        centers, smajors, sminors = spice.el2cgv_vector(ellipse=ellipses)

        # Get the angles to iterate through for the spanning vectors
        thetas = np.linspace(np.pi, -np.pi, num=num_theta, endpoint=False)[::-1]

        self.swath_beams = []

        # Iterate through each azimuthal beam
        for i in range(len(self.time_space)):

            # Create the vectors per temporal point
            vectors = [np.cos(theta) * smajors[i] + np.sin(theta) * sminors[i] for theta in thetas]

            # Get the intersects with the planet DSK
            intersects = self.planet.get_surface_intersects(vectors=vectors)

            # Make groundTarget objects with the intersects
            self.swath_beams.append([GroundTarget(x=intersect[0], y=intersect[1], z=intersect[2]) for
                                     intersect in intersects])

        # Calculate the relative positions of each GroundTarget for each beam
        relative_positions = self.point_positions_relative_to_satellite(satellite_positions=satellite_positions,
                                                                        satellite_velocities=satellite_velocities)

        # Calculate the look angle of each GroundTarget for each beam
        look_angles = self.get_angles_cartesian(times=self.time_space, satellite_positions=satellite_positions)

        # Calculate if the positions of the GroundTargets are past the limb point
        limb_checks = self.check_limbs(satellite_positions=satellite_positions)

        # Check the beams for off-limb, look_angles, and relative positions
        for i in range(len(self.swath_beams)):

            azimuthal_points = []
            for j in range(len(self.swath_beams[i])):

                # Get groundTarget
                ground = self.swath_beams[i][j]

                # Check if the look angle and relative position is correct
                if ((self.instrument.start_look_angle < look_angles[i][j] < self.instrument.end_look_angle)
                        and (relative_positions[i][j] == "right") and limb_checks[i][j] == 0):

                    # Append the groundTarget to azimuthal beam
                    azimuthal_points.append(ground)

            # Replace previous beam with all accepted points
            self.swath_beams[i] = azimuthal_points

    def check_limbs(self, satellite_positions):
        """
        Check if all the swath_beam rays from satellite to surface positions intersect the limb ellipse.
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
        for i in range(len(self.swath_beams)):

            # Get the surface positions per beam
            surface_positions = np.asarray([point.get_position() for point in self.swath_beams[i]])

            # Calculate vectors from satellite to point and the satellite's direction
            vectors_to_sat = satellite_positions[i] - surface_positions

            # Get intersect values
            intersect_values = spice.inrypl_vector(vertex=surface_positions, dir=vectors_to_sat, plane=planes[i])[0]

            # Append the intersects
            intersect_booleans_per_beam.append(intersect_values)

        # Return the intersections in a 2D array
        return intersect_booleans_per_beam

    def point_positions_relative_to_satellite(self, satellite_positions, satellite_velocities):
        """
        Determine whether a set of points on the surface of the Earth are to the left or right of the satellite's
        velocity vectors
        :param satellite_positions: Array of x, y, z, positions of the satellite [m]
        :param satellite_velocities: Array of x, y, z, velocity vectors of the satellite [m/s]
        :return: A list containing whether the points are to the left or right of their corresponding satellite velocity
        vector
        """

        # Iterate through the beams
        relative_positions_per_beam = []
        for i in range(len(self.swath_beams)):

            # Get surface positions per beam
            surface_positions = np.asarray([point.get_position() for point in self.swath_beams[i]])

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

    def get_angles_cartesian(self, times, satellite_positions):
        """
        Get the look angles between the GroundTargets and the satellite in degrees.
        :param times: times of observation
        :param satellite_positions: Array of x, y, z, positions of the satellite
        :return: Absolute look angle between the surface points and the corresponding satellite position in degrees
        """

        # Calculate the satellite intersects and heights with the shape
        intersects, satellite_heights = self.planet.get_sub_obs_points(times=times, instrument=self.instrument)

        # Calculate the distance between center and satellite intersects
        satellite_radii = np.asarray([np.linalg.norm(intersect - [0, 0, 0]) for intersect in intersects])

        # Iterate through the beams
        look_angles_per_beam = []
        for i in range(len(self.swath_beams)):

            # Get surface positions per beam
            surface_positions = np.asarray([point.get_position() for point in self.swath_beams[i]])

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

    def visualize(self, projection):
        """
        Visualize the ground swath
        :param projection: Cartopy projection
        """

        # Create empty list of positions
        positions = []

        # Get positions from swath beams
        for beam in self.swath_beams:
            for point in beam:
                positions.append(point.get_position())

        # Get planet axes
        a, b, c = self.planet.get_axes()

        # Create a coordinate conversion object
        convert = pet.ext.conversions(name="conversions", a=a, b=b, c=c)

        # Convert positions to geodetic coordinates
        geodetic_coordinates = convert.geodetic(cartesian_coordinates=positions)[:, :3]

        # Get the fig, ax, globe from the instrument orbit
        fig, ax, globe = self.instrument.plot_orbit(projection=projection,
                                                    start_time=self.start_time,
                                                    end_time=self.end_time,
                                                    time_interval=self.time_interval, return_fig=True)
        # Get the coordinates
        coordinates = [(long, lat, height) for lat, long, height in geodetic_coordinates]
        longitudes, latitudes, heights = zip(*coordinates)

        # Plot points on the map
        ax.scatter(longitudes, latitudes, transform=ccrs.PlateCarree(globe=globe),
                   color='red', marker='o', s=0.1, alpha=0.25)

        # Add labels and legend
        ax.set_title('Swath')

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/swath.png', format='png', dpi=500)

# end of file
