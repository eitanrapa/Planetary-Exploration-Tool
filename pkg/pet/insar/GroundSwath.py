#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import numpy as np
import cspyce as spice
import rasterio
from matplotlib import pyplot
from .GroundTarget import GroundTarget


class GroundSwath(pet.component):
    """

    """

    start_time = pet.properties.str()
    start_time.doc = "start time of swath"

    end_time = pet.properties.str()
    end_time.doc = "end time of swath"

    time_interval = pet.properties.float()
    time_interval.doc = "time between observations (s)"

    ground_resolution = pet.properties.float()
    ground_resolution.doc = "spatial distance between observations (km)"

    def __init__(self, name, locator, implicit):
        super().__init__(name, locator, implicit)
        self.time_space = None
        self.swath_beams = []
        self.raster_fp = ""

    def make_raster(self):
        z = np.empty((len(self.swath_beams), len(self.swath_beams[0])))
        for i in range(len(self.swath_beams)):
            for j in range(len(self.swath_beams[0])):
                z[i][j] = self.swath_beams[i][j].z

        self.raster_fp = '/tmp/new.tif'
        new_dataset = rasterio.open(fp='/tmp/new.tif', mode='w', driver='GTiff',
                                    height=z.shape[0], width=z.shape[1],
                                    count=1, dtype=z.dtype, crs='+proj=cart')

        new_dataset.write(z, 1)
        new_dataset.close()

    def calculate_swath(self, instrument, planet):
        """

        """

        # Create list of times to observe at
        self.time_space = np.arange(instrument.convert_time(time=self.start_time),
                                    instrument.convert_time(time=self.end_time) +
                                    self.time_interval, self.time_interval)

        # Get planet axes
        a, b, c = planet.get_axes()

        # Get the number of rays to intercept with DSK in accordance with spatial resolution
        num_theta = int((2 * np.pi / self.ground_resolution) * np.average((a, b, c)))

        # Iterate through the times
        satellite_positions, satellite_velocities = instrument.get_states(target_body_id=planet.body_id,
                                                                          times=self.time_space[:],
                                                                          reference_body=planet.reference_id)

        # Create a SPICE plane with the normal as the satellite velocity originating at planet center
        planes = spice.nvp2pl_vector(normal=satellite_velocities, point=[0, 0, 0])

        # Calculate SPICE ellipse that intercept the SPICE plane with planet triaxial ellipsoid.
        ellipses = spice.inedpl_vector(a=a, b=b, c=c, plane=planes)[0]

        centers, smajors, sminors = spice.el2cgv_vector(ellipse=ellipses)

        # Get the angles to iterate through for the spanning vectors
        thetas = np.linspace(np.pi, -np.pi, num=num_theta, endpoint=False)[::-1]

        self.swath_beams = []

        # Iterate through each azimuthal beam
        for i in range(len(self.time_space)):
            # Create the vectors per temporal point
            vectors = [np.cos(theta) * smajors[i] + np.sin(theta) * sminors[i] for theta in thetas]

            # Get the intersects with the planet DSK
            intersects = planet.get_surface_intersects(vectors=vectors)

            # Make groundTarget objects with the intersects
            self.swath_beams.append([GroundTarget(name="{}".format(i), x=intersect[0], y=intersect[1], z=intersect[2])
                                     for intersect in intersects])

        # Calculate and check position of beam points with respect to satellite
        relative_positions = self.point_positions_relative_to_satellite(satellite_positions=satellite_positions,
                                                                        satellite_velocities=satellite_velocities)

        # Look angles
        look_angles = self.get_angles_cartesian(instrument=instrument, times=self.time_space, planet=planet,
                                                satellite_positions=satellite_positions)

        # Check the beams for duplicates, look_angles, and relative positions
        for i in range(len(self.swath_beams)):
            azimuthal_points = []
            for j in range(len(self.swath_beams[i])):

                # Get groundTarget
                groundpoint = self.swath_beams[i][j]
                if ((instrument.start_look_angle < look_angles[i][j] < instrument.end_look_angle)
                        and (relative_positions[i][j] == "right")):

                    # Append the groundTarget to azimuthal beam
                    azimuthal_points.append(groundpoint)

                    # Get vector from groundTarget to satellite
                    new_vector = groundpoint.get_position() - satellite_positions[i]

                    # For all points in azimuthal beam, check there are no duplicates
                    for k in range(len(azimuthal_points) - 1):

                        # Get iteration position
                        iteration_position = azimuthal_points[k].get_position()

                        # Get iteration vector
                        iteration_vector = iteration_position - satellite_positions[i]

                        # Check if the vectors are on the same path
                        if np.abs(np.dot(new_vector / np.linalg.norm(new_vector),
                                         iteration_vector / np.linalg.norm(iteration_vector))) == 1:

                            # Check if the iteration_vector is closer than the new vector
                            if (np.linalg.norm(iteration_position - satellite_positions[i]) <
                                    np.linalg.norm(groundpoint.get_position() - satellite_positions[i])):
                                # Remove newest groundTarget
                                azimuthal_points.pop()

                                # Stop loop
                                break

            # Add all accepted points in azimuthal beam to swath
            self.swath_beams[i] = azimuthal_points

    def point_positions_relative_to_satellite(self, satellite_positions, satellite_velocities):
        """
        Determine whether a point on the surface of the Earth is to the left or right of the satellite's velocity vector
        :param satellite_positions:
        :param satellite_velocities:
        :return: Whether the point is to the left or right of the satellite's velocity vector
        """

        relative_positions_per_beam = []
        for i in range(len(self.swath_beams)):
            surface_positions = np.asarray([point.get_position() for point in self.swath_beams[i]])

            # Calculate vectors from satellite to point and the satellite's direction
            vectors_to_point = surface_positions - satellite_positions[i]

            # Calculate the cross product of the vectors
            # noinspection PyUnreachableCode
            cross_products = np.cross(satellite_velocities[i], vectors_to_point)

            # Project cross_product onto radial vector
            radial_vector = satellite_positions[i] / np.linalg.norm(satellite_positions[i])
            projections = np.dot(cross_products, radial_vector)

            # Determine the position of the point relative to the satellite
            relative_positions = []
            for projection in projections:

                if projection > 0:
                    relative_positions.append("left")

                elif projection < 0:
                    relative_positions.append("right")

            relative_positions_per_beam.append(relative_positions)

        return relative_positions_per_beam

    def get_angles_cartesian(self, instrument, times, planet, satellite_positions):
        """
        Get the angle between the pixel and the satellite in degrees.
        :param instrument: instrument observer
        :param times: time of observation
        :param planet: target body
        :param satellite_positions:
        :return: Absolute look angle between the surface point and the satellite in degrees
        """

        # Calculate the satellite intersect with the shape
        intersects, satellite_heights = planet.get_sub_obs_points(times=times, instrument_id=instrument.body_id)

        # Calculate the distance between center and satellite intersect
        satellite_radii = np.asarray([np.linalg.norm(intersect - [0, 0, 0]) for intersect in intersects])

        look_angles_per_beam = []
        for i in range(len(self.swath_beams)):
            surface_positions = np.asarray([point.get_position() for point in self.swath_beams[i]])

            # Get the distance between the satellite and the surface point
            distances = np.asarray([np.linalg.norm(surface_position - satellite_positions[i])
                                    for surface_position in surface_positions])

            # Calculate cosine of the angle
            z_plus_re = satellite_heights[i] + satellite_radii[i]
            altitude = np.asarray([np.linalg.norm(surface_position - [0, 0, 0])
                                   for surface_position in surface_positions])
            cosine_look_angle = (distances ** 2 + z_plus_re ** 2 - altitude ** 2) / (2 * distances * z_plus_re)

            # Use arc cosine to get the angle in radians
            look_angles_radians = np.arccos(cosine_look_angle)

            # Convert radians to degrees
            look_angles_degrees = np.degrees(look_angles_radians)

            # Append to list
            look_angles_per_beam.append(look_angles_degrees)

        return look_angles_per_beam

    def visualize(self):
        """

        """

        self.make_raster()
        src = rasterio.open(fp=self.raster_fp)
        pyplot.imshow(src.read(1), cmap='pink')
        pyplot.show()

# end of file
