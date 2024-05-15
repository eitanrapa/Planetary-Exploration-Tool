#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import numpy as np


class GroundTarget(pet.component):
    """
    Class that represents a target point on the ground
    """

    x = pet.properties.float()
    x.doc = "x coordinate"

    y = pet.properties.float()
    y.doc = "y coordinate"

    z = pet.properties.float()
    z.doc = "z coordinate"

    def get_position(self):
        return self.x, self.y, self.z

    def point_position_relative_to_satellite(self, instrument, time, planet):
        """
        Determine whether a point on the surface of the Earth is to the left or right of the satellite's velocity vector
        :param instrument: instrument observer
        :param time: time of observation
        :param planet: target body
        :return: Whether the point is to the left or right of the satellite's velocity vector
        """

        # Pack the coordinates into vector
        surface_position = np.asarray((self.x, self.y, self.z))

        # Get the satellite position
        satellite_position = instrument.get_state(target_body_id=planet.body_id, time=time,
                                                  reference_body=planet.reference_id)[0]
        satellite_velocity = instrument.get_state(target_body_id=planet.body_id, time=time,
                                                  reference_body=planet.reference_id)[1]

        # Calculate vectors from satellite to point and the satellite's direction
        vector_to_point = surface_position - satellite_position

        # Calculate the cross product of the vectors
        # noinspection PyUnreachableCode
        cross_product = np.cross(satellite_velocity, vector_to_point)

        # Project cross_product onto radial vector
        radial_vector = satellite_position / np.linalg.norm(satellite_position)
        projection = np.dot(cross_product, radial_vector)

        # Determine the position of the point relative to the satellite
        if projection > 0:

            return "left"

        elif projection < 0:

            return "right"

    def get_angles_cartesian(self, instrument, time, planet):
        """
        Get the angle between the pixel and the satellite in degrees.
        :param instrument: instrument observer
        :param time: time of observation
        :param planet: target body
        :return: Absolute look angle between the surface point and the satellite in degrees
        """

        # Pack the coordinates into vector
        surface_position = np.asarray((self.x, self.y, self.z))

        # Get the satellite position
        satellite_position = instrument.get_state(target_body_id=planet.body_id, time=time,
                                                  reference_body=planet.reference_id)[0]

        # Get the distance between the satellite and the surface point
        distance = np.linalg.norm(surface_position - satellite_position)

        # Calculate the satellite intersect with the shape
        intersect, satellite_height = planet.get_sub_obs_point(time=time, instrument_id=instrument.body_id)

        # Calculate the distance between center and satellite intersect
        satellite_radius = np.linalg.norm(intersect - [0, 0, 0])

        # Calculate cosine of the angle
        z_plus_re = satellite_height + satellite_radius
        altitude = np.linalg.norm(surface_position - [0, 0, 0])
        cosine_look_angle = (distance ** 2 + z_plus_re ** 2 - altitude ** 2) / (2 * distance * z_plus_re)

        # Use arc cosine to get the angle in radians
        look_angle_radians = np.arccos(cosine_look_angle)

        # Convert radians to degrees
        look_angle_degrees = np.degrees(look_angle_radians)

        return look_angle_degrees

# end of file
