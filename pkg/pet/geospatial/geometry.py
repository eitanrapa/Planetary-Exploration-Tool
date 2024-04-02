#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import numpy as np


def get_distance_cartesian(x1, y1, z1, x2, y2, z2):
    """
    Get the distance between two points in Cartesian coordinates.
    :param x1: x coordinate of the first point
    :param y1: y coordinate of the first point
    :param z1: z coordinate of the first point
    :param x2: x coordinate of the second point
    :param y2: y coordinate of the second point
    :param z2: z coordinate of the second point
    :return: Cartesian distance between the two points
    """

    # Calculate the distance between the two points using the Pythagorean theorem
    distance = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

    return distance


def point_position_relative_to_satellite(satellite_position, surface_position, satellite_velocity):
    """
    Determine whether a point on the surface of the Earth is to the left or right of the satellite's velocity vector.
    :param satellite_position: Vector containing the satellite's position in Cartesian coordinates
    :param surface_position: Vector containing the point's position in Cartesian coordinates
    :param satellite_velocity: Vector containing the satellite's velocity in Cartesian coordinates
    :return: Whether the point is to the left or right of the satellite's velocity vector
    """

    # Calculate vectors from satellite to point and the satellite's direction
    vector_to_point = surface_position - satellite_position

    # Calculate the cross product of the vectors
    # noinspection PyUnreachableCode
    cross_product = np.cross(satellite_velocity, vector_to_point)

    # Project cross_product onto radial vector
    radial_vector = sat_position / np.linalg.norm(sat_position)
    projection = np.dot(cross_product, radial_vector)

    # Determine the position of the point relative to the satellite
    if projection > 0:

        return "left"

    elif projection < 0:

        return "right"


def get_look_angle(surface_position, satellite_position, shape):
    """
    Get the angle between the pixel and the satellite in degrees.
    :param satellite_position: Vector containing the satellite's position in Cartesian coordinates
    :param surface_position: Vector containing the point's position in Cartesian coordinates
    :param shape: Shape class to get the intersect function
    :return: Absolute angle between the pixel and the satellite in degrees
    """

    # Unpack the positions
    surface_x, surface_y, surface_z = surface_position
    satellite_x, satellite_y, satellite_z = satellite_position

    # Get the distance between the satellite and the surface point
    distance = get_distance_cartesian(surface_x, surface_y, surface_z, satellite_x, satellite_y, satellite_z)

    # Calculate the satellite interesect with the shape
    intersect_x, intersect_y, intersect_z = shape.intersect(satellite_position)

    # Calculate the distance between center and satellite intersect
    satellite_radius = get_distance_cartesian(intersect_x, intersect_y, intersect_z, 0, 0, 0)

    # Calculate the distance between the satellite intersect and the satellite
    satellite_height = get_distance_cartesian(intersect_x, intersect_y, intersect_z, satellite_x,
                                              satellite_y, satellite_z)

    # Calculate cosine of the angle
    z_plus_re = satellite_height + satellite_radius
    altitude = get_distance_cartesian(surface_x, surface_y, surface_z, 0, 0, 0)
    cosine_look_angle = (distance ** 2 + z_plus_re ** 2 - altitude ** 2) / (2 * distance * z_plus_re)

    # Use arc cosine to get the angle in radians
    look_angle_radians = np.arccos(cosine_look_angle)

    # Convert radians to degrees
    look_angle_degrees = np.degrees(look_angle_radians)

    return look_angle_degrees

# end of file
