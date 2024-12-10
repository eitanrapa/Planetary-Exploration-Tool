#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet

# Create a conversions object
coordinate_conversions = pet.conversions.coordinateConversions(name="conversions", a=1, b=1, c=1)

# Make some cartesian coordinates
cartesian_coordinates = [10, 10, 10]

# Convert some coordinates to geodetic
geodetic_coordinates = coordinate_conversions.geodetic(cartesian_coordinates=cartesian_coordinates)

# Convert back to Cartesian
geodetic_coordinates_converted = coordinate_conversions.cartesian(geodetic_coordinates=geodetic_coordinates)

print(cartesian_coordinates, geodetic_coordinates_converted)

cartesian_coordinates = [[10, 10, 10], [5, 5, 5]]

# Convert some coordinates to geodetic
geodetic_coordinates = coordinate_conversions.geodetic(cartesian_coordinates=cartesian_coordinates)

# Convert back to Cartesian
cartesian_coordinates_converted = coordinate_conversions.cartesian(geodetic_coordinates=geodetic_coordinates)

print(cartesian_coordinates, cartesian_coordinates_converted)

# end of file
