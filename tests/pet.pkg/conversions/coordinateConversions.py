#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet

# Create a conversions object
coordinate_conversions = pet.conversions.coordinateConversions(name="conversions", a="1*km", b="1*km", c="1*km")

# Make some cartesian coordinates
cartesian_coordinates = [10*pet.units.SI.m, 10*pet.units.SI.m, 10*pet.units.SI.m]

# Convert some coordinates to geodetic
geodetic_coordinates = coordinate_conversions.geodetic(cartesian_coordinates=cartesian_coordinates)

# Convert back to Cartesian
cartesian_coordinates_converted = coordinate_conversions.cartesian(geodetic_coordinates=geodetic_coordinates)

print(cartesian_coordinates[0], cartesian_coordinates[1], cartesian_coordinates[2],
      cartesian_coordinates_converted[0], cartesian_coordinates_converted[1],
      cartesian_coordinates_converted[2])

cartesian_coordinates = [[10*pet.units.SI.m, 10*pet.units.SI.m, 10*pet.units.SI.m],
                         [5*pet.units.SI.m, 5*pet.units.SI.m, 5*pet.units.SI.m]]

# Convert some coordinates to geodetic
geodetic_coordinates = coordinate_conversions.geodetic(cartesian_coordinates=cartesian_coordinates)

# Convert back to Cartesian
cartesian_coordinates_converted = coordinate_conversions.cartesian(geodetic_coordinates=geodetic_coordinates)

print([[x.value, y.value, z.value] for x, y, z in cartesian_coordinates],
      [[x.value, y.value, z.value] for x, y, z in cartesian_coordinates_converted])

# end of file
