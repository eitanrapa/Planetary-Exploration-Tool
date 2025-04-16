#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import pet

# Create vector conversion object
vector_conversions = pet.conversions.vectorConversions()

# Define ENU vector and local coordinates
enu_vectors = [1, 2, 3]
latitudes = 40
longitudes = 70

# Convert vectors
uvw_vectors = vector_conversions.enu_to_cartesian_vector(enu_vectors=enu_vectors, latitudes=latitudes,
                                                         longitudes=longitudes)

# Convert back
enu_vectors_converted = vector_conversions.cartesian_to_enu_vector(uvw_vectors=uvw_vectors, latitudes=latitudes,
                                                                   longitudes=longitudes)

print(enu_vectors)
print(enu_vectors_converted)

# Define ENU vectors and local coordinates
enu_vectors = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
latitudes = [40, 50, 60]
longitudes = [70, 80, 90]

# Convert vectors
uvw_vectors = vector_conversions.enu_to_cartesian_vector(enu_vectors=enu_vectors, latitudes=latitudes,
                                                         longitudes=longitudes)

# Convert back
enu_vectors_converted = vector_conversions.cartesian_to_enu_vector(uvw_vectors=uvw_vectors, latitudes=latitudes,
                                                                   longitudes=longitudes)

print(enu_vectors)
print(enu_vectors_converted)

# end of file
