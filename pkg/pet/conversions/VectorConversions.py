# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# california institute of technology
# (c) 2023-2024 all rights reserved
#

import pet
import numpy as np


class VectorConversions(pet.component, family="pet.conversions.vectorConversions",
                        implements=pet.protocols.conversion):
    """
    Class that encapsulates conversions between ENU and Cartesian vectors
    """

    def _enu_to_cartesian_vector(self, enu_vectors, latitudes, longitudes):
        """
        Function that converts between local east, north, up vectors and ellipsoidal center x, y, z vectors.
        Also needs local coordinates.
        :param enu_vectors: local east, north, up vectors
        :param latitudes: local latitudes of points
        :param longitudes: local longitudes of points
        :return: u, v, w vectors
        """

        # Convert longitude and latitude to radians
        longitudes = np.deg2rad(longitudes)
        latitudes = np.deg2rad(latitudes)

        # Get easts, norths, ups
        east = enu_vectors[:, 0]
        north = enu_vectors[:, 1]
        up = enu_vectors[:, 2]

        # Calculate vector transformation
        t = np.cos(latitudes) * up - np.sin(latitudes) * north
        u = np.cos(longitudes) * t - np.sin(longitudes) * east
        v = np.sin(longitudes) * t + np.cos(longitudes) * east
        w = np.sin(latitudes) * up + np.cos(latitudes) * north

        # Return the vector
        return np.asarray([u, v, w]).T

    def enu_to_cartesian_vector(self, enu_vectors, latitudes, longitudes):
        """
        Helper function for _enu_to_cartesian_vector
        """

        # Make sure inputs are numpy arrays
        enu_vectors = np.asarray(enu_vectors)
        latitudes = np.asarray(latitudes)
        longitudes = np.asarray(longitudes)

        # Make sure dimensions is 2
        single_return = False
        if enu_vectors.ndim == 1:
            single_return = True
            enu_vectors = np.asarray([enu_vectors])

        # Make sure dimensions is 1
        if latitudes.ndim == 0:
            latitudes = np.asarray([latitudes])

        if longitudes.ndim == 0:
            longitudes = np.asarray([longitudes])

        # Call function
        cartesian_vectors = self._enu_to_cartesian_vector(enu_vectors=enu_vectors, latitudes=latitudes,
                                                          longitudes=longitudes)

        # Strip array if ndim is 1
        if single_return:
            return cartesian_vectors[0]
        else:
            return cartesian_vectors

    def _cartesian_to_enu_vector(self, uvw_vectors, latitudes, longitudes):
        """
        Convert Cartesian coordinates (u, v, w) back to ENU coordinates (east, north, up).
        Also needs local coordinates.
        :param uvw_vectors: local cartesian vectors
        :param latitudes: local latitudes of points
        :param longitudes: local longitudes of points
        :return: ENU vectors
        """

        # Convert degrees to radians
        longitudes = np.deg2rad(longitudes)
        latitudes = np.deg2rad(latitudes)

        # Get us, vs, ws
        u = uvw_vectors[:, 0]
        v = uvw_vectors[:, 1]
        w = uvw_vectors[:, 2]

        # Reverse the transformations step by step
        t = np.cos(longitudes) * u + np.sin(longitudes) * v
        east = -np.sin(longitudes) * u + np.cos(longitudes) * v
        up = np.cos(latitudes) * t + np.sin(latitudes) * w
        north = -np.sin(latitudes) * t + np.cos(latitudes) * w

        return np.asarray([east, north, up]).T

    def cartesian_to_enu_vector(self, uvw_vectors, latitudes, longitudes):
        """
        Helper function for _enu_to_cartesian_vector
        """

        # Make sure inputs are numpy arrays
        uvw_vectors = np.asarray(uvw_vectors)
        latitudes = np.asarray(latitudes)
        longitudes = np.asarray(longitudes)

        # Make sure dimensions is 2
        single_return = False
        if uvw_vectors.ndim == 1:
            single_return = True
            uvw_vectors = np.asarray([uvw_vectors])

        # Make sure dimensions is 1
        if latitudes.ndim == 0:
            latitudes = np.asarray([latitudes])

        if longitudes.ndim == 0:
            longitudes = np.asarray([longitudes])

        # Call function
        enu_vectors = self._cartesian_to_enu_vector(uvw_vectors=uvw_vectors, latitudes=latitudes, longitudes=longitudes)

        # Strip array if ndim is 1
        if single_return:
            return enu_vectors[0]
        else:
            return enu_vectors

    def enu_to_cartesian_cubes(self, east_cube, north_cube, up_cube, latitudes, longitudes):
        """
        Function that converts between local east, north, up vectors and ellipsoidal center x, y, z vectors.
        Also needs local coordinates.
        :param east_cube: local east pointing cube
        :param north_cube: local north pointing cube
        :param up_cube: local up pointing cube
        :param latitudes: local latitude of points
        :param longitudes: local longitude of points
        :return: u, v, w vector
        """

        # Convert longitude and latitude to radians
        longitude = np.deg2rad(longitudes)
        latitude = np.deg2rad(latitudes)

        # Reshape longitudes and latitudes for broadcasting
        longitude = longitude[:, np.newaxis, np.newaxis]
        latitude = latitude[np.newaxis, :, np.newaxis]

        # Calculate vector transformation
        t = np.cos(latitude) * up_cube - np.sin(latitude) * north_cube
        w = np.sin(latitude) * up_cube + np.cos(latitude) * north_cube
        u = np.cos(longitude) * t - np.sin(longitude) * east_cube
        v = np.sin(longitude) * t + np.cos(longitude) * east_cube

        # Return the vector
        return np.asarray([u, v, w])

# end of file
