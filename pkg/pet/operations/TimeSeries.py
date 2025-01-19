#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import numpy as np
from scipy.interpolate import griddata


class TimeSeries(pet.component):
    """
    Class that is able to create time series information from loaded interferograms
    """

    def __init__(self, planet, conops, **kwargs):
        super().__init__(**kwargs)
        self.planet = planet
        self.conops = conops
        self.data = []

    def create_1d_time_series(self, interferograms, spatial_points):
        """
        Return the amplitudes, phases, and topographic errors of a set of points in 1 look direction
        :param interferograms: interferograms to analyze, in 1 look direction
        :param spatial_points: points to do the inverse problem on
        :return: amplitudes, phases, topographic errors
        """

        # Set up a matrix of observations and a problem G matrix
        d_vectors = {i: [] for i in range(len(spatial_points))}  # List to collect rows for each point
        g_matrices = {i: [] for i in range(len(spatial_points))}

        # Define the tidal cycle, and therefore the angular frequency
        tidal_cycle = self.planet.tidal_cycle
        omega = 2 * np.pi / tidal_cycle

        # Iterate through interferograms
        for interferogram in interferograms:

            # Load the data
            dataarray = interferogram.data

            # Extract data
            latitudes = dataarray["latitude"].values
            longitudes = dataarray["longitude"].values
            points = np.column_stack((latitudes, longitudes))

            # Interpolate each desired variable
            results = {var: [] for var in ["time1", "time2", "los_displacements", "psi", "x", "y", "z"]}
            results["valid_points"] = []  # To track valid points

            # Loop through target points
            for point in spatial_points:
                target_point = np.array([point[0], point[1]])

                # Interpolate each desired variable
                for var in ["time1", "time2", "los_displacements", "psi", "x", "y", "z"]:
                    values = dataarray[var].values
                    interpolated_value = griddata(points, values, target_point, method="linear")
                    results[var].append(interpolated_value)

                # Track the valid point
                results["valid_points"].append((point[0], point[1]))

            # Extract interpolated values
            time1_interp = (np.array(results["time1"])[0] % tidal_cycle)
            time2_interp = (np.array(results["time2"])[0] % tidal_cycle)
            los_displacement_interp = np.array(results["los_displacements"])[0]
            psis_interp = np.array(results["psi"])[0]

            # Append to d_vectors and G_matrices for each point using list comprehension
            [d_vectors[i].append(los_displacement_interp[i]) for i in range(len(spatial_points))]

            [g_matrices[i].append([np.cos(omega * time2_interp[i]) - np.cos(omega * time1_interp[i]),
                                   np.sin(omega * time2_interp[i]) - np.sin(omega * time1_interp[i]),
                                   psis_interp[i]]) for i in range(len(spatial_points))]

        # Convert lists to numpy arrays and reshape
        d_vectors = {i: np.array(d_vectors[i]) for i in d_vectors}  # Shape: (num_times, d_dim)
        g_matrices = {i: np.array(g_matrices[i]) for i in g_matrices}  # Shape: (num_times, G_dim)

        # Setup a data covariance matrix
        c_d = np.eye(len(interferograms))

        # Make a list of amplitudes, phases, topography errors
        amplitudes = []
        phis = []
        zs = []

        for i in range(len(spatial_points)):
            # Set a small regularization parameter, equivalent to Cm
            lambda_reg = 1e-5

            # Do the inverse problem
            m = np.linalg.inv(g_matrices[i].T @ np.linalg.inv(c_d) @ g_matrices[i] +
                              lambda_reg * np.eye(7)) @ (g_matrices[i].T @ np.linalg.inv(c_d) @ d_vectors[i])

            # Retrieve the parameters
            c1, c2, c3, s1, s2, s3, z = m

            # Calculate amplitudes and phases
            amplitudes.append([np.sqrt(c1 ** 2 + s1 ** 2), np.sqrt(c2 ** 2 + s2 ** 2), np.sqrt(c3 ** 2 + s3 ** 2)])
            phis.append([np.arctan(c1 / s1), np.arctan(c2 / s2), np.arctan(c3 / s3)])
            zs.append(z)

        # Return the results
        return amplitudes, phis, zs

    def create_3d_time_series(self, interferograms, spatial_points):
        """
        Return the amplitudes, phases, and topographic errors of a set of points in 3 or more look directions
        :param interferograms: interferograms to analyze
        :param spatial_points: points to do the inverse problem on
        :return: amplitudes, phases, topographic errors
        """

        # Set up a matrix of observations and a problem G matrix
        d_vectors = {i: [] for i in range(len(spatial_points))}  # List to collect rows for each point
        g_matrices = {i: [] for i in range(len(spatial_points))}

        # Define the tidal cycle, and therefore the angular frequency
        tidal_cycle = self.planet.tidal_cycle
        omega = 2 * np.pi / tidal_cycle

        # Create a conversion object
        vector_conversions = pet.conversions.vectorConversions()

        # Iterate through the interferograms
        for interferogram in interferograms:

            # Load the data
            dataarray = interferogram.data

            # Extract data
            latitudes = dataarray["latitude"].values
            longitudes = dataarray["longitude"].values
            points = np.column_stack((latitudes, longitudes))

            # Interpolate each desired variable
            results = {var: [] for var in ["time1", "time2", "los_displacements", "sat_pos_time", "psi", "x", "y", "z"]}
            results["valid_points"] = []  # To track valid points

            # Loop through target points
            for point in spatial_points:
                target_point = np.array([point[0], point[1]])

                # Interpolate each desired variable
                for var in ["time1", "time2", "los_displacements", "sat_pos_time", "psi", "x", "y", "z"]:
                    values = dataarray[var].values
                    interpolated_value = griddata(points=points, values=values, xi=target_point, method="nearest")
                    results[var].extend(interpolated_value)

                # Track the valid point
                results["valid_points"].append((point[0], point[1]))

            # Extract interpolated values
            sat_times = np.asarray(results["sat_pos_time"])
            time1_interps = (np.asarray(results["time1"]) % tidal_cycle)
            time2_interps = (np.asarray(results["time2"]) % tidal_cycle)
            x_interps = np.asarray(results["x"])
            y_interps = np.asarray(results["y"])
            z_interps = np.asarray(results["z"])
            los_displacement_interps = np.asarray(results["los_displacements"])
            psis_interps = np.asarray(results["psi"])

            # Get the satellite and surface positions
            satellite_positions, satellite_velocities = self.conops.get_states(times=sat_times)
            satellite_positions = satellite_positions.T
            surface_positions = np.asarray([x_interps, y_interps, z_interps])

            # Calculate vectors from satellite to point and the satellite's direction
            vectors_to_sat = (satellite_positions - surface_positions)
            vectors_to_sat = vectors_to_sat.T
            vectors_to_sat = [vector / np.linalg.norm(vector) for vector in vectors_to_sat]

            # Append to d_vectors and G_matrices for each point using list comprehension
            [d_vectors[i].append(los_displacement_interps[i]) for i in range(len(spatial_points))]

            # Construct the G matrix, passing through vector conversion to ENU
            [g_matrices[i].append([
                *[arr.squeeze() for arr
                  in vector_conversions.cartesian_to_enu_vector(
                        uvw_vectors=
                        [vectors_to_sat[i][0] * (np.cos(omega * time2_interps[i]) - np.cos(omega * time1_interps[i])),
                         vectors_to_sat[i][1] * (np.cos(omega * time2_interps[i]) - np.cos(omega * time1_interps[i])),
                         vectors_to_sat[i][2] * (np.cos(omega * time2_interps[i]) - np.cos(omega * time1_interps[i]))],
                        latitudes=[spatial_points[i, 0]], longitudes=[spatial_points[i, 1]]
                    )],  # First 3 elements: vector0
                *[arr.squeeze() for arr
                  in vector_conversions.cartesian_to_enu_vector(
                        uvw_vectors=
                        [vectors_to_sat[i][0] * (np.sin(omega * time2_interps[i]) - np.sin(omega * time1_interps[i])),
                         vectors_to_sat[i][1] * (np.sin(omega * time2_interps[i]) - np.sin(omega * time1_interps[i])),
                         vectors_to_sat[i][2] * (np.sin(omega * time2_interps[i]) - np.sin(omega * time1_interps[i]))],
                        latitudes=[spatial_points[i, 0]], longitudes=[spatial_points[i, 1]]
                    )],  # Next 3 elements: vector1
                psis_interps[i],  # Final element
            ]
            ) for i in range(len(spatial_points))]

        # Convert lists to numpy arrays and reshape
        d_vectors = {i: np.array(d_vectors[i]) for i in d_vectors}  # Shape: (num_times, d_dim)
        g_matrices = {i: np.array(g_matrices[i]) for i in g_matrices}  # Shape: (num_times, G_dim)

        # Setup a data covariance matrix
        c_d = np.eye(len(interferograms))

        # Make a list of amplitudes, phases, topography errors
        amplitudes = []
        phis = []
        zs = []

        for i in range(len(spatial_points)):
            # Set a small regularization parameter, equivalent to Cm
            lambda_reg = 1e-5

            # Do the inverse problem
            m = np.linalg.inv(g_matrices[i].T @ np.linalg.inv(c_d) @ g_matrices[i] +
                              lambda_reg * np.eye(7)) @ (g_matrices[i].T @ np.linalg.inv(c_d) @ d_vectors[i])

            # Retrieve the parameters
            c1, c2, c3, s1, s2, s3, z = m

            # Calculate amplitudes and phases
            amplitudes.append([np.sqrt(c1 ** 2 + s1 ** 2), np.sqrt(c2 ** 2 + s2 ** 2), np.sqrt(c3 ** 2 + s3 ** 2)])
            phis.append([np.arctan2(c1, s1), np.arctan2(c2, s2), np.arctan2(c3, s3)])
            zs.append(z)

        # Return the results
        return amplitudes, phis, zs

# end of file
