#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved
import pet
import numpy as np
from scipy.interpolate import griddata
from tqdm import tqdm
import alphashape
from shapely.geometry import Point
import xarray as xr


class TimeSeries_1D(pet.component):
    """
    Class that is able to create time series information from loaded interferograms
    """

    def __init__(self, planet, conops, instrument, **kwargs):
        super().__init__(**kwargs)
        self.planet = planet
        self.conops = conops
        self.instrument = instrument
        self.data = []

    @classmethod
    def load_from_file(cls, planet, conops, instrument, filename):
        """

        """

        # Open the HDF5 file in read mode
        data = xr.open_dataset(filename_or_obj=filename)

        obj = cls(planet=planet, conops=conops, instrument=instrument)
        obj.data = data  # Restore computed result
        return obj

    @classmethod
    def from_files(cls, planet, conops, instrument, file_list):
        """

        """

        # Load all the files
        return [cls.from_file(planet=planet, conops=conops, instrument=instrument, filename=file) for file in file_list]

    def save(self, file_name):
        """
        Save the track to an HDF5 file
        :return: Nothing returned
        """

        # Open HDF5 file
        self.data.to_netcdf(file_name)

    def create_data_array(self, cartesian_coordinates, geodetic_coordinates, amplitudes, amplitude_uncertainties,
                          phases, phase_uncertainties, topography_corrections):
        """
        Create a xarray with the input data
        :param cartesian_coordinates: x, y, z coordinates to be saved
        :param geodetic_coordinates: lat, long, height coordinates to be saved
        :param amplitudes: amplitudes to be saved
        :param amplitude_uncertainties: amplitude uncertainties to be saved
        :param phases: phases to be saved
        :param phase_uncertainties: phase uncertainties to be saved
        :param topography_corrections: topography corrections to be saved
        """

        # Create the xarray Dataset
        da = xr.DataArray(
            data=amplitudes,
            dims=["points"],
            coords={
                "amplitude_uncertainties": ("points", np.asarray(amplitude_uncertainties)),
                "phases": ("points", np.asarray(phases)),
                "phase_uncertainties": ("points", np.asarray(phase_uncertainties)),
                "topography_corrections": ("points", np.asarray(topography_corrections)),
                "x": ("points", np.asarray([point[0] for point in cartesian_coordinates])),
                "y": ("points", np.asarray([point[1] for point in cartesian_coordinates])),
                "z": ("points", np.asarray([point[2] for point in cartesian_coordinates])),
                "latitude": ("points", np.asarray([point[0] for point in geodetic_coordinates])),
                "longitude": ("points", np.asarray([point[1] for point in geodetic_coordinates])),
                "height": ("points", np.asarray([point[2] for point in geodetic_coordinates]))},
            name="amplitudes",
            attrs=dict(
                body_id=self.conOps.body_id,
                start_time=self.start_time,
                end_time=self.end_time
            ),
        )

        # Save xarray to object
        self.data = da

    def calculate_uncertainty(self, uncertainty_c, uncertainty_s, a, phase):

        sigma_a = np.sqrt((uncertainty_c * np.sin(phase) ** 2 - uncertainty_s * np.cos(phase) ** 2) /
                          (np.sin(phase) ** 4 - np.cos(phase) ** 4))
        sigma_phase = np.sqrt((-uncertainty_c * np.cos(phase) ** 2 + uncertainty_s * np.sin(phase) ** 2) /
                              (a ** 2 * (np.sin(phase) ** 4 - np.cos(phase) ** 4)))
        return sigma_a, sigma_phase

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

        line_of_sight_vectors = [[] for _ in range(len(spatial_points))]

        # Iterate through interferograms
        for interferogram in interferograms:

            # Load the data
            dataarray = interferogram.data

            # Extract data
            latitudes = dataarray["latitude"].values
            longitudes = dataarray["longitude"].values
            points = np.column_stack((latitudes, longitudes))

            # Interpolate each desired variable
            results = {var: [] for var in ["time1", "time2", "los_displacements", "sat_pos_time", "psi", "x", "y", "z"]}

            # Loop through target points
            for i in range(len(spatial_points)):

                point = spatial_points[i]
                target_point = np.array([point[0], point[1]])

                # Interpolate each desired variable
                for var in ["time1", "time2", "los_displacements", "sat_pos_time", "psi", "x", "y", "z"]:
                    values = dataarray[var].values
                    interpolated_value = griddata(points=points, values=values, xi=target_point, method="cubic")
                    results[var].extend(interpolated_value)

            # Extract interpolated values
            sat_times = np.asarray(results["sat_pos_time"])
            time1_interps = (np.asarray(results["time1"]) % tidal_cycle)
            time2_interps = (np.asarray(results["time2"]) % tidal_cycle)
            x_interps = np.asarray(results["x"])
            y_interps = np.asarray(results["y"])
            z_interps = np.asarray(results["z"])
            los_displacement_interps = np.asarray(results["los_displacements"])
            psis_interps = np.asarray(results["psi"])

            alpha = 0.3  # Adjust as needed for a tighter fit
            polygon = alphashape.alphashape(points, alpha)
            spatial_points_to_check = [Point(point[0], point[1]) for point in spatial_points]
            is_inside = polygon.contains(spatial_points_to_check)

            # Get the satellite and surface positions
            satellite_positions, satellite_velocities = self.conops.get_states(times=sat_times)
            surface_positions = np.asarray([x_interps, y_interps, z_interps]).T

            # Calculate vectors from satellite to point and the satellite's direction
            vectors_to_sat = (satellite_positions - surface_positions)
            vectors_to_sat = [vector / np.linalg.norm(vector) for vector in vectors_to_sat]

            # Save indices of rows that are added
            valid_point_indices = []
            for i in range(len(spatial_points)):
                # Check if the look angle and relative position is correct
                if is_inside[i]:
                    valid_point_indices.append(i)
                    if not any(np.allclose(vectors_to_sat[i], x, atol=0.1) for x in line_of_sight_vectors[i]):
                        line_of_sight_vectors[i].append(vectors_to_sat[i])

            # Append to d_vectors and G_matrices for each point using list comprehension
            [d_vectors[i].append(los_displacement_interps[i]) for i in valid_point_indices]

            [g_matrices[i].append([np.cos(omega * time2_interps[i]) - np.cos(omega * time1_interps[i]),
                                   np.sin(omega * time2_interps[i]) - np.sin(omega * time1_interps[i]),
                                   psis_interps[i]]) for i in valid_point_indices]

        # Convert lists to numpy arrays and reshape
        d_vectors = {i: np.array(d_vectors[i]) for i in d_vectors}  # Shape: (num_times, d_dim)
        g_matrices = {i: np.array(g_matrices[i]) for i in g_matrices}  # Shape: (num_times, G_dim)

        # Setup data covariance matrices as identity matrices of size amount of rows in g_matrices
        c_ds = {i: np.eye(len(g_matrices[i])) for i in g_matrices}

        # Make a list of amplitudes, phases, topography errors
        amplitudes = []
        amplitude_uncertainties = []
        phases = []
        phase_uncertainties = []
        topography_corrections = []

        for i in tqdm(range(len(spatial_points))):

            if len(line_of_sight_vectors) > 1:
                raise Exception("More than 1 look direction is not supported in 1D inversions")

            # Do the inverse problem
            m = np.linalg.inv(g_matrices[i].T @ np.linalg.inv(c_d) @ g_matrices[i] +
                              np.eye(3)) @ (g_matrices[i].T @ np.linalg.inv(c_d) @ d_vectors[i])

            # Retrieve the parameters
            c, s, z = m

            # Calculate amplitudes and phases
            covariance_matrix = np.linalg.inv(g_matrices[i].T @ np.linalg.inv(c_d) @ g_matrices[i] + np.eye(3))
            a = np.sqrt(c ** 2 + s ** 2)
            phase = np.arctan2(s, c)
            main_diagonal = np.diag(covariance_matrix).tolist()
            amplitudes.append(a)
            phases.append(phase)
            uncertainty_a1, uncertainty_phase1 = self.calculate_uncertainty(uncertainty_c=main_diagonal[0],
                                                                            uncertainty_s=main_diagonal[1],
                                                                            a=a, phase=phase)
            amplitude_uncertainties.append(uncertainty_a1)
            phase_uncertainties.append(uncertainty_phase1)
            topography_corrections.append(z)

        # Return the results
        return amplitudes, amplitude_uncertainties, phases, phase_uncertainties, topography_corrections

# end of file
