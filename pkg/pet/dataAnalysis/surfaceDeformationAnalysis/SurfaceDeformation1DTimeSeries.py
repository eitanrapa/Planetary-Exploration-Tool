#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from scipy.interpolate import griddata
from tqdm import tqdm
import alphashape
from shapely.geometry import Point
import xarray as xr


class SurfaceDeformation1DTimeSeries(
    pet.component, family="pet.dataAnalysis.surfaceDeformationAnalysis.surfaceDeformation1DTimeSeries",
    implements=pet.protocols.dataAnalysis.surfaceDeformationAnalysis):

    """
    Class that is able to create time series information from loaded interferograms
    """

    planet = pet.protocols.planet()
    planet.doc = "target planet"

    campaign = pet.protocols.campaigns.orbiter()
    campaign.doc = "campaign (must be an oribiter)"

    instrument = pet.protocols.instruments.inSAR()
    instrument.doc = "observation instrument"

    data = None

    @classmethod
    def from_file(cls, planet, campaign, instrument, file_name):
        """

        """

        # Open the HDF5 file in read mode
        data = xr.open_dataset(filename_or_obj=file_name)

        obj = cls(planet=planet, campaign=campaign, instrument=instrument)
        obj.data = data  # Restore computed result
        return obj

    @classmethod
    def from_files(cls, planet, campaign, instrument, file_list):
        """

        """

        # Load all the files
        return [cls.from_file(planet=planet, campaign=campaign, instrument=instrument, file_name=file)
                for file in file_list]

    def save(self, file_name):
        """
        Save the track to an HDF5 file
        :return: Nothing returned
        """

        # Open HDF5 file
        self.data.to_netcdf(file_name)

    def create_data_array(self, geodetic_coordinates, amplitudes, amplitude_uncertainties,
                          phases, phase_uncertainties, topography_corrections):
        """
        Create a xarray with the input data
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
                "latitude": ("points", np.asarray([point[0] for point in geodetic_coordinates])),
                "longitude": ("points", np.asarray([point[1] for point in geodetic_coordinates]))},
            name="amplitudes",
            attrs=dict(
                body_id=self.campaign.body_id,
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

    def create_1d_time_series(self, interferograms, spatial_points, alpha=0.3):
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
        for interferogram in tqdm(interferograms, desc="Interpolating interferograms..."):

            # Load the data
            data_array = interferogram.data

            # Extract data
            latitudes = data_array["latitude"].values
            longitudes = data_array["longitude"].values
            points = np.column_stack((latitudes, longitudes))

            # Interpolate each desired variable
            results = {var: [] for var in ["time1", "time2", "los_displacements", "sat_pos_time", "psi", "x", "y", "z"]}

            # Interpolate each desired variable
            for var in ["time1", "time2", "los_displacements", "sat_pos_time", "psi", "x", "y", "z"]:
                values = data_array[var].values
                interpolated_value = griddata(points=points, values=values, xi=spatial_points, method="cubic")
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
            psis_interps = psis_interps * data_array.attrs["baseline"]
            los_displacement_interps_with_noise = (los_displacement_interps +
                                                   psis_interps * self.planet.topography_uncertainty)

            polygon = alphashape.alphashape(points, alpha)
            spatial_points_to_check = [Point(point[0], point[1]) for point in spatial_points]
            is_inside = polygon.contains(spatial_points_to_check)

            # Get the satellite and surface positions
            satellite_positions, satellite_velocities = self.campaign.get_states(times=sat_times)
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
            [d_vectors[i].append(los_displacement_interps_with_noise[i]) for i in valid_point_indices]

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

        for i in tqdm(range(len(spatial_points)), desc="Calculating inverse problems..."):

            if len(line_of_sight_vectors[i]) < 1 or len(line_of_sight_vectors[i]) > 1:
                amplitudes.append(np.nan)
                amplitude_uncertainties.append(np.nan)
                phases.append(np.nan)
                phase_uncertainties.append(np.nan)
                topography_corrections.append(np.nan)

            else:
                # Do the inverse problem
                m = (np.linalg.inv(g_matrices[i].T @ np.linalg.inv(c_ds[i]) @ g_matrices[i])
                     @ (g_matrices[i].T @ np.linalg.inv(c_ds[i]) @ d_vectors[i]))

                # Retrieve the parameters
                c, s, z = m

                # Calculate amplitudes and phases
                covariance_matrix = np.linalg.inv(g_matrices[i].T @ np.linalg.inv(c_ds[i]) @ g_matrices[i])
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

    def visualize_time_series_amplitudes(self, projection, fig=None, globe=None, ax=None, return_fig=False):
        """

        """

        if fig is None:
            # Get the projection
            fig, ax, globe = projection.proj(planet=self.planet)

        # Load the values to plot
        amplitudes = self.data.values
        longitudes = self.data["longitude"].values
        latitudes = self.data["latitude"].values

        # Make the colormap cyclical
        cm = plt.cm.get_cmap('hsv')

        # Plot the amplitudes
        im = ax.scatter(longitudes, latitudes, c=amplitudes, cmap=cm, transform=ccrs.PlateCarree(globe=globe))

        # return fig, ax, globe if necessary
        if return_fig:
            return fig, ax, globe

        # Add a colorbar
        plt.colorbar(im, ax=ax, orientation='vertical', shrink=0.25)

        # Add labels and legend
        plt.title('Displacement amplitudes', pad=20)

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/' + 'time_series_amplitudes_' +
                          self.deformation_map.pyre_name + '_' + str(self.campaign.body_id) + '.png', format='png',
                    dpi=500)

    def visualize_time_series_phases(self, projection, fig=None, globe=None, ax=None, return_fig=False):
        """

        """

        if fig is None:
            # Get the projection
            fig, ax, globe = projection.proj(planet=self.planet)

        # Load the values to plot
        phases = self.data["phases"].values
        longitudes = self.data["longitude"].values
        latitudes = self.data["latitude"].values

        # Make the colormap cyclical
        cm = plt.cm.get_cmap('hsv')

        # Plot the amplitudes
        im = ax.scatter(longitudes, latitudes, c=phases, cmap=cm, transform=ccrs.PlateCarree(globe=globe))

        # return fig, ax, globe if necessary
        if return_fig:
            return fig, ax, globe

        # Add a colorbar
        plt.colorbar(im, ax=ax, orientation='vertical', shrink=0.25)

        # Add labels and legend
        plt.title('Displacement phases', pad=20)

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/' + 'time_series_phases_' +
                          self.deformation_map.pyre_name + '_' + str(self.campaign.body_id) + '.png', format='png',
                    dpi=500)

# end of file
