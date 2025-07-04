#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import pet
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from scipy.interpolate import griddata
from tqdm import tqdm
import alphashape
from shapely.geometry import Point
import xarray as xr
from multiprocessing import Pool


class SurfaceDeformation1DTimeSeries(pet.component,
                                     family="pet.dataAnalysis.surfaceDeformationAnalysis."
                                            "surfaceDeformation1DTimeSeries",
                                     implements=pet.protocols.dataAnalysis.surfaceDeformationAnalysis):
    """
    Class that is able to create time series information from loaded interferograms
    """

    planet = pet.protocols.planet()
    planet.doc = "target planet"

    campaign = pet.protocols.campaigns.orbiter()
    campaign.doc = "campaign (must be an orbiter)"

    instrument = pet.protocols.instruments.inSAR()
    instrument.doc = "observation instrument"

    data = None

    @classmethod
    def from_file(cls, planet, campaign, instrument, file_name):
        """
        Load a track from an HDF5 file
        :param planet: Planet object
        :param campaign: Campaign object
        :param instrument: Instrument object
        :param file_name: Name of the file to load
        :return: A track object
        """

        # Open the HDF5 file in read mode
        data = xr.open_dataset(filename_or_obj=file_name)

        # Create the object
        obj = cls(name="time_series" + str(np.random.rand()), planet=planet, campaign=campaign, instrument=instrument)

        obj.data = data  # Restore computed result

        return obj

    @classmethod
    def from_files(cls, planet, campaign, instrument, file_list):
        """
        Load a list of tracks from HDF5 files
        :param planet: Planet object
        :param campaign: Campaign object
        :param instrument: Instrument object
        :param file_list: List of file names to load
        :return: A list of track objects
        """

        # Load all the files
        return [cls.from_file(planet=planet, campaign=campaign, instrument=instrument,
                              file_name=file) for file in file_list]

    def save(self, file_name):
        """
        Save the track to an HDF5 file
        :param file_name: Name of the file to save
        :return: Nothing returned
        """

        # Open HDF5 file
        self.data.to_netcdf(file_name, engine="netcdf4")

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
        :return: Nothing returned
        """

        # Create the xarray Dataset
        da = xr.DataArray(
            data=topography_corrections,
            dims=["points"],
            coords={
                "amplitudes": ("points", np.asarray(amplitudes)),
                "amplitude_uncertainties":
                    ("points", np.asarray(amplitude_uncertainties)),
                "phases": ("points", np.asarray(phases)),
                "phase_uncertainties":
                    ("points", np.asarray(phase_uncertainties)),
                "latitude": ("points", np.asarray([point[0] for point in geodetic_coordinates])),
                "longitude": ("points", np.asarray([point[1] for point in geodetic_coordinates]))},
            name="topography_corrections",
            attrs=dict(
                body_id=self.campaign.body_id,
            ),
        )

        # Save xarray to object
        self.data = da

    @staticmethod
    def calculate_uncertainty(uncertainty_c, uncertainty_s, a, phase):
        """
        Calculate the uncertainty of the amplitude and phase
        :param uncertainty_c: uncertainty of the c parameter
        :param uncertainty_s: uncertainty of the s parameter
        :param a: amplitude
        :param phase: phase
        :return: uncertainty of the amplitude, uncertainty of the phase
        """

        # Calculate the uncertainties
        sigma_a = np.sqrt((uncertainty_c * np.sin(phase) ** 2 - uncertainty_s * np.cos(phase) ** 2) /
                          (np.sin(phase) ** 4 - np.cos(phase) ** 4))
        sigma_phase = np.sqrt((-uncertainty_c * np.cos(phase) ** 2 + uncertainty_s * np.sin(phase) ** 2) /
                              (a ** 2 * (np.sin(phase) ** 4 - np.cos(phase) ** 4)))

        return sigma_a, sigma_phase

    @staticmethod
    def parallel_interpolate_interferogram(args):
        """
        Use parallel processing to interpolate an interferogram
        :param args: interferogram, spatial_points, alpha, shared_data, progress_queue,
         (line_of_sight_vectors, d_vectors, g_matrices)
        :return: Nothing returned
        """

        forward_data_array, inverse_data_array, spatial_points, alpha, tidal_cycle, wavelength = args

        # Extract data
        latitudes = inverse_data_array["latitude"].values
        longitudes = inverse_data_array["longitude"].values
        points = np.column_stack((latitudes, longitudes))

        # Create alpha shape (polygon)
        polygon = alphashape.alphashape(points, alpha)

        # Check which spatial points are inside the polygon
        spatial_points_to_check = np.array([Point(point[0], point[1]) for point in spatial_points])
        is_inside = np.array([polygon.contains(pt) for pt in spatial_points_to_check])

        # Extract valid spatial points
        valid_points = spatial_points[is_inside]

        # Initialize results dictionary
        results = {var: np.full(len(spatial_points), np.nan) for var in
                   ["time", "los_displacement_ref", "kz_k", "x", "y", "z", "sigma_phase", "nlooks",
                                                                                          "correlation"]}

        # Interpolate each desired variable
        for var in ["time", "los_displacement_ref", "kz_k", "x", "y", "z"]:
            values = inverse_data_array[var].values

            # Perform interpolation only for valid points
            interpolated_values = griddata(points=points, values=values, xi=valid_points, method="cubic")

            # Assign interpolated values back to the corresponding positions
            results[var][is_inside] = interpolated_values

        # Interpolate each desired variable
        for var in ["sigma_phase", "nlooks", "correlation"]:
            values = forward_data_array[var].values

            # Perform interpolation only for valid points
            interpolated_values = griddata(points=points, values=values, xi=valid_points, method="cubic")

            # Assign interpolated values back to the corresponding positions
            results[var][is_inside] = interpolated_values

        # Extract interpolated values
        times = np.asarray(results["time"])
        x_interps = np.asarray(results["x"])
        y_interps = np.asarray(results["y"])
        z_interps = np.asarray(results["z"])
        sigma_phase = np.asarray(results["sigma_phase"])
        nlooks = np.asarray(results["nlooks"])
        correlation = np.asarray(results["correlation"])
        phase_noise = np.where(
            np.isnan(sigma_phase),
            np.nan,
            np.random.normal(loc=0, scale=np.abs(sigma_phase))
        )
        los_displacement_interps_with_noise = (np.asarray(results["los_displacement_ref"])
                                               + phase_noise / (4 * np.pi / wavelength))
        psis = np.asarray(results["kz_k"]) / (4 * np.pi) * wavelength

        times1 = (times + inverse_data_array.attrs["time_offset_first_acquisition"]) % tidal_cycle
        times2 = (times + inverse_data_array.attrs["time_offset_second_acquisition"]) % tidal_cycle

        # Prepare data for return
        valid_point_indices = np.where(is_inside)[0]

        return (times, times1, times2, valid_point_indices, los_displacement_interps_with_noise, psis, x_interps,
                y_interps, z_interps, correlation, nlooks)

    def create_1d_time_series(self, interferograms, spatial_points, alpha=1, processors=1):
        """
        Return the amplitudes, phases, and topographic errors of a set of points in 1 look direction
        :param interferograms: interferograms to analyze, in 1 look direction
        :param spatial_points: points to do the inverse problem on
        :param alpha: alpha parameter for the alpha shape
        :param processors: number of processors to use
        :return: amplitudes, phases, topographic errors
        """

        # Set up a matrix of observations and a problem G matrix
        d_vectors = {i: [] for i in range(len(spatial_points))}  # List to collect rows for each point
        g_matrices = {i: [] for i in range(len(spatial_points))}

        # Set up data covariance matrices
        cd_s = {i: [] for i in range(len(spatial_points))}

        # Define the tidal cycle, and therefore the angular frequency
        tidal_cycle = self.planet.tidal_cycle
        omega = 2 * np.pi / tidal_cycle

        line_of_sight_vectors = [[] for _ in range(len(spatial_points))]

        if processors == 1:

            # Iterate through interferograms
            for interferogram in tqdm(interferograms, desc="Interpolating interferograms..."):

                # Load the data
                forward_data_array = interferogram.forward_data
                inverse_data_array = interferogram.inverse_data

                # Extract data
                latitudes = inverse_data_array["latitude"].values
                longitudes = inverse_data_array["longitude"].values
                points = np.column_stack((latitudes, longitudes))

                # Create alpha shape (polygon)
                polygon = alphashape.alphashape(points, alpha)

                # Check which spatial points are inside the polygon
                spatial_points_to_check = np.array([Point(point[0], point[1]) for point in spatial_points])
                is_inside = np.array([polygon.contains(pt) for pt in spatial_points_to_check])

                # Extract valid spatial points
                valid_points = spatial_points[is_inside]

                # Initialize results dictionary
                results = {var: np.full(len(spatial_points), np.nan) for var in
                           ["time", "los_displacement_ref", "kz_k", "x", "y", "z", "sigma_phase", "nlooks",
                                                                                                  "correlation"]}

                # Interpolate each desired variable
                for var in ["time", "los_displacement_ref", "kz_k", "x", "y", "z"]:

                    values = inverse_data_array[var].values

                    # Perform interpolation only for valid points
                    interpolated_values = griddata(points=points, values=values, xi=valid_points, method="cubic")

                    # Assign interpolated values back to the corresponding positions
                    results[var][is_inside] = interpolated_values

                # Interpolate each desired variable
                for var in ["sigma_phase", "nlooks","correlation"]:

                    values = forward_data_array[var].values

                    # Perform interpolation only for valid points
                    interpolated_values = griddata(points=points, values=values, xi=valid_points, method="cubic")

                    # Assign interpolated values back to the corresponding positions
                    results[var][is_inside] = interpolated_values

                # Extract interpolated values
                times = np.asarray(results["time"])
                x_interps = np.asarray(results["x"])
                y_interps = np.asarray(results["y"])
                z_interps = np.asarray(results["z"])
                sigma_phase = np.asarray(results["sigma_phase"])
                nlooks = np.asarray(results["nlooks"])
                correlation = np.asarray(results["correlation"])
                phase_noise = np.where(
                    np.isnan(sigma_phase),
                    np.nan,
                    np.random.normal(loc=0, scale=np.abs(sigma_phase))
                )
                los_displacement_interps_with_noise = (np.asarray(results["los_displacement_ref"])
                                            + phase_noise / (4 * np.pi / self.instrument.wavelength))
                psis = np.asarray(results["kz_k"]) / (4 * np.pi) * self.instrument.wavelength
                times1 = (times + inverse_data_array.attrs["time_offset_first_acquisition"]) % tidal_cycle
                times2 = (times + inverse_data_array.attrs["time_offset_second_acquisition"]) % tidal_cycle

                # Get the satellite and surface positions
                satellite_positions, _ = self.campaign.get_states(times=times)
                surface_positions = np.asarray([x_interps, y_interps, z_interps]).T

                # Calculate vectors from satellite to point and the satellite's direction
                vectors_to_sat = (satellite_positions - surface_positions)
                vectors_to_sat = [vector / np.linalg.norm(vector) for vector in vectors_to_sat]

                valid_point_indices = np.where(is_inside)[0]

                for i in valid_point_indices:
                    if not any(np.allclose(vectors_to_sat[i], x, atol=0.1) for x in line_of_sight_vectors[i]):
                        line_of_sight_vectors[i].append(vectors_to_sat[i])

                    cd_s[i].append(np.sqrt(self.instrument.wavelength / (4 * np.pi * nlooks[i] *
                                                                          (1 - correlation[i]**2) /
                                                                          correlation[i]**2)))
                    d_vectors[i].append(los_displacement_interps_with_noise[i])
                    g_matrices[i].append([np.cos(omega * times2[i]) - np.cos(omega * times1[i]),
                                          np.sin(omega * times2[i]) - np.sin(omega * times1[i]),
                                          psis[i]])

        else:

            with Pool(processors) as pool:

                results = list(tqdm(pool.imap(self.parallel_interpolate_interferogram,
                                              [(interferogram.forward_data, interferogram.inverse_data,
                                                spatial_points, alpha, tidal_cycle, self.instrument.wavelength)
                                               for interferogram in interferograms]),
                                    total=len(interferograms), desc="Interpolating interferograms..."))

            for (times, times1, times2, valid_point_indices, los_displacement_interps_with_noise, psis,
                 x_interps, y_interps, z_interps, correlation, nlooks) in results:

                # Get satellite and surface positions
                satellite_positions, _ = self.campaign.get_states(times=times)
                surface_positions = np.asarray([x_interps, y_interps, z_interps]).T

                # Compute vectors to satellite
                vectors_to_sat = [(satellite_positions[i] - surface_positions[i]) / np.linalg.norm(
                    satellite_positions[i] - surface_positions[i])
                                  for i in range(len(surface_positions))]

                for i in valid_point_indices:
                    if not any(np.allclose(vectors_to_sat[i], x, atol=0.1) for x in line_of_sight_vectors[i]):
                        line_of_sight_vectors[i].append(vectors_to_sat[i])

                    cd_s[i].append(np.sqrt(self.instrument.wavelength / (4 * np.pi * nlooks[i] *
                                                                         (1 - correlation[i] ** 2) /
                                                                         correlation[i] ** 2)))
                    d_vectors[i].append(los_displacement_interps_with_noise[i])
                    g_matrices[i].append([np.cos(omega * times2[i]) - np.cos(omega * times1[i]),
                                          np.sin(omega * times2[i]) - np.sin(omega * times1[i]),
                                          psis[i]])

        # Convert lists to numpy arrays and reshape
        d_vectors = {i: np.array(d_vectors[i]) for i in d_vectors}  # Shape: (num_times, d_dim)
        g_matrices = {i: np.array(g_matrices[i]) for i in g_matrices}  # Shape: (num_times, G_dim)
        cd_s = {i: np.diag(cd_s[i]) for i in cd_s}  # Shape: (num_times, num_times)

        # Make a list of amplitudes, phases, topography errors
        amplitudes = []
        amplitude_uncertainties = []
        phases = []
        phase_uncertainties = []
        topography_corrections = []

        for i in tqdm(range(len(spatial_points)), desc="Calculating inverse problems..."):

            # Check if there is exactly one line of sight vector
            if len(line_of_sight_vectors[i]) < 1 or len(line_of_sight_vectors[i]) > 1:
                amplitudes.append(np.nan)
                amplitude_uncertainties.append(np.nan)
                phases.append(np.nan)
                phase_uncertainties.append(np.nan)
                topography_corrections.append(np.nan)

            else:
                # Do the inverse problem
                m = (np.linalg.inv(g_matrices[i].T @ np.linalg.inv(cd_s[i]) @ g_matrices[i])
                     @ (g_matrices[i].T @ np.linalg.inv(cd_s[i]) @ d_vectors[i]))

                # Retrieve the parameters
                c, s, z = m

                # Calculate amplitudes and phases
                covariance_matrix = np.linalg.inv(g_matrices[i].T @ np.linalg.inv(cd_s[i]) @ g_matrices[i])
                a = np.sqrt(c ** 2 + s ** 2)

                # If a is negative, flip the sign and add pi to the phase
                if a < 0:
                    a = -a
                    phase = np.arctan2(c, s) + np.pi
                else:
                    phase = np.arctan2(c, s)

                # Append the results
                amplitudes.append(a)
                phases.append(phase)
                topography_corrections.append(z)

                # Calculate uncertainties
                main_diagonal = np.diag(covariance_matrix).tolist()
                uncertainty_a1, uncertainty_phase1 = self.calculate_uncertainty(uncertainty_c=main_diagonal[0],
                                                                                uncertainty_s=main_diagonal[1],
                                                                                a=a, phase=phase)

                # Append the uncertainties
                amplitude_uncertainties.append(uncertainty_a1)
                phase_uncertainties.append(uncertainty_phase1)

        # Create the data array
        self.create_data_array(geodetic_coordinates=spatial_points, amplitudes=amplitudes,
                               amplitude_uncertainties=amplitude_uncertainties, phases=phases,
                               phase_uncertainties=phase_uncertainties, topography_corrections=topography_corrections)

    def visualize_time_series_amplitudes(self, projection, direction, fig=None, globe=None, ax=None,
                                         return_fig=False):
        """
        Visualize the amplitudes of the time series
        :param projection: Cartopy projection
        :param direction: east, north, or up
        :param fig: matplotlib figure
        :param globe: cartopy globe
        :param ax: matplotlib ax
        :param return_fig: Whether to return fig, globe, ax
        :return: fig, ax, globe if return_fig is True
        """

        if fig is None:
            # Get the projection
            fig, ax, globe = projection.proj(planet=self.planet)

        # Load the values to plot
        if direction == "east":
            amplitudes = self.data["amplitudes_east"].values
        elif direction == "north":
            amplitudes = self.data["amplitudes_north"].values
        elif direction == "up":
            amplitudes = self.data["amplitudes_up"].values
        else:
            raise Exception("direction must be east, north, or up")

        longitudes = self.data["longitude"].values
        latitudes = self.data["latitude"].values

        # Plot the amplitudes
        im = ax.scatter(longitudes, latitudes, c=amplitudes, transform=ccrs.PlateCarree(globe=globe))

        # return fig, ax, globe if necessary
        if return_fig:
            return fig, ax, globe

        # Add a colorbar
        plt.colorbar(im, ax=ax, orientation='vertical', shrink=0.25)

        # Add labels and legend
        plt.title('Displacement amplitudes', pad=20)

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/' + 'time_series_amplitudes_' + str(direction) +
                          str(self.campaign.body_id) + '.png', format='png',
                    dpi=500)

    def visualize_time_series_phases(self, projection, direction, fig=None, globe=None, ax=None, return_fig=False):
        """
        Visualize the phases of the time series
        :param projection: Cartopy projection
        :param direction: east, north, or up
        :param fig: matplotlib figure
        :param globe: cartopy globe
        :param ax: matplotlib ax
        :param return_fig: Whether to return fig, globe, ax
        :return: fig, ax, globe if return_fig is True
        """

        if fig is None:
            # Get the projection
            fig, ax, globe = projection.proj(planet=self.planet)

        # Load the values to plot
        if direction == "east":
            phases = self.data["phases_east"].values
        elif direction == "north":
            phases = self.data["phases_north"].values
        elif direction == "up":
            phases = self.data["phases_up"].values
        else:
            raise Exception("direction must be east, north, or up")
        longitudes = self.data["longitude"].values
        latitudes = self.data["latitude"].values

        # Plot the amplitudes
        im = ax.scatter(longitudes, latitudes, c=phases, transform=ccrs.PlateCarree(globe=globe))

        # return fig, ax, globe if necessary
        if return_fig:
            return fig, ax, globe

        # Add a colorbar
        plt.colorbar(im, ax=ax, orientation='vertical', shrink=0.25)

        # Add labels and legend
        plt.title('Displacement phases', pad=20)

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/' + 'time_series_phases_' + str(direction) +
                          str(self.campaign.body_id) + '.png', format='png',
                    dpi=500)

    # end of file
