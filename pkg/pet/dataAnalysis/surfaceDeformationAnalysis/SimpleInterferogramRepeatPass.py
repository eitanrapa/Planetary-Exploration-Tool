#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import pet
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import sys


class SimpleInterferogramRepeatPass(pet.component,
                                    family="pet.dataAnalysis.surfaceDeformationAnalysis.simpleInterferogram",
                                    implements=pet.protocols.dataAnalysis.surfaceDeformationAnalysis):
    """
    Class that creates a single interferogram between two acquisitions given an instrument and orbit
    """


    planet = pet.protocols.planet()
    planet.doc = "target planet"

    campaign = pet.protocols.campaigns.orbiter()
    campaign.doc = "campaign (must be an orbiter)"

    instrument = pet.protocols.instruments.inSAR()
    instrument.doc = "observation instrument"

    deformation_map = pet.protocols.natureSimulations.geophysicalModel()
    deformation_map.doc = "tidal deformation model"

    track = pet.protocols.dataAcquisition()
    track.doc = "ground track"

    time_offset_first_acquisition = pet.properties.float()
    time_offset_first_acquisition.default = 0
    time_offset_first_acquisition.doc = "time offset between satellite position time and the first acquisition time"

    time_offset_second_acquisition = pet.properties.float()
    time_offset_second_acquisition.default = 0
    time_offset_second_acquisition.doc = "time offset between satellite position time and the second acquisition time"

    perpendicular_baseline = pet.properties.float()
    perpendicular_baseline.doc = "perpendicular baseline between the two acquisitions [m]"

    data = None

    @classmethod
    def from_file(cls, planet, instrument, campaign, deformation_map, file_name):
        """
        Load the interferogram from an HDF5 file
        :param planet: Planet object
        :param instrument: Instrument object
        :param campaign: Campaign object
        :param deformation_map: Deformation map object
        :param file_name: Name of the file to load
        :return: Interferogram object
        """

        # Open the HDF5 file in read mode
        datatree = xr.open_datatree(filename_or_obj=file_name)
        track_da = datatree["track"].to_dataset()["track"]
        data = datatree["interferogram"]["phases"]

        # Create the track object
        track = pet.dataAcquisition.track.from_data_array(planet=planet, campaign=campaign, instrument=instrument,
                                                          data=track_da)

        # Create the interferogram object
        obj = cls(name="igram" + str(np.random.rand()),
                  planet=planet, instrument=instrument, campaign=campaign,
                  deformation_map=deformation_map, track=track,
                  time_offset_first_acquisition=data.attrs["time_offset_first_acquisition"],
                  time_offset_second_acquisition=data.attrs["time_offset_second_acquisition"],
                  perpendicular_baseline=data.attrs["perpendicular_baseline"])

        obj.data = data  # Restore computed result

        return obj

    @classmethod
    def from_files(cls, planet, instrument, campaign, deformation_map, file_list):
        """
        Load a list of interferograms from a list of HDF5 files
        :param planet: Planet object
        :param instrument: Instrument object
        :param campaign: Campaign object
        :param deformation_map: Deformation map object
        :param file_list: List of files to load
        :return: List of interferogram objects
        """

        # Load all the files
        return [cls.from_file(planet=planet, instrument=instrument, campaign=campaign, deformation_map=deformation_map,
                              file_name=file) for file in file_list]

    def save(self, file_name):
        """
        Save the track to an HDF5 file
        :param file_name: Name of the file to save
        :return: Nothing returned
        """

        # Open HDF5 file
        self.data.to_netcdf(file_name, engine="netcdf4")

    def create_data_array(self, phases, psis, corrs, nlooks):
        """
        Create a xarray with the input data
        :param phases: forward phases of interferogram
        :param psis: Flattened values
        :return: Nothing returned
        """

        # Create the xarray datarray
        phases_da = xr.DataArray(
            data=phases,
            dims=["points"],
            coords={
                "time1": ("points", self.track.data["time"].values + self.time_offset_first_acquisition),
                "time2": ("points", self.track.data["time"].values + self.time_offset_second_acquisition),
                "psi_sub_baseline": ("points", psis),
                "correlation": ("points", corrs),
                "nlooks": ("points", nlooks),
                "x": ("points", self.track.data["x"].values),
                "y": ("points", self.track.data["y"].values),
                "z": ("points", self.track.data["z"].values),
                "latitude": ("points", self.track.data["latitude"].values),
                "longitude": ("points", self.track.data["longitude"].values),
                "height": ("points", self.track.data["height"].values)},
            name="phases",
            attrs=dict(
                deformation_map=self.deformation_map.pyre_name,
                body_id=self.campaign.body_id,
                perpendicular_baseline=self.perpendicular_baseline,
                time_offset_first_acquisition=self.time_offset_first_acquisition,
                time_offset_second_acquisition=self.time_offset_second_acquisition,
                start_time1=self.track.start_time + self.time_offset_first_acquisition,
                end_time1=self.track.end_time + self.time_offset_first_acquisition,
                start_time2=self.track.start_time + self.time_offset_second_acquisition,
                end_time2=self.track.end_time + self.time_offset_second_acquisition,
            ),
        )

        # Store them in a DataTree
        dt = xr.DataTree.from_dict({
            "interferogram": xr.Dataset({"phases": phases_da}),
            "track": xr.Dataset({"track": self.track.data})
        })

        self.data = dt

    def get_flattened_angles(self, time, satellite_position, flat_positions):
        """
        Calculate the flattened angle from a given set of positions at a time in the satellite orbit
        :param flat_positions: Flattened positions
        :param time: Time at which the satellite is in satellite_position
        :param satellite_position: Position of the satellite at a specific time
        :return: Flattened angles of satellite to ground positions
        """

        # Calculate the satellite intersect and height with the shape
        intersects, satellite_heights = self.planet.get_sub_obs_points(times=time, campaign=self.campaign)

        # Calculate the distance between center and satellite intersects
        satellite_radii = np.asarray([np.linalg.norm(intersect - [0, 0, 0]) for intersect in intersects])

        # Get the distance between the satellite and the surface points
        vectors = flat_positions - satellite_position
        distances = np.asarray([np.linalg.norm(vector) for vector in vectors])

        # Calculate cosines of the angles
        z_plus_re = satellite_heights + satellite_radii
        altitudes = np.asarray([np.linalg.norm(flat_position - [0, 0, 0])
                                for flat_position in flat_positions])
        cosine_look_angle = (distances ** 2 + z_plus_re ** 2 - altitudes ** 2) / (2 * distances * z_plus_re)

        # Use arc cosine to get the angles in radians
        look_angles_radians = np.arccos(cosine_look_angle)

        # Return the look angles 2-D array
        return look_angles_radians

    def calculate_igram(self):
        """
        Calculate the flattened phases between two swaths given a baseline
        :return: Nothing returned 
        """

        print("Getting surface displacements...", file=sys.stderr)

        # For speed, calculate all points at once
        u_displacements, v_displacements, w_displacements = (
            self.deformation_map.get_displacements(track=self.track,
                                                   time_difference=self.time_offset_first_acquisition))
        displacements_1 = np.array([u_displacements, v_displacements, w_displacements]).T

        # For speed, calculate all points at once
        u_displacements, v_displacements, w_displacements = (
            self.deformation_map.get_displacements(track=self.track,
                                                   time_difference=self.time_offset_second_acquisition))
        displacements_2 = np.array([u_displacements, v_displacements, w_displacements]).T

        # Access x, y, z values
        x = self.track.data["x"].values
        y = self.track.data["y"].values
        z = self.track.data["z"].values

        # Get the positions
        positions = np.asarray([x, y, z]).T

        print("Getting satellite positions...", file=sys.stderr)

        # Get satellite positions and velocities
        satellite_positions, satellite_velocities = \
            self.campaign.get_states(times=self.track.data["time"].values)

        # Get the line-of-sight vectors
        los_vectors = satellite_positions - positions

        # Calculate the second satellite positions, assuming a b_perp
        v_unit = [v_sat / np.linalg.norm(v_sat) for v_sat in satellite_velocities]
        los_unit = [los / np.linalg.norm(los) for los in los_vectors]

        cross_tracks = [np.cross(v_unit, los_unit) for v_unit, los_unit in zip(v_unit, los_unit)]
        cross_track_unit = np.asarray([cross_track / np.linalg.norm(cross_track) for cross_track in cross_tracks])

        # Calculate the second satellite positions
        satellite_positions_plus_pb = satellite_positions + self.perpendicular_baseline * cross_track_unit

        # Get the line-of-sight vectors
        los_vectors_plus_pb = satellite_positions_plus_pb - positions

        distances = np.asarray([np.linalg.norm(vector) for vector in los_vectors])
        distances_plus_pb = np.asarray([np.linalg.norm(vector) for vector in los_vectors_plus_pb])

        print("Getting flattened points on the ellipsoid...", file=sys.stderr)

        flattened_points = self.get_flattened_positions(positions=positions, satellite_position=satellite_positions)

        print("Getting flattened angles...", file=sys.stderr)

        # Get flattened angles for the satellite positions and ground points
        angles_0 = self.get_flattened_angles(time=self.track.data["time"].values,
                                             satellite_position=satellite_positions, flat_positions=flattened_points)

        # Get the distances from the flattened points to the satellites
        vectors = flattened_points - satellite_positions
        distances_0 = np.asarray([np.linalg.norm(vector) for vector in vectors])

        # Calculate the vectors from the ground points to the satellite positions
        ground_satellite_vectors = satellite_positions - positions

        print("Calculating line-of-sight displacements...", file=sys.stderr)

        # Calculate the line of sight displacements given the satellite positions for the first swath
        los_displacements1 = np.asarray([(np.dot(displacement, ground_satellite_vector) /
                                          np.linalg.norm(ground_satellite_vector))
                                         for displacement, ground_satellite_vector
                                         in zip(displacements_1, ground_satellite_vectors)])

        # Calculate the line of sight displacements given the satellite positions for the second swath
        los_displacements2 = np.asarray([(np.dot(displacement, ground_satellite_vector) /
                                          np.linalg.norm(ground_satellite_vector))
                                         for displacement, ground_satellite_vector
                                         in zip(displacements_2, ground_satellite_vectors)])

        # Calculate the noise
        sigma_noises, corrs, nlooks = self.instrument.get_instrument_noise(
            planet=self.planet, baseline=self.perpendicular_baseline,
            satellite_velocities=satellite_velocities, look_angle=look_angles,
            incidence_angles=incidence_angles, distances=distances)

        forward_phases = (4 * np.pi * (1 / self.instrument.wavelength) * ((los_displacements2 - los_displacements1) +
                                                                          (distances_plus_pb - distances)) +
                          np.random.normal(loc=0, scale=sigma_noises, size=len(los_displacements1)))

        # Calculate psis sub perpendicular baselines
        psis_sub_baseline = 1 / (distances_0 * np.sin(angles_0))

        # Create array and save the data
        self.create_data_array(phases=forward_phases, psis=psis_sub_baseline, corrs=corrs, nlooks=nlooks)

    def get_inverse_phases(self, position1_error, position2_error, topography_error):
        """

        """

        # Access x, y, z values
        x = self.track.data["x"].values
        y = self.track.data["y"].values
        z = self.track.data["z"].values

        # Get the positions
        positions = np.asarray([x, y, z]).T

        print("Getting satellite positions...", file=sys.stderr)

        # Get satellite positions and velocities
        satellite_positions, satellite_velocities = \
            self.campaign.get_states(times=self.track.data["time"].values)

        # Get the line-of-sight vectors
        los_vectors = satellite_positions - positions

        # Calculate the second satellite positions, assuming a b_perp
        v_unit = [v_sat / np.linalg.norm(v_sat) for v_sat in satellite_velocities]
        los_unit = [los / np.linalg.norm(los) for los in los_vectors]

        cross_tracks = [np.cross(v_unit, los_unit) for v_unit, los_unit in zip(v_unit, los_unit)]
        cross_track_unit = np.asarray([cross_track / np.linalg.norm(cross_track) for cross_track in cross_tracks])

        # Calculate the second satellite positions
        satellite_positions_plus_pb = satellite_positions + self.perpendicular_baseline * cross_track_unit

        # Get the line-of-sight vectors
        los_vectors_plus_pb = satellite_positions_plus_pb - positions

        delta_distances = (np.asarray([np.linalg.norm(vector) for vector in los_vectors]) + position1_error +
                           topography_error)
        delta_distances_plus_pb = (np.asarray([np.linalg.norm(vector) for vector in los_vectors_plus_pb]) +
                                   position2_error + topography_error)

        phases_reference = 4 * np.pi * (1 / self.instrument.wavelength) * (delta_distances_plus_pb - delta_distances)

        return self.data.values - phases_reference

    def visualize_interferogram(self, projection, fig=None, globe=None, ax=None, return_fig=False):
        """
        Visualize the interferogram
        :param projection: Projection object to use for plotting
        :param fig: matplotlib figure
        :param globe: cartopy globe
        :param ax: matplotlib ax
        :param return_fig: Whether to return fig, globe, ax
        :return: Nothing returned
        """

        if fig is None:
            # Get the projection
            fig, ax, globe = projection.proj(planet=self.planet)

        # Load the values to plot
        longitudes = self.data["longitude"].values
        latitudes = self.data["latitude"].values
        phases = self.data.values

        # Wrap interferogram
        phases = np.fmod(phases, 2 * np.pi)
        phases = [angle + 2 * np.pi if angle < 0 else angle for angle in phases]

        # Make the colormap cyclical
        cm = plt.cm.get_cmap('hsv')

        # Iterate through the interferogram beams
        im = ax.scatter(longitudes, latitudes, vmin=0, vmax=2 * np.pi, cmap=cm,
                        transform=ccrs.PlateCarree(globe=globe), c=phases, marker='o', s=0.1)

        # return fig, ax, globe if necessary
        if return_fig:
            return fig, ax, globe

        # Add colorbar
        plt.colorbar(im, label="Phase")

        # Add labels and legend
        plt.title('Interferogram', pad=20)

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/' + 'interferogram_' + self.deformation_map.pyre_name + '_' +
                    str(self.campaign.body_id) +
                    '_' + str(self.track.start_time + self.time_offset_first_acquisition) + '_' + str(
                        self.track.end_time + self.time_offset_first_acquisition) +
                    '_' + str(self.track.start_time + self.time_offset_second_acquisition) +
                    '_' + str(self.track.end_time + self.time_offset_second_acquisition) + '.png', format='png',
                    dpi=500)

    def visualize_displacements(self, projection, fig=None, globe=None, ax=None, return_fig=False):
        """
        Visualize the interferogram
        :param projection: Projection object to use for plotting
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
        longitudes = self.data["longitude"].values
        latitudes = self.data["latitude"].values
        phases = self.get_inverse_phases(position1_error=0, position2_error=0, topography_error=0)
        los_displacements = self.instrument.wavelength * phases / (4 * np.pi)

        im = ax.scatter(longitudes, latitudes,
                        transform=ccrs.PlateCarree(globe=globe),
                        c=los_displacements, marker='o', s=0.1)

        # return fig, ax, globe if necessary
        if return_fig:
            return fig, ax, globe

        # Add colorbar
        plt.colorbar(im, label="Displacement")

        # Add labels and legend
        plt.title('LOS displacements', pad=20)

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/' + 'displacements_' + self.deformation_map.pyre_name + '_' +
                    str(self.campaign.body_id) +
                    '_' + str(self.track.start_time + self.time_offset_first_acquisition) + '_' + str(
                        self.track.end_time + self.time_offset_first_acquisition) +
                    '_' + str(self.track.start_time + self.time_offset_second_acquisition) +
                    '_' + str(self.track.end_time + self.time_offset_second_acquisition) + '.png', format='png',
                    dpi=500)

# end of file
