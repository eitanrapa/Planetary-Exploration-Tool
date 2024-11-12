#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import xarray as xr
import git
import cspyce as spice
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


class Interferogram(pet.component):
    """
    Class that creates a single interferogram between two points of a displacement map given an instrument orbit
    """

    def __init__(self, planet, instrument, deformation_map, conops, track1, track2, baseline, **kwargs):
        super().__init__(**kwargs)
        self.planet = planet
        self.instrument = instrument
        self.conOps = conops
        self.deformation_map = deformation_map
        self.track1 = track1
        self.track2 = track2
        self.baseline = baseline
        self.data = None

    def save(self):
        """
        Save the interferogram to an HDF5 file
        :return: Nothing returned
        """

        body_id = self.track1.data.attrs["body_id"]
        start_time1 = self.track1.data.attrs["start_time"]
        end_time1 = self.track1.data.attrs["end_time"]
        start_time2 = self.track2.data.attrs["start_time"]
        end_time2 = self.track2.data.attrs["end_time"]

        # Open HDF5 file
        repo = git.Repo('.', search_parent_directories=True)
        self.data.to_netcdf(repo.working_tree_dir + '/files/igram_' + self.deformation_map.pyre_name + '_' +
                            str(body_id) + '_' + str(start_time1) + '_' + str(end_time1) +
                            str(start_time2) + '_' + str(end_time2) + ".nc")

    def load(self):
        """
        Load a hdf5 file containing a previous SimpleInterferogram object
        :return: Nothing returned
        """

        bodyd_id = self.track1.data.attrs.get("body_id")
        start_time1 = self.track1.data.attrs.get("start_time")
        end_time1 = self.track1.data.attrs.get("end_time")
        start_time2 = self.track2.data.attrs.get("start_time")
        end_time2 = self.track2.data.attrs.get("end_time")

        # Open the HDF5 file in read mode
        repo = git.Repo('.', search_parent_directories=True)
        data = xr.open_dataset(repo.working_tree_dir + '/files/igram_' + self.deformation_map.pyre_name + '_' +
                               str(bodyd_id) + '_' + str(start_time1) + '_' + str(end_time1) +
                               str(start_time2) + '_' + str(end_time2) + ".nc")
        self.data = data

    def create_data_array(self, phases, los_displacements):
        """

        """

        # Create the xarray Dataset
        phases_da = xr.DataArray(
            data=phases,
            dims=["points"],
            coords={
                "time": ("points", self.track1.data["time"].values),
                "x": ("points", self.track1.data["x"].values),
                "y": ("points", self.track1.data["y"].values),
                "z": ("points", self.track1.data["z"].values),
                "latitude": ("points", self.track1.data["latitude"].values),
                "longitude": ("points", self.track1.data["longitude"].values),
                "height": ("points", self.track1.data["height"].values)},
            attrs=dict(
                body_id=self.conOps.body_id,
                start_time1=self.track1.start_time,
                end_time1=self.track1.end_time,
                start_time2=self.track1.start_time,
                end_time2=self.track2.end_time,
            ),
        )

        # Create the xarray Dataset
        los_displacements_da = xr.DataArray(
            data=los_displacements,
            dims=["points"],
            coords={
                "time": ("points", self.track1.data["time"].values),
                "x": ("points", self.track1.data["x"].values),
                "y": ("points", self.track1.data["y"].values),
                "z": ("points", self.track1.data["z"].values),
                "latitude": ("points", self.track1.data["latitude"].values),
                "longitude": ("points", self.track1.data["longitude"].values),
                "height": ("points", self.track1.data["height"].values)},
            attrs=dict(
                body_id=self.conOps.body_id,
                start_time1=self.track1.start_time,
                end_time1=self.track1.end_time,
                start_time2=self.track1.start_time,
                end_time2=self.track2.end_time,
            ),
        )

        # Create a Dataset
        dataset = xr.Dataset({
            "phases": phases_da,
            "los_displacements": los_displacements_da,
        })

        self.data = dataset

    def get_flattened_angles(self, time, satellite_position, flat_positions):
        """
        Calculate the flattened angle from a given set of positions at a time in the satellite orbit
        :param flat_positions: Flattened positions
        :param time: Time at which the satellite is in satellite_position
        :param satellite_position: Position of the satellite at a specific time
        :return: Flattened angles of satellite to ground positions
        """

        # Calculate the satellite intersect and height with the shape
        intersects, satellite_heights = self.planet.get_sub_obs_points(times=time, conops=self.conOps)

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

    def get_flattened_positions(self, positions, satellite_position):

        # Get planet axes
        a, b, c = self.planet.get_axes()

        # Calculate vectors from the ground positions to the satellite
        positions = np.asarray(positions)
        ground_satellite_vectors = positions - satellite_position

        # Create plates normal to the vector from the ground to the satellite at the ground position
        planes = spice.nvp2pl_vector(normal=ground_satellite_vectors, point=positions)

        # Generate ellipses from these planes that intersect the planet ellipsoid
        ellipses = spice.inedpl_vector(a=a, b=b, c=c, plane=planes)[0]

        # Get the "flattened" points from the intersection of ellipses and ground positions
        flat_points = spice.npelpt_vector(point=positions, ellips=ellipses)[0]

        return flat_points

    def calculate_igram(self):
        """
        Calculate the flattened phases between two swaths given a baseline
        :return: Nothing returned 
        """

        # Get planet axes
        a, b, c = self.planet.get_axes()

        # For speed, calculate all points at once
        u_displacements, v_displacements, w_displacements = (
            self.deformation_map.get_displacements(track=self.track1))
        displacements_1 = np.array([u_displacements, v_displacements, w_displacements]).T

        # For speed, calculate all points at once
        u_displacements, v_displacements, w_displacements = (
            self.deformation_map.get_displacements(track=self.track2))
        displacements_2 = np.array([u_displacements, v_displacements, w_displacements]).T

        # Access x, y, z values
        x = self.track1.data["x"].values
        y = self.track1.data["y"].values
        z = self.track1.data["z"].values

        # Get the positions of the groundTargets
        positions = np.asarray([x, y, z]).T

        # Get satellite positions and velocities
        satellite_positions, sat_velocities = self.conOps.get_states(times=self.track1.data["time"].values)

        flattened_points = self.get_flattened_positions(positions=positions, satellite_position=satellite_positions)

        # Get flattened angles for the satellite positions and ground points
        angles = self.get_flattened_angles(time=self.track1.data["time"].values,
                                           satellite_position=satellite_positions, flat_positions=flattened_points)

        # Calculate geodetic heights of each of the ground points
        geodetic_heights = spice.nearpt_vector(positn=positions, a=a, b=b, c=c)[1]

        # Get the distances from the flattened points to the satellites
        vectors = positions - satellite_positions
        distances = np.asarray([np.linalg.norm(vector) for vector in vectors])

        # Calculate the vectors from the ground points to the satellite positions
        ground_satellite_vectors = satellite_positions - positions

        # Calculate the line of sight displacements given the satellite positions for the first swath
        los_displacements1 = np.asarray([(np.dot(displacement, ground_satellite_vector) /
                                          np.linalg.norm(ground_satellite_vector))
                                         for displacement, ground_satellite_vector
                                         in zip(displacements_1, ground_satellite_vectors)])

        # Calculate the line of sight displacements given the satellite positions for the second swath
        los_displacements2 = np.asanyarray([(np.dot(displacement, ground_satellite_vector) /
                                             np.linalg.norm(ground_satellite_vector))
                                            for displacement, ground_satellite_vector
                                            in zip(displacements_2, ground_satellite_vectors)])

        los_displacements = los_displacements2 - los_displacements1

        # Calculate the phases using the LOS displacements, baseline, geodetic heights, distances, angles
        phases = 4 * np.pi * (1 / self.instrument.wavelength) * (los_displacements -
                                                                 ((self.baseline * geodetic_heights) /
                                                                  (distances * np.sin(angles))))

        phases = np.fmod(phases, 2 * np.pi)
        phases = [angle + 2 * np.pi if angle < 0 else angle for angle in phases]

        self.create_data_array(phases=phases, los_displacements=los_displacements)

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

        longitudes = self.data["longitude"].values
        latitudes = self.data["latitude"].values
        phases = self.data["phases"].values

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
        ax.set_title('Interferogram', pad=20)

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/' + 'interferogram_no_baseline.png', format='png', dpi=500)

    def visualize_displacements(self, projection, fig=None, globe=None, ax=None, return_fig=False):
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

        longitudes = self.data["longitude"].values
        latitudes = self.data["latitude"].values
        los_displacements = self.data["los_displacements"].values

        im = ax.scatter(longitudes, latitudes,
                        transform=ccrs.PlateCarree(globe=globe),
                        c=los_displacements, marker='o', s=0.1)

        # return fig, ax, globe if necessary
        if return_fig:
            return fig, ax, globe

        # Add colorbar
        plt.colorbar(im, label="Phase")

        # Add labels and legend
        ax.set_title('Interferogram', pad=20)

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/' + 'displacements.png', format='png', dpi=500)

# end of file
