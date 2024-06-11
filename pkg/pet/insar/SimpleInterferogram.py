#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved
import pet
import numpy as np
import cspyce as spice
import snaphu
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from ..ext import conversions
from ..projections import plottingTools
import numpy as np
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
import os


class SimpleInterferogram(pet.component):
    """

    """
    pairing_one = pet.properties.int()
    pairing_one.doc = ""

    pairing_two = pet.properties.int()
    pairing_two.doc = ""

    track_number = pet.properties.int()
    track_number.doc = ""

    def __init__(self, name, locator, implicit):
        super().__init__(name, locator, implicit)
        self.igram = []

    def get_angles_cartesian(self, instrument, time, planet, satellite_position, flat_positions):
        """
        """

        # Calculate the satellite intersect and height with the shape
        intersect, satellite_height = planet.get_sub_obs_points(times=time, instrument=instrument)

        # Calculate the distance between center and satellite intersects
        satellite_radius = np.asarray(np.linalg.norm(intersect - [0, 0, 0]))

        # Get the distance between the satellite and the surface points
        distances = np.asarray([np.linalg.norm(flat_position - satellite_position)
                                for flat_position in flat_positions])

        # Calculate cosines of the angles
        z_plus_re = satellite_height + satellite_radius
        altitudes = np.asarray([np.linalg.norm(flat_position - [0, 0, 0])
                                for flat_position in flat_positions])
        cosine_look_angle = (distances ** 2 + z_plus_re ** 2 - altitudes ** 2) / (2 * distances * z_plus_re)

        # Use arc cosine to get the angles in radians
        look_angles_radians = np.arccos(cosine_look_angle)

        # Return the look angles 2-D array
        return look_angles_radians

    def get_flattened_angles(self, positions, instrument, time, satellite_position, planet):
        """

        """

        # Get planet axes
        a, b, c = planet.get_axes()

        positions = np.asarray(positions)
        ground_satellite_vectors = positions - satellite_position

        planes = spice.nvp2pl_vector(normal=ground_satellite_vectors, point=positions)

        ellipses = spice.inedpl_vector(a=a, b=b, c=c, plane=planes)[0]

        flat_points = spice.npelpt_vector(point=positions, ellips=ellipses)[0]

        flat_angles = self.get_angles_cartesian(instrument=instrument, time=time, planet=planet,
                                                satellite_position=satellite_position, flat_positions=flat_points)

        return flat_angles

    def make_raster(self, planet, phases, total_coordinates):
        """

        """
        planet_axes = planet.get_axes()

        # Define Enceladus globe
        img_globe = ccrs.Globe(semimajor_axis=planet_axes[0] * 1e3, semiminor_axis=planet_axes[2] * 1e3, ellipse=None)

        # Create a circular map using Cylindrical projection
        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Stereographic(central_latitude=-90, globe=img_globe)})

        ax.set_extent([-180, -60, -90, 90], crs=ccrs.PlateCarree(globe=img_globe))

        for i in range(len(phases)):
            coordinates = [(long, lat, height * 1e3) for lat, long, height in total_coordinates[i]]
            longitudes, latitudes, heights = zip(*coordinates)
            # Plot points on the map
            ax.scatter(longitudes, latitudes, transform=ccrs.PlateCarree(globe=img_globe),
                       c=phases[i], marker='o', s=0.1, alpha=0.25)

        # Add latitude and longitude lines
        gl = ax.gridlines(crs=ccrs.PlateCarree(globe=img_globe), linewidth=0.5, color='black', alpha=0.5,
                          linestyle='--', draw_labels=True)
        gl.top_labels = True
        gl.left_labels = True
        gl.right_labels = True
        gl.xlines = True
        gl.xlocator = mticker.FixedLocator([-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180])
        gl.ylocator = mticker.FixedLocator(
            [-90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 8, 'color': 'gray'}
        gl.ylabel_style = {'size': 8, 'color': 'grey'}

        # Get the current working directory
        current_dir = os.getcwd()

        # Construct the path three directories up
        path_three_dirs_up = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir, os.pardir))

        # Define the target directory within the three-up directory
        target_dir = os.path.join(path_three_dirs_up, 'figs')

        # Define the full path to save the plot
        full_path = os.path.join(target_dir, 'test.png')

        # Show the plot
        plt.savefig(fname=full_path, format='png', dpi=500)

        # grid_x, grid_y, grid_z = np.mgrid[np.min(positions[:, 0]):np.max(positions[:, 0]):200j,
        #                          np.min(positions[:, 1]):np.max(positions[:, 1]):200j,
        #                          np.min(positions[:, 2]):np.max(positions[:, 2]):200j]
        # phase_grid = griddata(positions, phases, (grid_x, grid_y, grid_z), method='cubic')
        # self.igram = phase_grid

    def calculate_flattened_phases(self, instrument, planet, ground_swath1, ground_swath2, baseline):
        """

        """

        a, b, c = planet.get_axes()

        total_coordinates = []
        total_phases = []
        for i in range(len(ground_swath1.swath_beams)):
            positions = np.asarray([point.get_position() for point in ground_swath1.swath_beams[i]])
            seen_values = np.asarray([point.seen for point in ground_swath1.swath_beams[i]])
            geodetic_coordinates = np.asarray(conversions.geodetic(planet=planet,
                                                                   cartesian_coordinates=positions))[:, :3]
            satellite_position, sat_velocity = instrument.get_states(planet=planet, times=ground_swath1.time_space[i])
            angles = self.get_flattened_angles(positions=positions, instrument=instrument,
                                               time=ground_swath1.time_space[i],
                                               satellite_position=satellite_position, planet=planet)
            geodetic_heights = spice.nearpt_vector(positn=positions, a=a, b=b, c=c)[1]
            distances = np.asarray([np.linalg.norm(position - [0, 0, 0]) for position in positions])

            ground_satellite_vectors = satellite_position - positions
            displacements_1 = [point.displacement for point in ground_swath1.swath_beams[i]]
            los_displacements1 = (
                np.asarray([np.dot(displacement, ground_satellite_vector) / np.linalg.norm(ground_satellite_vector)
                            for displacement, ground_satellite_vector
                            in zip(displacements_1, ground_satellite_vectors)]))

            displacements_2 = [point.displacement for point in ground_swath2.swath_beams[i]]
            los_displacements2 = (
                np.asarray([np.dot(displacement, ground_satellite_vector) / np.linalg.norm(ground_satellite_vector)
                            for displacement, ground_satellite_vector
                            in zip(displacements_2, ground_satellite_vectors)]))

            phases = (4 * np.pi * (1 / instrument.wavelength) *
                      ((los_displacements2 - los_displacements1) -
                       (np.random.normal(loc=0, scale=baseline) * geodetic_heights * 1e3 /
                        (distances * 1e3 * np.sin(angles)))))
            for k in range(len(phases)):
                if not seen_values[k]:
                    phases[k] = np.nan
            total_phases.append(phases)
            total_coordinates.append(geodetic_coordinates)

        self.make_raster(phases=total_phases, planet=planet, total_coordinates=total_coordinates)

    def calculate_igram(self, planet, instrument, displacements, time_interval=10, ground_resoluton=0.2, baseline=10):
        """

        """

        times = instrument.get_five_tracks(planet)

        orbit_cycle_time = instrument.orbit_cycle

        gs1 = pet.insar.groundSwath(name="swath_1", start_time=times[self.track_number - 1],
                                    end_time=times[self.track_number], time_interval=time_interval,
                                    ground_resolution=ground_resoluton)

        gs1.calculate_swath(instrument=instrument, planet=planet)

        gs2 = pet.insar.groundSwath(name="swath_2", start_time=times[self.track_number - 1],
                                    end_time=times[self.track_number], time_interval=time_interval,
                                    ground_resolution=ground_resoluton)

        gs2.swath_beams = gs1.swath_beams

        gs2.time_space = gs1.time_space

        displacements.attach(swath=gs1, planet=planet, use_mid_point=True)

        displacements.attach(swath=gs2, planet=planet, use_mid_point=True, time_displacement=orbit_cycle_time)

        self.calculate_flattened_phases(instrument=instrument, planet=planet, ground_swath1=gs1, ground_swath2=gs2,
                                        baseline=baseline)

    def unwrap(self):
        # Sample coherence for an interferogram with no noise.
        corr = np.ones(self.igram.shape, dtype=np.float32)

        # Unwrap using the 'SMOOTH' cost mode and 'MCF' initialization method.
        unw, conncomp = snaphu.unwrap(self.igram.shape, corr, nlooks=1.0, cost="smooth", init="mcf")

        return unw, conncomp

    def visualize(self):
        """

        """
        # Convert the irregular list of lists into a regular grid
        plt.figure(figsize=(8, 6))
        plt.imshow(self.igram, cmap='hsv', edgecolors='k')
        plt.colorbar(label='Phase (radians)')
        plt.show()

# end of file
