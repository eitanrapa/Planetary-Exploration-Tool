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
import cspyce as spice
from scipy.optimize import root
import matplotlib.pyplot as plt
from tqdm import tqdm
import sys


class SimpleInterferogramSinglePass(pet.component,
                                    family="pet.dataAnalysis.surfaceDeformationAnalysis.simpleInterferogram",
                                    implements=pet.protocols.dataAnalysis.surfaceDeformationAnalysis):
    """
    Class that creates an interferogram from a single pass given an instrument and orbit
    """

    planet = pet.protocols.planet()
    planet.doc = "target planet"

    campaign = pet.protocols.campaigns.orbiter()
    campaign.doc = "campaign (must be an orbiter)"

    instrument = pet.protocols.instruments.inSAR()
    instrument.doc = "observation instrument"

    track = pet.protocols.dataAcquisition()
    track.doc = "ground track"

    baseline = pet.properties.float()
    baseline.doc = "baseline between the two tracks [m]"

    baseline_uncertainty = pet.properties.float()
    baseline_uncertainty.doc = "baseline uncertainty [m]"

    roll = pet.properties.float()
    roll.doc = "roll angle [deg]"

    roll_uncertainty = pet.properties.float()
    roll_uncertainty.doc = "roll uncertainty [deg]"

    forward_data = None
    inverse_data = None

    @classmethod
    def from_file(cls, planet, instrument, campaign, forward_phase_file_name, inverse_phase_file_name):
        """
        Load the interferogram from an HDF5 file
        :param planet: Planet object
        :param instrument: Instrument object
        :param campaign: Campaign object
        :param forward_phase_file_name: Name of the forward phase file
        :param inverse_phase_file_name: Name of the inverse phase file (optional)
        :return: Interferogram object
        """

        # Open the HDF5 file in read mode
        datatree = xr.open_datatree(filename_or_obj=forward_phase_file_name)
        track_da = datatree["track"].to_dataset()["track"]
        data = datatree["interferogram"]["forward_phases"]

        # Check if there is an inverse phase data array
        if inverse_phase_file_name is not None:
            inverse_data = xr.open_dataarray(filename_or_obj=inverse_phase_file_name)
        else:
            inverse_data = None

        # Create the track object
        track = pet.dataAcquisition.track.from_data_array(planet=planet, campaign=campaign, instrument=instrument,
                                                          data=track_da)

        # Create the interferogram object
        obj = cls(name="igram" + str(np.random.rand()),
                  planet=planet, instrument=instrument, campaign=campaign, track=track,
                  baseline=data.attrs["baseline"], baseline_uncertainty=data.attrs["baseline_uncertainty"],
                  roll=data.attrs["roll"], roll_uncertainty=data.attrs["roll_uncertainty"])

        obj.forward_data = data # Restore computed result
        obj.inverse_data = inverse_data  # Restore computed result if available

        return obj

    @classmethod
    def from_files(cls, planet, instrument, campaign, forward_file_list, inverse_file_list=None):
        """
        Load a list of interferograms from a list of HDF5 files
        :param planet: Planet object
        :param instrument: Instrument object
        :param campaign: Campaign object
        :param forward_file_list: List of forward phase file names
        :param inverse_file_list: List of inverse phase file names (optional)
        :return: List of interferogram objects
        """

        if inverse_file_list is None:
            inverse_file_list = [None] * len(forward_file_list)

        if len(inverse_file_list) != len(forward_file_list):
            raise ValueError("The number of forward phase files must match the number of inverse phase files.")

        # Load all the files
        return [cls.from_file(planet=planet, instrument=instrument, campaign=campaign,
                              forward_phase_file_name=forward_file, inverse_phase_file_name=inverse_file) for
                forward_file, inverse_file in zip(forward_file_list, inverse_file_list)]

    def save_forward(self, file_name):
        """
        Save the track and forward calculations to an HDF5 file
        :param file_name: Name of the file to save
        :return: Nothing returned
        """

        # Open HDF5 file
        self.forward_data.to_netcdf(file_name, engine="netcdf4")

    def save_inverse(self, file_name):
        """
        Save the inverse calculations to an HDF5 file
        :param file_name: Name of the file to save
        :return: Nothing returned
        """

        # Open HDF5 file
        self.inverse_data.to_netcdf(file_name, engine="netcdf4")

    def create_forward_data_array(self, phases, corrs, nlooks, sigma_phase, sigma_height, sigma_disp,
                                  nesn, sigma0):
        """
        Create a xarray with the input data
        :param phases: forward phases of interferogram
        :param corrs: Correlation values
        :param nlooks: Number of looks
        :param sigma_phase: Standard deviation of the phase noise
        :param sigma_height: Standard deviation of the height noise
        :param sigma_disp: Standard deviation of the displacement noise
        :param nesn: Noise equivalent sigma nought
        :param sigma0: Standard deviation of the sigma nought
        :return: Nothing returned
        """

        # Create the xarray datarray
        forward_phases_da = xr.DataArray(
            data=phases,
            dims=["points"],
            coords={
                "time": ("points", self.track.data["time"].values),
                "correlation": ("points", corrs),
                "nlooks": ("points", nlooks),
                "sigma_phase": ("points", sigma_phase),
                "sigma_height": ("points", sigma_height),
                "sigma_disp": ("points", sigma_disp),
                "nesn": ("points", nesn),
                "sigma0": ("points", sigma0),
                "x": ("points", self.track.data["x"].values),
                "y": ("points", self.track.data["y"].values),
                "z": ("points", self.track.data["z"].values),
                "latitude": ("points", self.track.data["latitude"].values),
                "longitude": ("points", self.track.data["longitude"].values),
                "height": ("points", self.track.data["height"].values)},
            name="forward_phases",
            attrs=dict(
                body_id=self.campaign.body_id,
                baseline=self.baseline,
                baseline_uncertainty=self.baseline_uncertainty,
                roll=self.roll,
                roll_uncertainty=self.roll_uncertainty,
                start_time=self.track.start_time,
                end_time=self.track.end_time
            ),
        )

        # Store them in a DataTree
        dt = xr.DataTree.from_dict({
            "interferogram": xr.Dataset({"forward_phases": forward_phases_da}),
            "track": xr.Dataset({"track": self.track.data})
        })

        self.forward_data = dt

    def create_inverse_data_array(self, phases_ref, kz_k, recovered_dem):
        """
        Create a xarray with the input data
        :param phases_ref: Reference phases
        :param kz_k: Vertical interferometric wavenumber
        :param recovered_dem: Recovered DEM
        :return: Nothing returned
        """

        # Create the xarray datarray
        inverse_phases_da = xr.DataArray(
            data=phases_ref,
            dims=["points"],
            coords={
                "time": ("points", self.track.data["time"].values),
                "phase_ref": ("points", phases_ref),
                "kz_k": ("points", kz_k),
                "recovered_x": ("points", recovered_dem[:, 0]),
                "recovered_y": ("points", recovered_dem[:, 1]),
                "recovered_z": ("points", recovered_dem[:, 2]),
                "x": ("points", self.track.data["x"].values),
                "y": ("points", self.track.data["y"].values),
                "z": ("points", self.track.data["z"].values),
                "latitude": ("points", self.track.data["latitude"].values),
                "longitude": ("points", self.track.data["longitude"].values),
                "height": ("points", self.track.data["height"].values)},
            name="inverse_phases",
            attrs=dict(
                body_id=self.campaign.body_id,
                baseline=self.baseline,
                baseline_uncertainty=self.baseline_uncertainty,
                roll=self.roll,
                roll_uncertainty=self.roll_uncertainty,
                start_time=self.track.start_time,
                end_time=self.track.end_time
            ),
        )

        self.inverse_data = inverse_phases_da

    @staticmethod
    def compute_local_coordsystem(p_sat, v_sat):
        """
        Compute the local coordinate system
        :param p_sat: Satellite position
        :param v_sat: Satellite velocity
        :return: ux, uy, uz: Local coordinate system
        """

        n = len(p_sat[:, 0])
        uz = p_sat / np.sqrt(np.sum(p_sat ** 2, axis=1)).reshape([n, 1])
        uy = np.cross(v_sat / np.sqrt(np.sum(v_sat ** 2, axis=1)).reshape([n, 1]), uz, axis=1)
        ux = np.cross(uz, uy, axis=1)

        return ux, uy, uz

    def add_constant_baseline(self, p_sat, v_sat, dy0, dz0):
        """
        Add a constant baseline
        :param p_sat: Satellite position
        :param v_sat: Satellite velocity
        :param dy0: Baseline in y direction
        :param dz0: Baseline in z direction
        :return: pSat, vSat: Satellite position and velocity with the baseline added
        """

        # baseline in local coordinates
        dx = 0.
        dy = dy0
        dz = dz0
        ddx = 0.
        ddy = 0.
        ddz = 0.

        ux, uy, uz = self.compute_local_coordsystem(p_sat=p_sat, v_sat=v_sat)
        p_sat = p_sat + dx * ux + dy * uy + dz * uz
        v_sat = v_sat + ddx * ux + ddy * uy + ddz * uz
        return p_sat, v_sat

    @staticmethod
    def solve_insar(vars, phi, wavelength, p_m, p_s, v_m, rho, doppler):
        """
        Solve the InSAR equations to get geocoded dem
        :param vars: positions
        :param phi: reference phases
        :param wavelength: wavelength of the instrument
        :param p_m: position with uncertainty of the first antenna
        :param p_s: position with uncertainty of the second antenna
        :param v_m: velocity of the first antenna
        :param rho: distance to the first antenna
        :param doppler: Doppler frequency
        :return: residuals of the equations
        """

        p = np.array(vars)
        return [(4 * np.pi / wavelength * (np.linalg.norm(p - p_s) - np.linalg.norm(p - p_m))) - phi,
                np.linalg.norm(p - p_m) - rho,
                (-2 / wavelength * np.dot(v_m, (p - p_m)) / np.linalg.norm((p - p_m))) - doppler]

    @staticmethod
    def solve_geocode_ellipsoid(vars, wavelength, p_m, v_m, rho, doppler, a, b, c):
        """
        Solve the geocode equations for an ellipsoid
        :param vars: positions
        :param wavelength: wavelength of the instrument
        :param p_m: position with uncertainty of the first antenna
        :param v_m: velocity of the first antenna
        :param rho: distance to the first antenna
        :param doppler: Doppler frequency
        :param a: semi-major axis of the ellipsoid
        :param b: semi-minor axis of the ellipsoid
        :param c: semi-minor axis of the ellipsoid
        :return: residuals of the equations
        """

        p = np.array(vars)
        return [(p[0] ** 2 / a ** 2 + p[1] ** 2 / b ** 2 + p[2] ** 2 / c ** 2) * 2e3 - 2e3,
                np.linalg.norm(p - p_m) - rho,
                (-2 / wavelength * np.dot(v_m, (p - p_m)) / np.linalg.norm((p - p_m))) * 2e3 - doppler * 2e3]

    def solve_geocode_dem(self, vars, wavelength, p_m, v_m, rho, doppler, dem_file):
        """
        Solve the geocode equations for a DEM
        :param vars: positions
        :param wavelength: wavelength of the instrument
        :param p_m: position with uncertainty of the first antenna
        :param v_m: velocity of the first antenna
        :param rho: distance to the first antenna
        :param doppler: Doppler frequency
        :param dem_file: DEM file to use for the geocoding
        """

        p = np.array(vars)
        return [self.get_dem_distance_from_surface(dem_file=dem_file, point=p),
                np.linalg.norm(p - p_m) - rho,
                (-2 / wavelength * np.dot(v_m, (p - p_m)) / np.linalg.norm((p - p_m))) * 2e3 - doppler * 2e3]

    @staticmethod
    def get_dem_incidence_angles(dem_file, satellite_position, los):
        """
        Get the incidence angles for an arbitrary DEM file
        :param dem_file: DEM file to use
        :param satellite_position: Position of the satellite
        :param los: Lines of sight from the satellite to the ground targets
        :return: Incidence angles in degrees
        """

        # Create a file manager
        fm = pet.spiceTools.fileManager(folder_path="")

        # Furnish the DEM file
        fm.furnsh(names_list=[dem_file])

        # Retrieve the DSK handle for a second loaded dsk file
        handle = spice.kdata(which=1, kind="dsk")[3]

        # Retrieve the DSK DLA
        dla = spice.dlabfs(handle=handle)

        total_plates = []

        # Use the SPICE toolkit to calculate the intersects
        for raydir in los:
            # Check if intersect exists
            try:
                plate_ids, _ = spice.dskx02(handle=handle, dladsc=dla,
                                                     vertex=satellite_position*1e-3, raydir=raydir)[:2]
                total_plates.append(plate_ids)
            except ValueError:
                # If not, do nothing
                pass

        # Calculate incidence angles
        normals = [spice.dskn02(handle, dla, plate_id) for plate_id in total_plates]
        normals = np.asarray(normals)

        # Project normal vector on incidence plane
        line_of_sight_norms = los / np.linalg.norm(los, axis=-1)[:, np.newaxis]

        # project normal vector on incidence plane
        plane_normal_vectors = np.cross(-satellite_position, line_of_sight_norms)
        plane_normal_vectors = plane_normal_vectors / np.linalg.norm(plane_normal_vectors, axis=-1)[:, np.newaxis]
        normal_perps = np.sum(normals * plane_normal_vectors, axis=-1)[:, np.newaxis] * plane_normal_vectors
        normal_projs = normals - normal_perps
        normal_projs = normal_projs / np.linalg.norm(plane_normal_vectors, axis=-1)[:, np.newaxis]

        # get local incident angle
        incidence_angles = np.arccos(np.sum(-line_of_sight_norms * normal_projs, axis=1))
        incidence_angles = np.degrees(incidence_angles)

        # Return the intersects, incidence angles, and look angles
        return incidence_angles

    @staticmethod
    def get_dem_distance_from_surface(dem_file, point):
        """
        Get the distance from a point to the surface of a DEM file
        :param dem_file: DEM file to use
        :param point: Point in Cartesian coordinates (x, y, z) in meters
        :return: Distance from the point to the surface of the DEM file in meters
        """

        # Create a file manager
        fm = pet.spiceTools.fileManager(folder_path="")

        # Furnish the DEM file
        fm.furnsh(names_list=[dem_file])

        # Retrieve the DSK handle for a second loaded dsk file
        handle = spice.kdata(which=1, kind="dsk")[3]

        # Retrieve the DSK DLA
        dla = spice.dlabfs(handle=handle)

        # # Use the SPICE toolkit to calculate the intersects
        intersect = spice.dskx02(handle=handle, dladsc=dla, vertex=point*1e-3*1.1, raydir= -1*point)[1]

        # Convert to meters
        intersect = np.asanyarray(intersect) * 1e3

        return np.linalg.norm(intersect - point)

    def get_satellite_positions(self):
        """
        Get the satellite positions and velocities for the two antennas
        :return: satellite_positions1, satellite_velocity1, satellite_positions2, satellite_velocity2,
                 positions, lines_of_sight, f_doppler
        """

        print("Starting the interferogram calculation for the topography...", file=sys.stderr)

        # Access x, y, z values
        x = self.track.data["x"].values
        y = self.track.data["y"].values
        z = self.track.data["z"].values

        # Get the positions of the groundTargets
        positions = np.asarray([x, y, z]).T

        print("     Getting satellite positions of orbit 1...", file=sys.stderr)

        # Get satellite positions and velocities of first orbit
        satellite_positions1, satellite_velocity1 = \
            self.campaign.get_states(times=self.track.data["time"].values)

        # Get lines of sight and Doppler frequency
        lines_of_sight = positions - satellite_positions1
        f_doppler = -2 * np.array([np.dot(v1, v2) for v1, v2 in zip(satellite_velocity1, lines_of_sight)]) \
                    / self.instrument.wavelength / np.linalg.norm(lines_of_sight, axis=1)

        # Make roll angle and uncertainty in radians
        roll = np.deg2rad(self.roll)

        print("     Getting satellite positions of orbit 2 by applying offset...", file=sys.stderr)
        base_rad = self.baseline * np.sin(roll)
        base_hor = self.baseline * np.cos(roll)

        # Get satellite positions and velocities of second orbit
        satellite_positions2, satellite_velocity2 = self.add_constant_baseline(p_sat=satellite_positions1,
                                                                               v_sat=satellite_velocity1,
                                                                               dy0=base_hor, dz0=base_rad)

        return (satellite_positions1, satellite_velocity1, satellite_positions2, satellite_velocity2,
                positions, lines_of_sight, f_doppler)

    def get_satellite_positions_knowledge(self, satellite_positions1, satellite_velocity1):
        """
        Get the satellite positions and velocities for the two antennas with knowledge uncertainties
        :param satellite_positions1: Satellite positions of the first orbit
        :param satellite_velocity1: Satellite velocities of the first orbit
        :return: satellite_positions_k1, satellite_velocity_k1, satellite_positions_k2, satellite_velocity_k2
        """

        print("     Getting satellite position knowledge of the orbits...", file=sys.stderr)

        # Make roll angle and uncertainty in radians
        roll = np.deg2rad(self.roll)
        roll_uncertainty = np.deg2rad(self.roll_uncertainty)

        # Get satellite positions and velocity knowledge of first orbit
        satellite_positions_k1, satellite_velocity_k1 = satellite_positions1, satellite_velocity1

        base_rad_k = (self.baseline + self.baseline_uncertainty) * np.sin(roll + roll_uncertainty)
        base_hor_k = (self.baseline + self.baseline_uncertainty) * np.cos(roll + roll_uncertainty)

        # Get satellite positions and velocity knowledge of second orbit
        satellite_positions_k2, satellite_velocity_k2 = self.add_constant_baseline(p_sat=satellite_positions1,
                                                                                 v_sat=satellite_velocity1,
                                                                                 dy0=base_hor_k, dz0=base_rad_k)

        return satellite_positions_k1, satellite_velocity_k1, satellite_positions_k2, satellite_velocity_k2

    @staticmethod
    def get_perpendicular_baseline(satellite_positions1, satellite_positions2, lines_of_sight):
        """
        Get the perpendicular baseline between two satellite positions
        :param satellite_positions1: Satellite positions of the first antenna
        :param satellite_positions2: Satellite positions of the second antenna
        :param lines_of_sight: Lines of sight from the satellite to the ground targets
        :return: b_perp: Perpendicular baseline
        """

        b_vec = satellite_positions2 - satellite_positions1

        lines_of_sight = lines_of_sight / np.linalg.norm(lines_of_sight, axis=-1)[:, np.newaxis]
        b_perp_vec = b_vec - np.sum(b_vec * lines_of_sight, axis=-1)[:, np.newaxis] * lines_of_sight
        b_perp = np.linalg.norm(b_perp_vec, axis=-1)

        return b_perp


    def calculate_igram(self):
        """
        Calculate the unwrapped phases between given a baseline
        :return: Nothing returned
        """

        (satellite_positions1, satellite_velocity1, satellite_positions2, satellite_velocity2,
        positions, lines_of_sight, _) = self.get_satellite_positions()

        print("     Forward computation of the interferometric phase of the topography...", file=sys.stderr)

        # Get distances and topo phase for forward simulation
        distances1 = np.linalg.norm(positions - satellite_positions1, axis=-1)
        distances2 = np.linalg.norm(positions - satellite_positions2, axis=-1)

        phases_forward = 4 * np.pi * (1 / self.instrument.wavelength) * (distances2 - distances1)

        # compute performance and noise-----------------------------------------
        print("     Get performance and noise...", file=sys.stderr)

        # get b_perp
        b_perp = self.get_perpendicular_baseline(satellite_positions1, satellite_positions2, lines_of_sight)

        # get performance
        sigma_phase, sigma_height, sigma_disp, corr_tot, n_looks, nesn, sigma0 = self.instrument.get_instrument_noise(
            planet=self.planet,
            perpendicular_baseline=b_perp,
            satellite_velocities=satellite_velocity1,
            look_angles=self.track.data.values,
            incidence_angles=self.track.data["incidence_angle"].values,
            non_projected_incidence_angles=self.track.data["non_projected_incidence_angle"].values,
            distances=distances1, variable_backscatter=True)

        noise = np.random.normal(loc=0, scale=sigma_phase, size=len(phases_forward))

        print("        Noise mean meas:" + str(np.mean(noise)) + "...", file=sys.stderr)
        print("        Noise std meas:" + str(np.std(noise)) + "...", file=sys.stderr)

        self.create_forward_data_array(phases=phases_forward,
                                       sigma_phase=sigma_phase, sigma_height=sigma_height, sigma_disp=sigma_disp,
                                       corrs=corr_tot, nlooks=n_looks, nesn=nesn, sigma0=sigma0)

    def get_igram_inversion(self, use_dem="exact", dem_file=None):
        """
        Calculate the interferogram inversion to get the topography and reference phases
        :param use_dem: Method to use for the DEM ('exact', 'reference', 'ellipsoid')
        :param dem_file: DEM file to use for the geocoding (if use_dem is 'reference')
        :return: Nothing returned
        """

        (satellite_positions1, satellite_velocity1, _, _,
        positions, lines_of_sight, f_doppler) = self.get_satellite_positions()

        distances1 = np.linalg.norm(positions - satellite_positions1, axis=-1)

        satellite_positions_k1, satellite_velocity_k1, satellite_positions_k2, satellite_velocity_k2 = \
        self.get_satellite_positions_knowledge(satellite_positions1=satellite_positions1,
                                               satellite_velocity1=satellite_velocity1)

        phases_forward = self.forward_data["interferogram"]["forward_phases"].values

        print("     Start InSAR inversion...", file=sys.stderr)
        print("     .....................", file=sys.stderr)

        print("        Getting reference phase...", file=sys.stderr)

        # Get distances and phase for inversion (i.e., reference phase)
        if use_dem == 'exact':

            print("          Using exact DEM...", file=sys.stderr)

            distances1_ref = np.linalg.norm(positions - satellite_positions_k1, axis=-1)
            distances2_ref = np.linalg.norm(positions - satellite_positions_k2, axis=-1)

            phases_ref = 4 * np.pi * (1 / self.instrument.wavelength) * (distances2_ref - distances1_ref)

        elif use_dem == 'reference':

            print("          Using reference DEM...", file=sys.stderr)
            print("            Getting points on reference DEM (aka backgeo)...", file=sys.stderr)

            positions_reference = np.zeros(positions.shape)

            for ii in tqdm(range(len(distances1))):
                res = root(self.solve_geocode_dem, positions[ii, :] + np.random.normal(0, 1e3, 3),
                           args=(self.instrument.wavelength, satellite_positions1[ii, :],
                                 satellite_velocity_k1[ii, :], distances1[ii], f_doppler[ii], dem_file))
                positions_reference[ii, :] = res.x

            distances1_ref = np.linalg.norm(positions_reference - satellite_positions_k1, axis=-1)
            distances2_ref = np.linalg.norm(positions_reference - satellite_positions_k2, axis=-1)

            phases_ref = 4 * np.pi * (1 / self.instrument.wavelength) * (distances2_ref - distances1_ref)

        elif use_dem == 'ellipsoid':

            print("          Using ellipsoid...", file=sys.stderr)
            print("            Getting points on ellipsoid (aka backgeo)...", file=sys.stderr)

            positions_ellipsoid = np.zeros(positions.shape)
            a, b, c = self.planet.get_axes()
            for ii in tqdm(range(len(distances1))):
                res = root(self.solve_geocode_ellipsoid, positions[ii, :],
                           args=(self.instrument.wavelength, satellite_positions1[ii, :],
                                 satellite_velocity_k1[ii, :], distances1[ii], f_doppler[ii],
                                 a, b, c))
                positions_ellipsoid[ii, :] = res.x
            distances1_ref = np.linalg.norm(positions_ellipsoid - satellite_positions_k1, axis=-1)
            distances2_ref = np.linalg.norm(positions_ellipsoid - satellite_positions_k2, axis=-1)

            phases_ref = 4 * np.pi * (1 / self.instrument.wavelength) * (distances2_ref - distances1_ref)

        else:
            raise ValueError("useDEM must be 'exact', 'reference' or 'ellipsoid'!")

        print("        Do exact InSAR height computation...", file=sys.stderr)

        # do InSAR Geocoding-------------------------------
        recovered_dem = np.zeros(positions.shape)
        for ii in tqdm(range(len(phases_forward))):
            # define variables
            res = root(self.solve_insar, positions[ii, :], args=(phases_ref[ii],
                                                           self.instrument.wavelength, satellite_positions_k1[ii, :],
                                                           satellite_positions_k2[ii, :], satellite_velocity_k1[ii, :],
                                                           distances1[ii], f_doppler[ii]))
            recovered_dem[ii, :] = res.x

        # -------------------------------------------------
        print("     .....................", file=sys.stderr)
        print("     Finished InSAR inversion!", file=sys.stderr)

        print("     Getting some relevant InSAR parameters (based on orbit and topo knowledge)...", file=sys.stderr)

        # get perpendicular baseline
        b_vec = satellite_positions_k2 - satellite_positions_k1
        los_k = positions - satellite_positions_k1
        los_k = los_k / np.linalg.norm(los_k, axis=-1)[:, np.newaxis]
        b_perp_vec = b_vec - np.sum(b_vec * los_k, axis=-1)[:, np.newaxis] * los_k
        b_perp_k = np.linalg.norm(b_perp_vec, axis=-1)

        if use_dem=='exact':

            print("        Using exact DEM...", file=sys.stderr)

            # get incident angle (just read from track here)
            inc_ang_k = self.track.data["incidence_angle"].values

        elif use_dem == 'reference':

            print("        Using reference DEM...", file=sys.stderr)

            # get incident angle
            inc_ang_k = self.get_dem_incidence_angles(dem_file=dem_file, satellite_position=satellite_positions_k1,
                                                      los=los_k)

        elif use_dem == 'ellipsoid':

            print("        Assuming ellipsoid...", file=sys.stderr)

            # get incident angle
            normals = np.asanyarray([spice.surfnm(a, b, c, pos) for pos in positions_ellipsoid])
            inc_ang_k = np.arccos(np.sum(normals * (-los_k), axis=-1))

        # get vertical interferometric wavenumber kz
        kz_k = 4 * np.pi / self.instrument.wavelength * b_perp_k / distances1_ref / np.sin(inc_ang_k)

        self.create_inverse_data_array(phases_ref=phases_ref,
                                       kz_k=kz_k, recovered_dem=recovered_dem)

        print("Finished the interferogram calculations!", file=sys.stderr)

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
        longitudes = self.forward_data["longitude"].values
        latitudes = self.forward_data["latitude"].values
        forward_phases = self.forward_data.values
        inverse_phases = self.inverse_data.values

        phases = forward_phases - inverse_phases

        # Wrap interferogram
        phases = np.mod(phases, 2 * np.pi)

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
        plt.savefig(fname=projection.folder_path + '/' + 'interferogram_single_pass' + '_' +
                    str(self.campaign.body_id) +
                    '_' + str(self.track.start_time) + '_' + str(
                    self.track.end_time) + '.png', format='png',
                    dpi=500)

# end of file
