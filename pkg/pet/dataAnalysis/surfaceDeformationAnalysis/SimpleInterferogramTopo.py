#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import pet
import xarray as xr
import cspyce as spice
import numpy as np
import cartopy.crs as ccrs
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.optimize import root
import sys
import pdb


def computeLocalCoordSystem(pSat, vSat):
    N = len(pSat[:, 0])
    uz = pSat / np.sqrt(np.sum(pSat**2, axis = 1)).reshape([N, 1])
    uy = np.cross(vSat / np.sqrt(np.sum(vSat**2, axis = 1)).reshape([N, 1]), uz, axis = 1)
    ux = np.cross(uz, uy, axis = 1)
    return ux, uy, uz
    
def addConstantBaseline(pSat, vSat, dy0, dz0):
    " adds a constant baseline "
    # baseline in local coordinates
    Dx  = 0.
    Dy  = dy0
    Dz  = dz0
    dDx = 0.
    dDy = 0.
    dDz = 0.

    ux, uy, uz = computeLocalCoordSystem(pSat, vSat)
    pSat = pSat + Dx * ux +  Dy * uy +  Dz * uz
    vSat = vSat + dDx * ux + dDy * uy + dDz * uz
    return pSat, vSat

def solveInSAR(vars, phi, wl, p_m, p_s, v_m, rho, doppler):
    p = np.array(vars)
    return[(4*np.pi/wl * (np.linalg.norm(p-p_s) - np.linalg.norm(p-p_m))) - phi,
           np.linalg.norm(p-p_m) - rho,
           (-2/wl * np.dot(v_m,(p-p_m)) / np.linalg.norm((p-p_m))) - doppler]

def solveGeocodeEllipsoid(vars, wl, p_m, v_m, rho, doppler, a, b, c):
    p = np.array(vars)
    return[(p[0]**2 / a**2 + p[1]**2 / b**2 + p[2]**2 / c**2)*2e3 - 2e3,
           np.linalg.norm(p-p_m) - rho,
           (-2/wl * np.dot(v_m,(p-p_m)) / np.linalg.norm((p-p_m)))*2e3 - doppler*2e3]

class SimpleInterferogramTopo(pet.component, family="pet.dataAnalysis.surfaceDeformationAnalysis.simpleInterferogram",
                          implements=pet.protocols.dataAnalysis.surfaceDeformationAnalysis):
    """
    Class that creates a single interferogram between two points of a displacement map given an instrument orbit
    """

    planet = pet.protocols.planet()
    planet.doc = "target planet"

    campaign = pet.protocols.campaigns.orbiter()
    campaign.doc = "campaign (must be an orbiter)"

    instrument = pet.protocols.instruments.inSAR()
    instrument.doc = "observation instrument"

    track = pet.protocols.dataAcquisition()
    track.doc = "track"

    baseline = pet.properties.float()
    baseline.doc = "baseline between the two tracks [m]"

    baseline_uncertainty = pet.properties.float()
    baseline_uncertainty.doc = "baseline uncertainty [m]"

    roll = pet.properties.float()
    roll.doc = "roll angle [rad]"

    roll_uncertainty = pet.properties.float()
    roll_uncertainty.doc = "roll uncertainty [rad]"

    data = None

    def save(self, file_name):
        """
        Save the track to an HDF5 file
        :param file_name: Name of the file to save
        :return: Nothing returned
        """

        # Open HDF5 file
        self.data.to_netcdf(file_name, engine="netcdf4")
    
    def create_data_array(self, phase_absolut, phase_ref, kz, incAng):
        """
        Create a xarray with the input data
        :param los_displacements: displacement values measured in the LOS
        :return: Nothing returned
        """

        # Create the xarray datarray
        phase_absolut = xr.DataArray(
            data=phase_absolut,
            dims=["points"],
            coords={
                "sat_pos_time": ("points", self.track.data["sat_pos_time"].values),
                "time1": ("points", self.track.data["time"].values),
                "time2": ("points", self.track.data["time"].values),
                "x": ("points", self.track.data["x"].values),
                "y": ("points", self.track.data["y"].values),
                "z": ("points", self.track.data["z"].values),
                "latitude": ("points", self.track.data["latitude"].values),
                "longitude": ("points", self.track.data["longitude"].values),
                "height": ("points", self.track.data["height"].values)},
            name="interferogram",
            attrs=dict(
                body_id=self.campaign.body_id,
                baseline=self.baseline,
                baseline_uncertainty=self.baseline_uncertainty,
                roll=self.roll,
                roll_uncertainty=self.roll_uncertainty,
                start_time=self.track.start_time,
                end_time=self.track.end_time,
            ),
        )

        phase_ref = xr.DataArray(
            data=phase_ref,
            dims=["points"],
            coords={
                "sat_pos_time": ("points", self.track.data["sat_pos_time"].values),
                "time1": ("points", self.track.data["time"].values),
                "time2": ("points", self.track.data["time"].values),
                "x": ("points", self.track.data["x"].values),
                "y": ("points", self.track.data["y"].values),
                "z": ("points", self.track.data["z"].values),
                "latitude": ("points", self.track.data["latitude"].values),
                "longitude": ("points", self.track.data["longitude"].values),
                "height": ("points", self.track.data["height"].values)},
            name="interferogram",
            attrs=dict(
                body_id=self.campaign.body_id,
                baseline=self.baseline,
                baseline_uncertainty=self.baseline_uncertainty,
                roll=self.roll,
                roll_uncertainty=self.roll_uncertainty,
                start_time=self.track.start_time,
                end_time=self.track.end_time,
            ),
        )

        kz = xr.DataArray(
            data=kz,
            dims=["points"],
            coords={
                "sat_pos_time": ("points", self.track.data["sat_pos_time"].values),
                "time1": ("points", self.track.data["time"].values),
                "time2": ("points", self.track.data["time"].values),
                "x": ("points", self.track.data["x"].values),
                "y": ("points", self.track.data["y"].values),
                "z": ("points", self.track.data["z"].values),
                "latitude": ("points", self.track.data["latitude"].values),
                "longitude": ("points", self.track.data["longitude"].values),
                "height": ("points", self.track.data["height"].values)},
            name="interferogram",
            attrs=dict(
                body_id=self.campaign.body_id,
                baseline=self.baseline,
                baseline_uncertainty=self.baseline_uncertainty,
                roll=self.roll,
                roll_uncertainty=self.roll_uncertainty,
                start_time=self.track.start_time,
                end_time=self.track.end_time,
            ),
        )

        incAng = xr.DataArray(
            data=incAng,
            dims=["points"],
            coords={
                "sat_pos_time": ("points", self.track.data["sat_pos_time"].values),
                "time1": ("points", self.track.data["time"].values),
                "time2": ("points", self.track.data["time"].values),
                "x": ("points", self.track.data["x"].values),
                "y": ("points", self.track.data["y"].values),
                "z": ("points", self.track.data["z"].values),
                "latitude": ("points", self.track.data["latitude"].values),
                "longitude": ("points", self.track.data["longitude"].values),
                "height": ("points", self.track.data["height"].values)},
            name="interferogram",
            attrs=dict(
                body_id=self.campaign.body_id,
                baseline=self.baseline,
                baseline_uncertainty=self.baseline_uncertainty,
                roll=self.roll,
                roll_uncertainty=self.roll_uncertainty,
                start_time=self.track.start_time,
                end_time=self.track.end_time,
            ),
        )

        # Store them in a DataTree
        dt = xr.DataTree.from_dict({
            "phase_absolut": xr.Dataset({"phase_absolut": phase_absolut}),
            "phase_ref": xr.Dataset({"phase_ref": phase_ref}),
            "kz": xr.Dataset({"kz": kz}),
            "incAng": xr.Dataset({"incAng": incAng}),
            "track": xr.Dataset({"track": self.track.data}),
        })
        self.data = dt

    
    def calculate_igram(self):
        """
        Calculate the flattened phases between two swaths given a baseline
        :return: Nothing returned 
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
        satellite_positions1, satellite_velocity1 = self.campaign.get_states(times=self.track.data["sat_pos_time"].values)
        # Get lines of sight and Doppler frequency
        LoS = positions - satellite_positions1
        fDoppler = -2 * np.array([np.dot(v1, v2)for v1, v2 in zip(satellite_velocity1, LoS)])\
                     / self.instrument.wavelength / np.linalg.norm(LoS, axis=1)
        
        print("     Getting satellite positions of orbit 2 by applying offset...", file=sys.stderr)
        base_rad = self.baseline * np.sin(self.roll)
        base_hor = self.baseline * np.cos(self.roll)
        # Get satellite positions and velocities of second orbit
        satellite_positions2, satellite_velocity2 = addConstantBaseline(satellite_positions1,
                                                                             satellite_velocity1,
                                                                             base_hor, base_rad)
        
        print("     Getting satellite position knowledge of the orbits...", file=sys.stderr)
        # Get satellite positions and velocity knowledge of first orbit
        satellite_positionsK1, satellite_velocityK1 = satellite_positions1, satellite_velocity1

        base_rad_K = (self.baseline + self.baseline_uncertainty) * np.sin(self.roll + self.roll_uncertainty)
        base_hor_K = (self.baseline + self.baseline_uncertainty) * np.cos(self.roll + self.roll_uncertainty)
        # Get satellite positions and velocity knowledge of second orbit
        satellite_positionsK2, satellite_velocityK2 = addConstantBaseline(satellite_positions1,
                                                                             satellite_velocity1,
                                                                             base_hor_K, base_rad_K)

        print("     Forward computation of the interferometric phase of the topography...", file=sys.stderr)
        # Get distances and topo phase for forward simulation 
        distances1 = np.linalg.norm(positions - satellite_positions1, axis=-1)
        distances2 = np.linalg.norm(positions - satellite_positions2, axis=-1)
        phase_topo_forward = 4 * np.pi * (1 / self.instrument.wavelength) * (distances2 - distances1)

        #compute performance and noise-----------------------------------------
        print("     Get performance and noise...", file=sys.stderr)
        # get b_perp
        b_vec = satellite_positions2 - satellite_positions1
        LoS = positions - satellite_positions1
        LoS = LoS / np.linalg.norm(LoS, axis=-1)[:,np.newaxis]
        b_perp_vec = b_vec - np.sum(b_vec*LoS, axis=-1)[:,np.newaxis] * LoS
        b_perp = np.linalg.norm(b_perp_vec, axis=-1)
        # get look angle
        nadir_norm = -satellite_positions1/np.linalg.norm(satellite_positions1,axis=-1)[:,np.newaxis]
        look_ang = np.arccos(np.sum(LoS * nadir_norm, axis=-1))
        #get performance
        self.sigma_phase, self.corr_tot, self.Nlooks, self.NESN = self.instrument.get_instrument_noise(planet=self.planet,
                                             baseline=b_perp,
                                             satellite_velocity=satellite_velocity1,
                                             look_angles=look_ang,
                                             incidence_angles=self.track.data.values,
                                             distances=distances1)
        #----------------------------------------------------------------------


        print("     Start InSAR inversion...", file=sys.stderr)
        print("     .....................", file=sys.stderr)

        print("        Getting reference phase...", file=sys.stderr)
        # Get distances and phase for inversion (i.e., reference phase)
        useDEM = False
        if useDEM:
            print("          Using exact DEM...", file=sys.stderr)
            distances1_ref = np.linalg.norm(positions - satellite_positionsK1, axis=-1)
            distances2_ref = np.linalg.norm(positions - satellite_positionsK2, axis=-1)
        else:
            print("          Using ellipsoid...", file=sys.stderr)
            print("            Getting points on ellipsoid (aka backgeo)...", file=sys.stderr)
            positionsEllipsoid = np.zeros(positions.shape)
            a, b, c = self.planet.get_axes()
            for ii in tqdm(range(len(distances1))):
                res = root(solveGeocodeEllipsoid, positions[ii,:],
                            args=(self.instrument.wavelength, satellite_positions1[ii,:],
                            satellite_velocityK1[ii,:], distances1[ii], fDoppler[ii],
                            a, b, c))
                positionsEllipsoid[ii,:] = res.x
            distances1_ref = np.linalg.norm(positionsEllipsoid - satellite_positionsK1, axis=-1)
            distances2_ref = np.linalg.norm(positionsEllipsoid - satellite_positionsK2, axis=-1)
        phase_ref = 4 * np.pi * (1 / self.instrument.wavelength) * (distances2_ref - distances1_ref)
        
        print("        Flatten phase...", file=sys.stderr)
        phase_flat = phase_topo_forward - phase_ref

        #here would come the phase unwrapping....

        phase_absolut = phase_flat + phase_ref

        print("        Do exact InSAR height computation...", file=sys.stderr)
        #do InSAR Geocoding-------------------------------
        #Does only make sense without deformation
        xyz_meas = np.zeros(positions.shape)
        for ii in tqdm(range(len(phase_topo_forward))):
            #define variables
            res = root(solveInSAR, positions[ii,:], args=(phase_absolut[ii],
                            self.instrument.wavelength, satellite_positionsK1[ii,:],
                            satellite_positionsK2[ii,:], satellite_velocityK1[ii,:],
                            distances1[ii], fDoppler[ii]))
            xyz_meas[ii,:] = res.x
        self.xyz_meas = xyz_meas
        #-------------------------------------------------
        print("     .....................", file=sys.stderr)
        print("     Finished InSAR inversion!", file=sys.stderr)

        print("     Getting some relevant InSAR parameters (based on orbit and topo knowledge)...", file=sys.stderr)
        if useDEM:
            print("        Using exact DEM...", file=sys.stderr)
            # get perpendicular baseline
            b_vec = satellite_positionsK2 - satellite_positionsK1
            LoS_k = positions - satellite_positionsK1
            LoS_k = LoS_k / np.linalg.norm(LoS_k, axis=-1)[:,np.newaxis]
            b_perp_vec = b_vec - np.sum(b_vec*LoS_k, axis=-1)[:,np.newaxis] * LoS_k
            b_perp_k = np.linalg.norm(b_perp_vec, axis=-1)
            # get incident angle (just read from track here)
            incAng_k = self.track.data["incAng"].values
        else:
            print("        Assuming ellipsoid...", file=sys.stderr)
            # get perpendicular baseline
            b_vec = satellite_positionsK2 - satellite_positionsK1
            LoS_k = positionsEllipsoid - satellite_positionsK1
            LoS_k = LoS_k / np.linalg.norm(LoS_k, axis=-1)[:,np.newaxis]
            b_perp_vec = b_vec - np.sum(b_vec*LoS_k, axis=-1)[:,np.newaxis] * LoS_k
            b_perp_k = np.linalg.norm(b_perp_vec, axis=-1)
            # get incident angle
            normals = np.asanyarray([spice.surfnm(a, b, c, pos) for pos in positionsEllipsoid])
            incAng_k = np.arccos(np.sum(normals*(-LoS_k),axis=-1))
        # get vertical interferometric wavenumber kz
        kz_k = 4*np.pi/self.instrument.wavelength * b_perp_k/distances1_ref/np.sin(incAng_k)
            
        print("     Preparation of output data...", file=sys.stderr)
        self.create_data_array(phase_topo_forward, phase_ref, kz_k, incAng_k)
        
        print("Finished the interferogram calculations!", file=sys.stderr)


    def visualize_interferogram(self, phase, projection, cmap, title, fileName, fig=None, globe=None, ax=None, return_fig=False):
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

        longitudes = self.track.data["longitude"].values
        latitudes = self.track.data["latitude"].values

        # Make the colormap cyclical
        cm = plt.cm.get_cmap(cmap)

        # Iterate through the interferogram beams
        im = ax.scatter(longitudes, latitudes, cmap=cm,
                        transform=ccrs.PlateCarree(globe=globe), c=phase, s=0.01, marker=",")
        # im = ax.pcolormesh(longitudes, latitudes, phase, cmap=cm,
        #                 transform=ccrs.PlateCarree(globe=globe))

        # return fig, ax, globe if necessary
        if return_fig:
            return fig, ax, globe
        # Add colorbar
        plt.colorbar(im, fraction=0.02, pad=0.1)

        # Add labels and legend
        plt.title(title, pad=12)

        # Save the plot
        plt.savefig(fileName, format='png',
                    dpi=300,bbox_inches="tight")

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
        los_displacements = self.data.values

        im = ax.scatter(longitudes, latitudes,
                        transform=ccrs.PlateCarree(globe=globe),
                        c=los_displacements, marker='.', s=0.0001)

        # return fig, ax, globe if necessary
        if return_fig:
            return fig, ax, globe

        # Add colorbar
        plt.colorbar(im, label="Displacement")

        # Add labels and legend
        plt.title('LOS displacements', pad=20)

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/' + 'displacements_' + self.deformation_map.pyre_name + '_' +
                    str(self.campaign.body_id) + '_' + str(self.track1.start_time) + '_' + str(
                        self.track1.end_time) +
                    '_' + str(self.track2.start_time) + '_' + str(self.track2.end_time) + '.png', format='png',
                    dpi=500)

# end of file
