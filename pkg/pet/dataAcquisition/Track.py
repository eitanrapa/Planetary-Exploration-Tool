#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import sys
import cartopy.crs as ccrs
import cspyce as spice
import matplotlib.pyplot as plt
import numpy as np
import pet
import xarray as xr
from tqdm import tqdm
from scipy.interpolate import griddata
import pdb

def create_LoS(position, velocity, look_ang, force_zero_Doppler=True, right_looking=True):
    """ Find the LoS vector corresponding to a certain position and velocity
        vector, taking into consideration the look angle and squint

        :param position: antenna origins (satellite's position) [m].
        :type position: float 2-D array [N x 3]
        :param velocity: velocity of the antenna [m/s].
        :type velocity: float 2-D array [N x 3]
        :param look_ang: look angle [radians]
        :param squint_a: antenna squint, defaults to zero [deg]

        returns: float 3-D array
        """
    Nla = look_ang.size
    if velocity.ndim == 1:
        ax = 0
        Nv = 1
    else:
        ax = 1
        Nv = velocity.shape[0]
    # Calculate velocity and positions versor
    v_ver = velocity/np.linalg.norm(velocity, axis=ax).reshape((Nv, 1))
    r_ver = position/np.linalg.norm(position, axis=ax).reshape((Nv, 1))
    n_ver = np.cross(v_ver, r_ver)    # cross product of versors
    if force_zero_Doppler:
        r_ver2 = np.cross(n_ver, v_ver)
    else:
        r_ver2 = r_ver
    if not right_looking:
        look_ang = -look_ang
    LoS = (-np.cos(look_ang.reshape(1, Nla, 1))*r_ver2.reshape(Nv, 1, 3) +
           np.sin(look_ang.reshape(1, Nla, 1))*n_ver.reshape(Nv, 1, 3))
    
    return LoS

class Track(pet.component, family="pet.dataAcquisition.track", implements=pet.protocols.dataAcquisition):
    """
    An object containing the swath of a satellite pass. Contains the start, end, and temporal_resolution of the
    satellite pass
    as well as given ground resolution to calculate vectors with. Also contains the satellite time vector and a
    2-D array of the azimuthal beams by row and GroundTargets populating the rows.
    """

    start_time = pet.properties.float()
    start_time.doc = "start time of track"

    end_time = pet.properties.float()
    end_time.doc = "end time of track"

    planet = pet.protocols.planet()
    planet.doc = "target planet"

    campaign = pet.protocols.campaigns.orbiter()
    campaign.doc = "campaign (must be an orbiter)"

    instrument = pet.protocols.instruments.inSAR()
    instrument.doc = "observation instrument"

    temporal_resolution = pet.properties.float()
    temporal_resolution.default = 10
    temporal_resolution.doc = "temporal resolution of orbit [s]"

    spatial_resolution = pet.properties.float()
    spatial_resolution.default = 200
    spatial_resolution.doc = "spatial resolution of ground track [m]"

    interpol2FinalRes = pet.properties.int()
    interpol2FinalRes.default = 0
    interpol2FinalRes.doc = "flag for interpolation to final grid"
    
    interferogram_resolution_az = pet.properties.float()
    interferogram_resolution_az.default = 300
    interferogram_resolution_az.doc = "spatial resolution of final interferogram in azimuth [m]"

    interferogram_resolution_rg = pet.properties.float()
    interferogram_resolution_rg.default = 300
    interferogram_resolution_rg.doc = "spatial resolution of final interferogram in range [m]"
    

    data = None

    @classmethod
    def from_file(cls, planet, campaign, instrument, file_name):
        """
        Load a track from an HDF5 file
        :param planet: Planet object
        :param campaign: Campaign object
        :param instrument: Instrument object
        :param file_name: File name of HDF5 file
        :return: Track object
        """

        # Open the HDF5 file in read mode
        data = xr.open_dataarray(filename_or_obj=file_name)

        # Get the start and end times
        start_time = data.attrs["start_time"]
        end_time = data.attrs["end_time"]

        # Create the object
        obj = cls(name="track" + str(np.random.rand()), start_time=start_time, end_time=end_time,
                  planet=planet, campaign=campaign, instrument=instrument,
                  temporal_resolution=data.attrs["temporal_resolution"],
                  spatial_resolution=data.attrs["spatial_resolution"])

        obj.data = data  # Restore computed result

        return obj

    @classmethod
    def from_data_array(cls, planet, campaign, instrument, data):
        """
        Create a track object from a data array
        :param planet: Planet object
        :param campaign: Campaign object
        :param instrument: Instrument object
        :param data: Data array
        :return: Track object
        """

        # Get the start and end times
        start_time = data.attrs["start_time"]
        end_time = data.attrs["end_time"]

        # Create the object
        obj = cls(name="track" + str(start_time) + str(end_time),
                  start_time=start_time, end_time=end_time,
                  planet=planet, campaign=campaign, instrument=instrument,
                  temporal_resolution=data.attrs["temporal_resolution"],
                  spatial_resolution=data.attrs["spatial_resolution"])

        obj.data = data  # Restore computed result

        return obj

    @classmethod
    def from_files(cls, planet, campaign, instrument, file_list):
        """
        Load a list of tracks from HDF5 files
        :param planet: Planet object
        :param campaign: Campaign object
        :param instrument: Instrument object
        :param file_list: List of file names
        :return: List of track objects
        """

        # Load all the files
        return [cls.from_file(planet=planet, campaign=campaign, instrument=instrument, file_name=file) for
                file in file_list]

    def save(self, file_name):
        """
        Save the track to an HDF5 file
        :return: Nothing returned
        """

        # Open HDF5 file
        self.data.to_netcdf(file_name, engine="netcdf4")

    def create_data_array(self, cartesian_coordinates, geodetic_coordinates, inc_angles, times):
        """
        Create a xarray with the input data
        :param cartesian_coordinates: x, y, z coordinates to be saved
        :param geodetic_coordinates: lat, long, height coordinates to be saved
        :param look_angles: Calculated look angles from satellite
        :param times: Times of observation for points
        :return: Nothing returned
        """

        # Create the xarray Dataset
        da = xr.DataArray(
            data=inc_angles,
            dims=["points"],
            coords={
                "sat_pos_time": ("points", times),
                "time": ("points", times),
                "x": ("points", np.asarray([point[0] for point in cartesian_coordinates])),
                "y": ("points", np.asarray([point[1] for point in cartesian_coordinates])),
                "z": ("points", np.asarray([point[2] for point in cartesian_coordinates])),
                "latitude": ("points", np.asarray([point[0] for point in geodetic_coordinates])),
                "longitude": ("points", np.asarray([point[1] for point in geodetic_coordinates])),
                "height": ("points", np.asarray([point[2] for point in geodetic_coordinates]))},
            name="inc_angles",
            attrs=dict(
                body_id=self.campaign.body_id,
                start_time=self.start_time,
                end_time=self.end_time,
                temporal_resolution=self.temporal_resolution,
                spatial_resolution=self.spatial_resolution
            ),
        )

        # Save xarray to object
        self.data = da

    def modify_time(self, time):
        """
        Modify the time of the track
        :param time: Time to add to the track
        :return: Nothing returned
        """

        self.start_time = self.start_time + time
        self.end_time = self.end_time + time
        self.data["time"].values = self.data["time"].values + time
        self.data.attrs["start_time"] = self.data.attrs["start_time"] + time
        self.data.attrs["end_time"] = self.data.attrs["end_time"] + time

    def convert_positions(self, cartesian_positions):
        """
        Convert flattened points to geodetic coordinates
        :param cartesian_positions: Points to convert
        :return: Points in geodetic coordinates
        """

        # Get planet axes
        a, b, c = self.planet.get_axes()

        # Create a coordinate conversion object
        coordinate_conversions = pet.conversions.coordinateConversions(name="conversions", a=a, b=b, c=c)

        # Get the corresponding latitude and longitude with a triaxial ellipsoid conversion
        flattened_swath_coordinates = np.asarray(coordinate_conversions.geodetic(
            cartesian_coordinates=cartesian_positions))

        return flattened_swath_coordinates


    def calculate_ground_swath_new(self):
        """
        Vectorized calculation of the ground swath.
        :return: Nothing returned
        """
        print("Starting Swath Computation...", file=sys.stderr)
        # Get the time space
        time_space = np.arange(self.start_time, self.end_time + self.temporal_resolution,
                            self.temporal_resolution)

        # Get planet axes
        a, b, c = self.planet.get_axes()

        print("     Getting satellite positions...", file=sys.stderr)
        # Get positions and velocities of the satellite for the observation times
        satellite_positions, satellite_velocities = self.campaign.get_states(times=time_space)

        print("     Getting antenna beam sampling...", file=sys.stderr)
        # Get the number of rays to intercept with DSK in accordance with spatial resolution (rough approximation)
        angSpan = np.deg2rad(self.instrument.end_look_angle - self.instrument.start_look_angle)
        meanHeight = np.mean(np.linalg.norm(satellite_positions,axis=1) - np.average((a, b, c)))
        meanRange = meanHeight/np.cos(np.mean(angSpan))
        num_theta = int((angSpan / self.spatial_resolution) * meanRange)
        # Get look angle array
        thetas = np.linspace(self.instrument.start_look_angle, self.instrument.end_look_angle,num_theta)
        thetas = np.deg2rad(thetas)

        print("     Getting lines of sight...", file=sys.stderr)
        LoS = create_LoS(satellite_positions, satellite_velocities,thetas, right_looking=True)

        print("     Getting intercept points and local incident angles...", file=sys.stderr)
        icp = np.zeros(LoS.shape)
        incAng = np.zeros((LoS.shape[0], LoS.shape[1]))
        for pp in tqdm(range(satellite_positions.shape[0])):
            icp[pp,:,:], incAng[pp,:] = self.planet.get_surface_intercept_LoS(pSat=satellite_positions[pp,:]*1e-3, LoSs=LoS[pp,:,:])
        
        # if wanted, do interpolation to interferogram grid
        if self.interpol2FinalRes:
            print("     Interpolation to final interferogram grid...", file=sys.stderr)
            ranges = np.linalg.norm(icp - satellite_positions[:,np.newaxis,:], axis=-1)
            # Get resolution of final product
            res_az = self.interferogram_resolution_az
            res_rg = self.interferogram_resolution_rg
            v_ground = np.mean(np.linalg.norm(satellite_velocities,axis=-1)) * np.average((a, b, c)) / \
                            (np.average((a, b, c)) + np.mean(np.linalg.norm(satellite_positions,axis=-1)) - np.average((a, b, c))) # ground velocity (is an approximation... can be improved later)
            
            # get interpolation grid
            duration = time_space[-1] - time_space[0]
            aimuth_vec = np.arange(0, duration*v_ground + res_az, res_az)
            range_vec = np.arange(np.min(ranges), np.max(ranges) + res_rg, res_rg)
            grid_azimuth, grid_range = np.meshgrid(aimuth_vec, range_vec, indexing='ij')
            # get data and coordinates
            cartesian_positions_aux = np.reshape(icp, (icp.shape[0]*icp.shape[1], 3))
            range_coord = np.reshape(ranges, ranges.size)
            time_coord = np.reshape(np.tile(time_space[:, np.newaxis], (1, icp.shape[1])), icp.shape[0]*icp.shape[1])
            azimuth_coord = time_coord * v_ground
            azimuth_coord -= azimuth_coord[0]
            # do interpolation
            gridded_x =  griddata(points=(azimuth_coord, range_coord), values=cartesian_positions_aux[:,0],
                                  xi=(grid_azimuth, grid_range), method='cubic')
            gridded_y =  griddata(points=(azimuth_coord, range_coord), values=cartesian_positions_aux[:,1],
                                  xi=(grid_azimuth, grid_range), method='cubic')
            gridded_z =  griddata(points=(azimuth_coord, range_coord), values=cartesian_positions_aux[:,2],
                        xi=(grid_azimuth, grid_range), method='cubic')
            gridded_inc =  griddata(points=(azimuth_coord, range_coord), values=np.reshape(incAng, incAng.size),
                        xi=(grid_azimuth, grid_range), method='cubic')
            mask = np.zeros(gridded_x.shape, dtype=int)
            for aa in range(len(aimuth_vec)):
                # pdb.set_trace()
                idxClose = np.argmin(abs((time_space-time_space[0])*v_ground - aimuth_vec[aa]))
                mask_aux = np.where((range_vec<ranges[idxClose,0]) | (range_vec>ranges[idxClose,-1]))[0]
                mask[aa,mask_aux] = 1
                # pdb.set_trace()
            gridded_x[mask==1] = np.nan
            gridded_y[mask==1] = np.nan
            gridded_z[mask==1] = np.nan
            gridded_inc[mask==1] = np.nan
            # pdb.set_trace()
            
            # reformat
            cartesian_positions = np.stack((gridded_x, gridded_y, gridded_z), axis=-1)
            cartesian_positions = np.reshape(cartesian_positions, (cartesian_positions.shape[0]*cartesian_positions.shape[1], 3))
            flat_incAng = np.reshape(gridded_inc, gridded_inc.size)
            flat_times = np.reshape(grid_azimuth, grid_azimuth.size)/v_ground + time_space[0]
            #remove nan
            validIdx = np.where(~np.isnan(flat_incAng))[0]
            cartesian_positions = cartesian_positions[validIdx,:]
            cartesian_positions = self.planet.get_closest_point(cartesian_positions)
            flat_incAng = flat_incAng[validIdx]
            flat_times = flat_times[validIdx]

        else:
            print("     No interpolation to final interferogram grid...", file=sys.stderr)
            cartesian_positions = np.reshape(icp, (icp.shape[0]*icp.shape[1], 3))
            flat_incAng = np.reshape(incAng, incAng.size)
            flat_times = np.reshape(np.tile(time_space[:, np.newaxis], (1, icp.shape[1])), icp.shape[0]*icp.shape[1])
        
        print("     Coordinate conversion...", file=sys.stderr)
        swath_coordinates = self.convert_positions(cartesian_positions=cartesian_positions)        
        
        
        self.create_data_array(cartesian_coordinates=cartesian_positions, geodetic_coordinates=swath_coordinates,
                            inc_angles=flat_incAng, times=flat_times)
        print("Swath calculation finished!", file=sys.stderr)

    def calculate_ground_swath(self):
        """
        Vectorized calculation of the ground swath.
        :return: Nothing returned
        """

        # Get the time space
        time_space = np.arange(self.start_time, self.end_time + self.temporal_resolution,
                               self.temporal_resolution)

        # Get planet axes
        a, b, c = self.planet.get_axes()

        # Get the number of rays to intercept with DSK in accordance with spatial resolution
        num_theta = int((2 * np.pi / self.spatial_resolution) * np.average((a, b, c)))

        print("Getting satellite positions...", file=sys.stderr)

        # Get positions and velocities of the satellite for the observation times
        satellite_positions, satellite_velocities = self.campaign.get_states(times=time_space)

        # Normalize the velocity and position vectors
        vhats = np.asarray([satellite_velocity / np.linalg.norm(satellite_velocity)
                            for satellite_velocity in satellite_velocities])
        rhats = np.asarray([sat_position / np.linalg.norm(sat_position)
                            for sat_position in satellite_positions])

        # Get the azimuthal unit vectors
        # Use radial direction instead of projection (more stable)
        u1s = [-rhat for rhat in rhats]  # radially inward

        if self.instrument.look_direction == "right":
            u2s = [-np.cross(vhat, u1) for vhat, u1 in zip(vhats, u1s)]
        elif self.instrument.look_direction == "left":
            u2s = [np.cross(vhat, u1) for vhat, u1 in zip(vhats, u1s)]
        else:
            raise ValueError("Invalid look direction. Must be 'left' or 'right'.")
        u2s = [u2 / np.linalg.norm(u2) for u2 in u2s]

        swath_beams = []

        # Iterate through each azimuthal beam
        for i in tqdm(range(len(time_space)), desc="Calculating planet surface intersects..."):

            ray_dirs = [np.cos(np.deg2rad(theta)) * u1s[i] + np.sin(np.deg2rad(theta)) * u2s[i]
                        for theta in np.linspace(self.instrument.start_look_angle,
                                                 self.instrument.end_look_angle, num_theta)]

            # Get the intersects with the planet DSK
            intersects = self.planet.get_surface_intersects(vertex=satellite_positions[i], raydirs=ray_dirs)

            # Make groundTarget objects with the intersects
            swath_beams.append([[intersect[0], intersect[1], intersect[2]] for intersect in intersects
                                if intersect[0] is not np.nan])

        # Calculate the look angle of each GroundTarget for each beam
        look_angles = self.get_angles_cartesian(swath_beams=swath_beams,
                                                times=time_space, satellite_positions=satellite_positions)

        cartesian_positions = np.asanyarray([coord for sublist in swath_beams for coord in sublist])
        flat_angles = np.asarray([angle for sublist in look_angles for angle in sublist])
        times = [[time] * len(beams) for time, beams in zip(time_space, swath_beams)]
        flat_times = np.asarray([time for sublist in times for time in sublist])
        swath_coordinates = self.convert_positions(cartesian_positions=cartesian_positions)

        self.create_data_array(cartesian_coordinates=cartesian_positions, geodetic_coordinates=swath_coordinates,
                               look_angles=flat_angles, times=flat_times)


    def get_angles_cartesian(self, swath_beams, times, satellite_positions):
        """
        Get the look angles between the GroundTargets and the satellite in degrees.
        :param swath_beams: Azimuthal beams of the swath
        :param times: Times of observation
        :param satellite_positions: Array of x, y, z, positions of the satellite
        :return: Absolute look angle between the surface points and the corresponding satellite position in degrees
        """

        # Calculate the satellite intersects and heights with the shape
        intersects, satellite_heights = self.planet.get_sub_obs_points(times=times, campaign=self.campaign)

        # Calculate the distance between center and satellite intersects
        satellite_radii = np.asarray([np.linalg.norm(intersect - [0, 0, 0]) for intersect in intersects])

        # Iterate through the beams
        look_angles_per_beam = []
        for i in tqdm(range(len(swath_beams)), desc="Calculating look angles..."):

            # Get surface positions per beam
            surface_positions = np.asarray([(point[0], point[1], point[2]) for point in swath_beams[i]])

            # Get the distance between the satellite and the surface points
            distances = np.asarray([np.linalg.norm(surface_position - satellite_positions[i])
                                    for surface_position in surface_positions])

            # Calculate cosines of the angles
            z_plus_re = satellite_heights[i] + satellite_radii[i]
            altitude = np.asarray([np.linalg.norm(surface_position - [0, 0, 0])
                                   for surface_position in surface_positions])
            cosine_look_angle = (distances ** 2 + z_plus_re ** 2 - altitude ** 2) / (2 * distances * z_plus_re)

            # Use arc cosine to get the angles in radians
            look_angles_radians = np.arccos(cosine_look_angle)

            # Convert radians to degrees
            look_angles_degrees = np.degrees(look_angles_radians)

            # Yield output
            look_angles_per_beam.append(look_angles_degrees)

        # Return the look angles 2-D array
        return look_angles_per_beam

    def visualize_swath(self, projection, fig=None, globe=None, ax=None, return_fig=False):
        """
        Visualize the ground swath
        :param projection: Cartopy projection
        :param fig: matplotlib figure
        :param globe: cartopy globe
        :param ax: matplotlib ax
        :param return_fig: Whether to return fig, globe, ax
        :return: fig, ax, globe if return_fig is True
        """

        if fig is None:

            # Get the projection
            fig, ax, globe = projection.proj(planet=self.planet)

        # Access longitude and latitude coordinates
        longitudes = self.data["longitude"].values
        latitudes = self.data["latitude"].values

        # Plot points on the map
        ax.scatter(longitudes, latitudes, transform=ccrs.PlateCarree(globe=globe),
                   color='red', marker='o', s=0.1, alpha=0.25)

        # return fig, ax, globe if necessary
        if return_fig:

            return fig, ax, globe

        # Add labels and legend
        plt.title('Swath', pad=20)

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/swath_' + str(self.campaign.body_id) + '_' +
                    str(self.start_time) + '_' + str(self.end_time) + '.png', format='png', dpi=500)

# end of file
