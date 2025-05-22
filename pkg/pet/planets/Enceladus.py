# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import pet
import cspyce as spice
import spiceypy as spice2
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from scipy.interpolate import griddata
import cartopy.crs as ccrs#
import pdb


class Enceladus(pet.component, family="pet.planets.enceladus", implements=pet.protocols.planet):
    """
    An encapsulation of the data and parameters used to represent the planetary body of Enceladus
    """

    body_id = pet.properties.str()
    body_id.default = "602"
    body_id.doc = "The NAIF ID of the planet Enceladus"

    reference_id = pet.properties.str()
    reference_id.default = "IAU_ENCELADUS"
    reference_id.doc = "The NAIF ID of the reference frame of Enceladus"

    tidal_cycle = pet.properties.float()
    tidal_cycle.default = 118386.8352
    tidal_cycle.doc = "The tidal cycle of Enceladus [s]"

    topography_uncertainty = pet.properties.float()
    topography_uncertainty.default = 100
    topography_uncertainty.doc = "The uncertainty in the topography of Enceladus [m]"

    surface_backscatter = pet.properties.float()
    surface_backscatter.default = 0
    surface_backscatter.doc = "The surface backscatter of Enceladus [dB]"

    radar_penetration_depth = pet.properties.float()
    radar_penetration_depth.default = 30
    radar_penetration_depth.doc = "The radar penetration depth of Enceladus [m]"

    surface_permittivity = pet.properties.float()
    surface_permittivity.default = 2.2
    surface_permittivity.doc = "Relative permittivity of the snow/ice volume [m]"

    @pet.export
    def get_axes(self):
        """
        Return the axes of the planet
        :return: a, b, c, Equatorial major and minor semi-axes and polar semi-axis
        """

        # Get the axes using the SPICE toolkit
        a, b, c = spice.bodvcd(bodyid=int(self.body_id), item="RADII")

        # Convert to meters
        a = a * 1e3
        b = b * 1e3
        c = c * 1e3

        # Return the axes
        return np.asarray([a, b, c])

    @pet.export
    def get_surface_intersects(self, vertex, raydirs):
        """
        Find the intersects of a set of vectors with the planet DSK
        :param vertex: Vectors for which to find intersects
        :param raydirs: The directions of the vectors
        :return: The x, y, z intersects of the vectors with the DSK [m]
        """

        # Retrieve the DSK handle
        handle = spice.kdata(which=0, kind="dsk")[3]

        # Retrieve the DSK DLA
        dla = spice.dlabfs(handle=handle)

        intersects = []

        # Use the SPICE toolkit to calculate the intersects
        for raydir in raydirs:
            # Check if intersect exists
            try:
                intersects.append(spice.dskx02(handle=handle, dladsc=dla, vertex=vertex*1e-3, raydir=raydir)[1])
            except ValueError:
                # If not, append NaN
                intersects.append([np.nan, np.nan, np.nan])

        # Convert to meters
        intersects = np.asanyarray(intersects) * 1e3

        # Return the intersects
        return intersects
    

    @pet.export
    def get_surface_intercept_LoS(self, pSat, LoSs):
        # Retrieve the DSK handle
        handle = spice.kdata(which=0, kind="dsk")[3]

        # Retrieve the DSK DLA
        dla = spice.dlabfs(handle=handle)

        # get intersects
        intersects = [spice.dskx02(handle=handle, dladsc=dla, vertex=pSat, raydir= LoS)[1] for
                      LoS in LoSs]
        intersects = np.asanyarray(intersects) * 1e3

        # get plate ids for normal vector and get normal vector of plate
        plate_ids = [spice.dskx02(handle=handle, dladsc=dla, vertex=pSat, raydir= LoS)[0] for
                      LoS in LoSs]
        normal = [spice.dskn02(handle, dla, plate_id) for plate_id in plate_ids]
        normal = np.asanyarray(normal)

        #project normal vector on incidence plane
        LoS_norm = LoSs / np.linalg.norm(LoSs, axis=-1)[:,np.newaxis]
        plane_normalVec = np.cross(-pSat, LoS_norm)
        plane_normalVec = plane_normalVec/np.linalg.norm(plane_normalVec, axis=-1)[:,np.newaxis]
        normal_perp = np.sum(normal*plane_normalVec, axis=-1)[:,np.newaxis] * plane_normalVec
        normal_proj = normal - normal_perp
        normal_proj = normal_proj/np.linalg.norm(plane_normalVec, axis=-1)[:,np.newaxis]

        # get local incident angle
        incAngles = np.arccos(np.sum(-LoS_norm * normal_proj, axis=1))

        # Return the intersects and incident angles
        return intersects, incAngles

    @pet.export
    def get_height_above_surface(self, points):
        """
        Find the intersects of a set of vectors with the planet DSK
        :param vectors: Vectors for which to find intersects
        :return: The x, y, z intersects of the vectors with the DSK [m]
        """

        # Retrieve the DSK handle
        handle = spice.kdata(which=0, kind="dsk")[3]

        # Retrieve the DSK DLA
        dla = spice.dlabfs(handle=handle)

        # Use the SPICE toolkit to calculate the intersects
        intersects = [spice.dskx02(handle=handle, dladsc=dla, vertex=point, raydir= -1*point)[1] for
                      point in points]

        # Convert to meters
        intersects = np.asanyarray(intersects) * 1e3

        # Return the intersects
        return np.linalg.norm(intersects-points,axis=1)
    
    @pet.export
    def get_closest_point(self, points):
        """
        Find the intersects of a set of vectors with the planet DSK
        :param vectors: Vectors for which to find intersects
        :return: The x, y, z intersects of the vectors with the DSK [m]
        """

        # Retrieve the DSK handle
        handle = spice.kdata(which=0, kind="dsk")[3]

        # Retrieve the DSK DLA
        dla = spice.dlabfs(handle=handle)

        # Use the SPICE toolkit to calculate the intersects
        intersects = [spice.dskx02(handle=handle, dladsc=dla, vertex=point, raydir= -1*point)[1] for
                      point in points]

        # Convert to meters
        intersects = np.asanyarray(intersects) * 1e3

        # Return the intersects
        return intersects

    @pet.export
    def get_sub_obs_points(self, times, campaign):
        """
        Get the nadir sub-observation points of the orbit with the DSK
        :param times: Times to calculate the sub-observation point in accordance to the satellite position [s]
        :param campaign: Concept of operation containing orbits
        :return: The x, y, z sub-observations points and the distance between the satellite and this observation points.
        """

        # Use the SPICE toolkit to calculate sub-observation points and vectors
        sub_points, surface_vectors = spice.subpnt_vector(method="nadir/dsk/unprioritized",
                                                          target=self.body_id, et=times, fixref=self.reference_id,
                                                          abcorr="None", obsrvr=str(campaign.body_id))[::2]

        # Convert to meters
        sub_points = np.asanyarray(sub_points) * 1e3
        distances = np.asanyarray(spice.vnorm_vector(v1=surface_vectors)) * 1e3

        # Return sub-observation points and distances of vectors
        return sub_points, distances

    @pet.export
    def visualize_topography(self, projection, fig=None, globe=None, ax=None, return_fig=False):
        """
        Creates a visualization of the surface of Enceladus using the planet DSK
        :param projection: Cartopy projection
        :param fig: matplotlib figure
        :param globe: cartopy globe
        :param ax: matplotlib ax
        :param return_fig: Whether to return fig, globe, ax
        :return: fig, ax, globe if return_fig is True
        """

        # Load the topography from DSK
        handle = spice.kdata(which=0, kind="dsk")[3]
        dla = spice.dlabfs(handle=handle)
        n_vert = spice.dskb02(handle=handle, dladsc=dla)[0]
        vertices = spice.dskv02(handle=handle, dladsc=dla, start=1, room=n_vert)
        x = vertices[:, 0] * 1e3  # Convert to meters
        y = vertices[:, 1] * 1e3  # Convert to meters
        z = vertices[:, 2] * 1e3  # Convert to meters

        # Get planet axes
        a, b, c = self.get_axes()

        # Create a coordinate conversion object
        coordinate_conversions = pet.conversions.coordinateConversions(name="conversions", a=a, b=b, c=c)

        # coordinate_conversions coordinate vertices
        geodetic_coordinates = coordinate_conversions.geodetic(cartesian_coordinates=np.asarray([x, y, z]).T)

        if fig is None:

            # Get the projection
            fig, ax, globe = projection.proj(planet=self)

        # Get reduced set of coordinates
        trimmed_coordinates = [(long, lat, height) for lat, long, height in
                               geodetic_coordinates if projection.south_extent < lat < projection.north_extent and
                               projection.west_extent < long < projection.east_extent]
        longitudes, latitudes, heights = zip(*trimmed_coordinates)

        # Make a grid to down-sample the topographical map for visualization
        n = 1000
        extent = [-180, 180, -90, 90]
        grid_x, grid_y = np.mgrid[extent[0]:extent[1]:n * 1j, extent[2]:extent[3]:n * 1j]

        # Interpolation grid
        grid = griddata((longitudes, latitudes), heights, (grid_x, grid_y), method='cubic', fill_value=0)

        # Make a light-source for topographic hill-shading
        cmap = plt.get_cmap('terrain')
        ls = LightSource(azdeg=315, altdeg=45)

        # Shade the grid with some arbitrary vertical exaggeration
        rgb = ls.shade(grid.T, cmap=cmap, blend_mode='soft', vert_exag=100)

        # Plot the grid and colorbar
        try:
            im = ax.imshow(rgb, extent=[extent[0], extent[1], extent[3], extent[2]],
                           transform=ccrs.PlateCarree(globe=globe), vmin=min(heights), vmax=max(heights), cmap=cmap)
        except ValueError:
            print("Exception: Try a smaller size plot")
            quit()

        # return fig, ax, globe if necessary
        if return_fig:
            return fig, ax, globe

        # Add colorbar
        plt.colorbar(im, label="Heights [m]")

        # Add labels and legend
        plt.title('Topography', pad=20)

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/enceladus_topography.png', format='png', dpi=500)

# end of file
