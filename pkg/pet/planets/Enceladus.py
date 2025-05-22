# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import pet
import cspyce as spice
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from scipy.interpolate import griddata
import cartopy.crs as ccrs


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
    radar_penetration_depth.default = 15
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
    def get_surface_intersects_local_angles(self, satellite_position, raydirs):
        """
        Find the intersects of a set of vectors with the planet DSK
        :param satellite_position: Satellite position in the body-fixed reference frame [m]
        :param raydirs: The directions of the vectors
        :return: The x, y, z intersects of the vectors with the DSK [m]
        """

        # Retrieve the DSK handle
        handle = spice.kdata(which=0, kind="dsk")[3]

        # Retrieve the DSK DLA
        dla = spice.dlabfs(handle=handle)

        total_intersects = []
        total_plates = []

        # Use the SPICE toolkit to calculate the intersects
        for raydir in raydirs:
            # Check if intersect exists
            try:
                plate_ids, intersects = spice.dskx02(handle=handle, dladsc=dla,
                                                     vertex=satellite_position*1e-3, raydir=raydir)[:2]
                total_plates.append(plate_ids)
                total_intersects.append(intersects)
            except ValueError:
                # If not, do nothing
                pass

        # Convert to meters
        intersects = np.asarray(total_intersects) * 1e3

        # Calculate incidence angles
        normals = [spice.dskn02(handle, dla, plate_id) for plate_id in total_plates]
        normals = np.asarray(normals)

        # Project normal vector on incidence plane
        line_of_sight_norms = raydirs / np.linalg.norm(raydirs, axis=-1)[:, np.newaxis]
        plane_normal_vectors = np.cross(-satellite_position, line_of_sight_norms)

        # Normalize
        plane_normal_vectors = plane_normal_vectors / np.linalg.norm(plane_normal_vectors, axis=-1)[:, np.newaxis]

        # Normal projections
        normal_perps = np.sum(normals * plane_normal_vectors, axis=-1)[:, np.newaxis] * plane_normal_vectors
        normal_projs = normals - normal_perps
        normal_projs = normal_projs / np.linalg.norm(plane_normal_vectors, axis=-1)[:, np.newaxis]

        # get local incident angle
        incidence_angles = np.arccos(np.sum(-line_of_sight_norms * normal_projs, axis=1))
        incidence_angles = np.degrees(incidence_angles)

        # Get the look angles by calculating the angle between the satellite position vector and the
        # satellite position to intersect vector
        vector_to_intersects = intersects - satellite_position
        vector_to_intersects = vector_to_intersects / np.linalg.norm(vector_to_intersects, axis=-1)[:, np.newaxis]
        look_angles = np.arccos(np.sum(-satellite_position * vector_to_intersects, axis=1))

        # Return the intersects, incidence angles, and look angles
        return intersects, incidence_angles, look_angles

    @pet.export
    def get_closest_point_to_surface(self, points):
        """
        Get the closest intersect of a point with the planet and the distance to the surface
        :param points: Points to calculate the closest intersect with the DSK [m]
        :return: The x, y, z intersects of the vectors with the DSK [m] and the distance to the surface [m]
        """

        # Retrieve the DSK handle
        handle = spice.kdata(which=0, kind="dsk")[3]

        # Retrieve the DSK DLA
        dla = spice.dlabfs(handle=handle)

        # Use the SPICE toolkit to calculate the intersects
        intersects = [spice.dskx02(handle=handle, dladsc=dla, vertex=point, raydir=-1 * point)[1] for
                      point in points]

        # Convert to meters
        intersects = np.asanyarray(intersects) * 1e3

        # Get the distance
        distances = np.linalg.norm(intersects - points, axis=1)

        # Return the intersects and distances
        return intersects, distances

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
