# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

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

    body_id = "602"
    reference_id = "IAU_ENCELADUS"

    @pet.export
    def get_axes(self):
        """
        Return the axes of the planet
        :return: a, b, c, Equatorial major and minor semiaxes and polar semiaxis [m]
        """

        # Get the axes using the SPICE toolkit
        a, b, c = spice.bodvcd(bodyid=int(self.body_id), item="RADII")

        # Convert to meters
        a = a * 1e3
        b = b * 1e3
        c = c * 1e3

        # Return the axes
        return a, b, c

    @pet.export
    def get_surface_intersects(self, vectors):
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
        intersects = [spice.dskx02(handle=handle, dladsc=dla, vertex=vector * 10, raydir=vector - (vector * 10))[1] for
                      vector in vectors]

        # Convert to meters
        intersects = np.asanyarray(intersects) * 1e3

        # Return the intersects
        return intersects

    @pet.export
    def get_sub_obs_points(self, times, instrument):
        """
        Get the nadir sub-observation points of the instrument with the DSK
        :param times: Times to calculate the sub-observation point in accordance to the satellite position [s]
        :param instrument: The instrument
        :return: The x, y, z sub-observations points and the distance between the satellite and this observation points.
        """

        # Use the SPICE toolkit to calculate sub-observation points and vectors
        sub_points, surface_vectors = spice.subpnt_vector(method="nadir/dsk/unprioritized",
                                                          target=self.body_id, et=times, fixref=self.reference_id,
                                                          abcorr="None", obsrvr=str(instrument.body_id))[::2]

        # Convert to meters
        sub_points = np.asanyarray(sub_points) * 1e3
        distances = np.asanyarray(spice.vnorm_vector(v1=surface_vectors)) * 1e3

        # Return sub-observation points and distances of vectors
        return sub_points, distances

    @pet.export
    def visualize_topography(self, projection, return_fig=False):
        """
        Creates a visualization of the surface of Enceladus using the planet DSK
        :param visualization: Visualization tool to use
        :param projection: Cartopy projection
        :param return_fig: Whether to return the fig, ax, globe objects
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
        convert = pet.ext.conversions(name="conversions", a=a, b=b, c=c)

        # Convert coordinate vertices
        geodetic_coordinates = convert.geodetic(cartesian_coordinates=np.asarray([x, y, z]).T)[:, :3]

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

        # Make a light-source for topographic hillshading
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
        ax.set_title('Topography')

        # Save the plot
        plt.savefig(fname=projection.folder_path + '/topography.png', format='png', dpi=500)

# end of file
