# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import cspyce as spice
from matplotlib import pyplot
import numpy as np
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource


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
        :return: a, b, c, Equatorial major and minor semiaxes and polar semiaxis
        """

        # Get the axes using the SPICE toolkit
        a, b, c = spice.bodvcd(bodyid=int(self.body_id), item="RADII")  # in km

        # Return the axes
        return a, b, c

    @pet.export
    def get_surface_intersects(self, vectors):
        """
        Find the intersects of a set of vectors with the planet DSK
        :param vectors: Vectors for which to find intersects
        :return: The x, y, z intersects of the vectors with the DSK
        """

        # Retrieve the DSK handle
        handle = spice.kdata(which=0, kind="dsk")[3]

        # Retrieve the DSK DLA
        dla = spice.dlabfs(handle=handle)

        # Use the SPICE toolkit to calculate the intersects
        intersects = [spice.dskx02(handle=handle, dladsc=dla, vertex=vector * 10, raydir=vector - (vector * 10))[1] for
                      vector in vectors]

        # Return the intersects
        return intersects

    @pet.export
    def get_sub_obs_points(self, times, instrument_id):
        """
        Get the nadir sub-observation points of the instrument with the DSK
        :param times: Times to calculate the sub-observation point in accordance to the satellite position
        :param instrument_id: The body ID of the instrument
        :return: The x, y, z sub-observations points and the distance between the satellite and this observation points.
        """

        # Use the SPICE toolkit to calculate sub-observation points and vectors
        sub_points, surface_vectors = spice.subpnt_vector(method="nadir/dsk/unprioritized",
                                                 target=self.body_id, et=times, fixref=self.reference_id,
                                                 abcorr="None", obsrvr=str(instrument_id))[::2]

        # Return sub-observation points and distances of vectors
        return sub_points, spice.vnorm_vector(v1=surface_vectors)

    @pet.export
    def visualize_topography(self, projection):
        """
        Creates a visualization of the surface of Enceladus using the planet DSK
        :return: None
        """

        # Get planet axes
        a, b, c = self.get_axes()

        # Define Enceladus globe
        img_globe = ccrs.Globe(semimajor_axis=a, semiminor_axis=c, ellipse=None)

        # Create a circular map using Cylindrical projection
        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.SouthPolarStereo(central_longitude, globe=img_globe)})

        # Plot the circular map
        ax.set_global()

        # Zoom in on South Pole
        ax.set_extent([-180, 180, -90, -30], crs=ccrs.PlateCarree(globe=img_globe))

        # Load the topography from DSK
        handle = spice.kdata(which=0, kind="dsk")[3]
        dla = spice.dlabfs(handle=handle)
        n_vert = spice.dskb02(handle=handle, dladsc=dla)[0]
        vertices = spice.dskv02(handle=handle, dladsc=dla, start=1, room=n_vert)
        x = vertices[:, 0]
        y = vertices[:, 1]
        z = vertices[:, 2]

        # Get a reduced set of topographical points for efficiency
        trimmed_coordinates = [(long, lat, height) for long, lat, height in geo_topo if lat < -30]
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
        im = ax.imshow(rgb, extent=[extent[0], extent[1], extent[3], extent[2]],
                       transform=ccrs.PlateCarree(globe=img_globe), vmin=min(heights), vmax=max(heights), cmap=cmap)

        # Add colorbar
        plt.colorbar(im, label="Heights [m]")

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

        # Add labels and legend
        ax.set_title('Topography')

        # Show the plot
        plt.savefig('/home/erapapor/kraken-bak/Enceladus-Exploration-Twin-files/coverage_maps/Topo_map',
                    format='png', dpi=500)

    def visualize_deformation(self, displacements):
        """
        Creates a visualization of the surface deformation of Enceladus at a specific time in its tidal cycle
        """


        return

# end of file
