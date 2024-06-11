# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import cspyce as spice
import numpy as np
import matplotlib.pyplot as plt
import os
from ..ext import conversions
from ..projections import plottingTools


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
    def get_sub_obs_points(self, times, instrument):
        """
        Get the nadir sub-observation points of the instrument with the DSK
        :param times: Times to calculate the sub-observation point in accordance to the satellite position
        :param instrument: The instrument
        :return: The x, y, z sub-observations points and the distance between the satellite and this observation points.
        """

        # Use the SPICE toolkit to calculate sub-observation points and vectors
        sub_points, surface_vectors = spice.subpnt_vector(method="nadir/dsk/unprioritized",
                                                          target=self.body_id, et=times, fixref=self.reference_id,
                                                          abcorr="None", obsrvr=str(instrument.body_id))[::2]

        # Return sub-observation points and distances of vectors
        return sub_points, spice.vnorm_vector(v1=surface_vectors)

    @pet.export
    def visualize_topography(self, projection, north_extent=90, south_extent=-90, east_extent=180, west_extent=-180,
                             return_fig=False):
        """
        Creates a visualization of the surface of Enceladus using the planet DSK
        :return: None
        """

        # Load the topography from DSK
        handle = spice.kdata(which=0, kind="dsk")[3]
        dla = spice.dlabfs(handle=handle)
        n_vert = spice.dskb02(handle=handle, dladsc=dla)[0]
        vertices = spice.dskv02(handle=handle, dladsc=dla, start=1, room=n_vert)
        x = vertices[:, 0]
        y = vertices[:, 1]
        z = vertices[:, 2]

        geodetic_coordinates = np.asarray(
            conversions.geodetic(planet=self, cartesian_coordinates=np.asarray([x, y, z]).T))

        fig, ax, globe = projection.proj(planet=self)

        fig, ax = plottingTools.grid_plot(fig=fig, ax=ax, globe=globe, geodetic_coordinates=geodetic_coordinates[:, :3],
                                          north_extent=north_extent, south_extent=south_extent,
                                          east_extent=east_extent, west_extent=west_extent)

        if return_fig:
            return fig, ax, globe

        # Add labels and legend
        ax.set_title('Topography')

        # Get the current working directory
        current_dir = os.getcwd()

        # Construct the path three directories up
        path_three_dirs_up = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir, os.pardir))

        # Define the target directory within the three-up directory
        target_dir = os.path.join(path_three_dirs_up, 'figs')

        # Define the full path to save the plot
        full_path = os.path.join(target_dir, 'topography.png')

        # Show the plot
        plt.savefig(fname=full_path, format='png', dpi=500)

# end of file
