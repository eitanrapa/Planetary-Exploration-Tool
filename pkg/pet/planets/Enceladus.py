# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import cspyce as spice
from matplotlib import pyplot
import numpy as np
from matplotlib.colors import LightSource
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata


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
    def visualize_topography(self):
        """
        Creates a visualization of the surface of Enceladus
        """
        handle = spice.kdata(which=0, kind="dsk")[3]
        dla = spice.dlabfs(handle=handle)
        n_vert = spice.dskb02(handle=handle, dladsc=dla)[0]
        vertices = spice.dskv02(handle=handle, dladsc=dla, start=1, room=n_vert)
        x = vertices[:, 0]
        y = vertices[:, 1]
        z = vertices[:, 2]

        # Define grid
        grid_x, grid_y = np.meshgrid(np.linspace(x.min(), x.max(), 1000),
                                     np.linspace(y.min(), y.max(), 1000))

        # Interpolate z values on the grid
        grid_z = griddata((x, y), z, (grid_x, grid_y), method='cubic')

        fig = pyplot.figure(figsize=(10, 7))
        ax = pyplot.axes(projection='3d')
        cmap = pyplot.get_cmap('terrain')
        surf = ax.plot_surface(grid_x, grid_y, grid_z, cmap=cmap, edgecolor='none')
        # ls = LightSource(azdeg=315, altdeg=45)
        # rgb = ls.shade(grid.T, cmap=cmap, blend_mode='soft', vert_exag=100)
        fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
        ax.set_title('Surface plot')
        pyplot.show()

        return

    def visualize_deformation(self, displacements):
        """
        Creates a visualization of the surface deformation of Enceladus at a specific time in its tidal cycle
        """


        return

# end of file
