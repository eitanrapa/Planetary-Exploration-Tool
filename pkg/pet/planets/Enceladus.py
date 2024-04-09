# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
from pyre.units.SI import meter, second
from pyre.units.angle import degree
import pickle
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


class Enceladus(pet.component, family="pet.planets.enceladus", implements=pet.protocols.planet):
    """
    An encapsulation of the data and parameters used to represent the planetary body of Enceladus
    """

    # geometry
    major_equatorial_semiaxis = 256.6 * 1e3 * meter
    minor_equatorial_semiaxis = 251.4 * 1e3 * meter
    polar_semiaxis = 248.3 * 1e3 * meter
    orbital_time = 118386.8352 * second
    inclination = 0.009 * degree

    @pet.export
    def topography(self, coordinate_system="geodetic"):
        """
        Provides a list of control points that defines the surface of Enceladus
        """
        if (coordinate_system is not "geodetic") or (coordinate_system is not "cartesian"):
            raise Exception("Coordinate system must be either 'geodetic' or 'cartesian'")
        if coordinate_system == "geodetic":
            # Load from file
            with open('/home/erapapor/kraken-bak/Enceladus-Exploration-Twin-files/model_files/geodetic_topography.bin',
                      'rb') as file:
                topography = pickle.load(file)
        else:
            with open('/home/erapapor/kraken-bak/Enceladus-Exploration-Twin-files/model_files/cartesian_topography.bin',
                      'rb') as file:
                topography = pickle.load(file)

        return topography

    @pet.export
    def visualize_topography(self):
        """
        Creates a visualization of the surface of Enceladus
        """
        return None

    @pet.export
    def surface_deformation(self):
        """
        Returns the crustal deformation of Enceladus at a specific time in its tidal cycle
        """
        return []

    @pet.export
    def visualize_deformation(self):
        """
        Creates a visualization of the surface deformation of Enceladus at a specific time in its tidal cycle
        """
        return None

# end of file
