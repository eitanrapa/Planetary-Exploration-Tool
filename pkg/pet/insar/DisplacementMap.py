#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved
import numpy as np
import pet
import h5py
from scipy.interpolate import RegularGridInterpolator
from ..ext import Conversions
import time


class DisplacementMap(pet.component):
    """

    """

    displacement_data_path = pet.properties.str()
    displacement_data_path.doc = "path to hdf5 file containing surface displacement"

    def __init__(self, name, locator, implicit, planet):
        super().__init__(name, locator, implicit)
        self.planet_axes = planet.get_axes()

    def read_displacements(self):
        """
        Returns the crustal deformation of Enceladus at a specific time in its tidal cycle
        """

        displacement_file = h5py.File(name=self.displacement_data_path, mode='r')
        east_cube = displacement_file["displacement_values/east_cube"][:]
        north_cube = displacement_file["displacement_values/north_cube"][:]
        up_cube = displacement_file["displacement_values/up_cube"][:]
        latitudes = displacement_file["latitudes"][:]
        longitudes = displacement_file["longitudes"][:]
        cycle_time = displacement_file["times"].attrs["Cycle Time"]
        times = displacement_file["times"][:]
        displacement_file.close()

        return latitudes, longitudes, times, cycle_time, east_cube, north_cube, up_cube

    def attach(self, swath):
        """

        """

        # Get list of times to observe at
        time_space = swath.time_space

        latitudes, longitudes, times, cycle_time, east_cube, north_cube, up_cube = self.read_displacements()

        # Create interpolating functions for each dimension
        east_interp = RegularGridInterpolator(points=(longitudes, latitudes, times),
                                              values=east_cube, method='linear', bounds_error=False, fill_value=None)
        north_interp = RegularGridInterpolator(points=(longitudes, latitudes, times),
                                               values=north_cube, method='linear', bounds_error=False, fill_value=None)
        up_interp = RegularGridInterpolator(points=(longitudes, latitudes, times),
                                            values=up_cube, method='linear', bounds_error=False, fill_value=None)

        for i in range(len(swath.swath_beams)):
            conversions = Conversions.Conversions(name="conversions")
            for point in swath.swath_beams[i]:
                cartesian_coordinates = point.x, point.y, point.z
                lat, long = conversions.geodetic(planet_axes=self.planet_axes,
                                                 cartesian_coordinates=cartesian_coordinates)[:2]
                modified_time = (time_space[i] % cycle_time) / cycle_time
                point.displacement = [east_interp((long, lat, modified_time)), north_interp((long, lat, modified_time)),
                                      up_interp((long, lat, modified_time))]

# end of file
