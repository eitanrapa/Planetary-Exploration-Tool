#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import h5py
from scipy.interpolate import RegularGridInterpolator


class DisplacementMap(pet.component):
    """

    """

    displacement_data_path = pet.properties.str()
    displacement_data_path.doc = "path to hdf5 file containing surface displacement"

    def read_displacements(self):
        """
        Returns the crustal deformation of Enceladus at a specific time in its tidal cycle
        """

        displacement_file = h5py.File(name=self.displacement_data_path, mode='r')
        east_cube = displacement_file["displacement_values/east_cube"]
        north_cube = displacement_file["displacement_values/north_cube"]
        up_cube = displacement_file["displacement_values/up_cube"]
        latitudes = displacement_file["latitudes"]
        longitudes = displacement_file["longitudes"]
        times = displacement_file["times"]
        cycle_time = times.attrs["Cycle Time"]

        return latitudes, longitudes, times, cycle_time, east_cube, north_cube, up_cube

    def attach(self, swath):
        """

        """

        # Get list of times to observe at
        time_space = swath.time_space

        latitudes, longitudes, times, cycle_time, east_cube, north_cube, up_cube = self.read_displacements()

        # Create interpolating functions for each dimension
        east_interp = RegularGridInterpolator(points=(latitudes, longitudes, times),
                                              values=east_cube, method='cubic')
        north_interp = RegularGridInterpolator(points=(latitudes, longitudes, times),
                                               values=north_cube, method='cubic')
        up_interp = RegularGridInterpolator(points=(latitudes, longitudes, times),
                                            values=up_cube, method='cubic')

        for i in range(len(swath.swath_beams)):
            for point in swath.swath_beams[i]:
                lat, long = geodesy.geodetic(point.x, point.y, point.z)
                modified_time = (time_space[i] % cycle_time) / cycle_time
                point.displacement = [east_interp((lat, long, modified_time)), north_interp((lat, long, modified_time)),
                                      up_interp((lat, long, modified_time))]

# end of file
