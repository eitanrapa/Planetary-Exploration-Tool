#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import numpy as np
import spiceypy as spice
import rasterio
from matplotlib import pyplot
from .GroundTarget import GroundTarget

class GroundSwath(pet.component):
    """

    """

    start_time = pet.properties.str()
    start_time.doc = "start time of swath"

    end_time = pet.properties.str()
    end_time.doc = "end time of swath"

    time_interval = pet.properties.float()
    time_interval.doc = "time between observations (s)"

    ground_resolution = pet.properties.float()
    ground_resolution.doc = "spatial distance between observations (km)"

    def __init__(self, name, locator, implicit):
        super().__init__(name, locator, implicit)
        self.time_space = None
        self.swath_beams = []
        self.raster_fp = ""

    def make_raster(self):
        z = np.empty((len(self.swath_beams), len(self.swath_beams[0])))
        for i in range(len(self.swath_beams)):
            for j in range(len(self.swath_beams[0])):
                z[i][j] = self.swath_beams[i][j].z

        self.raster_fp = '/tmp/new.tif'
        new_dataset = rasterio.open(fp='/tmp/new.tif', mode='w', driver='GTiff',
                                    height=z.shape[0], width=z.shape[1],
                                    count=1, dtype=z.dtype, crs='+proj=cart')

        new_dataset.write(z, 1)
        new_dataset.close()

    def calculate_swath(self, instrument, planet):
        """

        """

        # Create list of times to observe at
        self.time_space = np.arange(instrument.convert_time(time=self.start_time),
                                    instrument.convert_time(time=self.end_time) +
                                    self.time_interval, self.time_interval)

        # Get planet axes
        a, b, c = planet.get_axes()

        # Get the number of rays to intercept with DSK in accordance with spatial resolution
        num_theta = int((2 * np.pi / self.ground_resolution) * np.average((a, b, c)))

        # Iterate through the times
        for time in self.time_space:

            # Create list of points on azimuthal beam
            azimuthal_points = []

            # Get the state of the satellite
            sat_position, sat_velocity = instrument.get_state(target_body_id=planet.body_id, time=time,
                                                              reference_body=planet.reference_id)

            # Create a SPICE plane with the normal as the satellite velocity originating at planet center
            plane = spice.nvp2pl(normal=sat_velocity, point=[0, 0, 0])

            # Calculate SPICE ellipse that intercept the SPICE plane with planet triaxial ellipsoid.
            ellipse = spice.inedpl(a=a, b=b, c=c, plane=plane)

            # Get center of SPICE ellipse and spanning vectors
            center, smajor, sminor = spice.el2cgv(ellipse=ellipse)

            # Get the angles to iterate through for the spanning vectors
            thetas = np.linspace(np.pi, -np.pi, num=num_theta, endpoint=False)[::-1]

            # Iterate through the angles
            for theta in thetas:

                # Calculate the vector from the origin
                vector = np.cos(theta) * smajor + np.sin(theta) * sminor

                # Get the vector intersect with the planet DSK
                intersect = planet.get_surface_intersect(vector=vector)

                # Make a GroundTarget object with the x, y, z coordinates of intersect
                ground = GroundTarget(name="{}".format(theta),
                                                   x=intersect[0], y=intersect[1], z=intersect[2])

                # Calculate and check position of groundTarget with respect to satellite
                if ground.point_position_relative_to_satellite(instrument=instrument, time=time,
                                                               planet=planet) == "right":

                    # Calculate and check look angle of groundTarget with respect to satellite
                    if (instrument.start_look_angle <
                            ground.get_angles_cartesian(instrument=instrument, time=time, planet=planet) <
                            instrument.end_look_angle):

                        # Append the groundTarget to azimuthal beam
                        azimuthal_points.append(ground)

                        # Get vector from groundTarget to satellite
                        new_vector = ground.get_position() - sat_position

                        # For all points in azimuthal beam, check there are no duplicates
                        for i in range(len(azimuthal_points) - 1):

                            # Get iteration position
                            iteration_position = azimuthal_points[i].get_position()

                            # Get iteration vector
                            iteration_vector = iteration_position - sat_position

                            # Check if the vectors are on the same path
                            if np.abs(np.dot(new_vector / np.linalg.norm(new_vector),
                                             iteration_vector / np.linalg.norm(iteration_vector))) == 1:

                                # Check if the iteration_vector is closer than the new vector
                                if (np.linalg.norm(iteration_position - sat_position) <
                                        np.linalg.norm(ground.get_position() - sat_position)):
                                    # Remove newest groundTarget
                                    azimuthal_points.pop()

                                    # Stop loop
                                    break

            # Add all accepted points in azimuthal beam to swath
            self.swath_beams.append(azimuthal_points)

    def visualize(self):
        """

        """

        self.make_raster()
        src = rasterio.open(fp=self.raster_fp)
        pyplot.imshow(src.read(1), cmap='pink')
        pyplot.show()

# end of file
