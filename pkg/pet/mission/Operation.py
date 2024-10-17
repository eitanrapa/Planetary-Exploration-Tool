#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
from .Track import Track

class Operation(pet.component):
    """

    """

    def __init__(self, planet, instrument, deformation_map, acquisition_cadence, temporal_resolution=10,
                 spatial_resolution=1000, **kwargs):
        super().__init__(**kwargs)
        self.planet = planet
        self.instrument = instrument
        self.deformation_map = deformation_map
        self.acquisition_cadence = acquisition_cadence
        self.spatial_resolution = spatial_resolution
        self.temporal_resolution = temporal_resolution

    def get_interferograms(self, number_of_interferograms):
        """

        """

        for i in range(number_of_interferograms):

            # Get the times defining the first five tracks
            times = self.instrument.get_five_tracks()

            # Get the orbit cycle time of the instrument
            orbit_cycle_time = self.instrument.orbit_cycle

            # Create a groundSwath object from the track times, time interval, ground resolution, planet, and instrument
            track1 = Track(start_time, end_time, planet=self.planet, instrument=self.instrument,
                           temporal_resolution=self.temporal_resolution, spatial_resolution=self.spatial_resolution)
            track1.calculate_ground_swath()

            # Create a groundSwath object from the track times, time interval, ground resolution, planet, and instrument
            track2 = Track(start_time, end_time, planet=self.planet, instrument=self.instrument,
                           temporal_resolution=self.temporal_resolution, spatial_resolution=self.spatial_resolution)
            track2.time_space + orbit_cycle_time * self.pairing_two
            track2.swath = track1.swath
#
#
#             interferogram =
#
#                 start_time=times[self.acquisition_cadence.track_number],
#                                            end_time=times[self.track_number + 1], time_interval=time_interval,
#                                            ground_resolution=ground_resolution, planet=self.planet,
#                                            instrument=self.instrument)
#
#             # Create the second groundSwath object
#             gs2 = GroundSwath(start_time=self.groundSwath.start_time,
#                               end_time=self.groundSwath.end_time, time_interval=self.groundSwath.time_interval,
#                               ground_resolution=self.groundSwath.ground_resolution, planet=self.planet,
#                               instrument=self.instrument,
#                               copy=True)
#
#             gs2.swath_beams = copy.deepcopy(self.groundSwath.swath_beams)
#
#             # Set the correct time space for the second GroundSwath given the second pairing number
#             gs2.time_space = self.base_time_space + orbit_cycle_time * self.pairing_two
#
#             # Attach the displacements of both GroundSwaths
#             self.displacements.attach(swath1=self.groundSwath, swath2=gs2, use_mid_point=True)
#
# # end of file
