#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import numpy as np
import pdb
import pet


class ChirpChirp(pet.component, family="pet.instruments.inSAR.chirpChirp", implements=pet.protocols.instruments.inSAR):
    """
    Defines the Nightingale ChirpChirp instrument
    """

    wavelength = pet.properties.float()
    wavelength.default = 0.13
    wavelength.doc = "radar instrument wavelength [m]"

    look_angle = pet.properties.float()
    look_angle.default = 23.5
    look_angle.doc = "look angle of the radar instrument"

    antenna_elevation_width = pet.properties.float()
    antenna_elevation_width.default = 0.5
    antenna_elevation_width.doc = "elevation width of the antenna [m]"

    look_direction = pet.properties.str()
    look_direction.default = "right"
    look_direction.doc = "look direction of the radar instrument (right or left)"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Calculate the beamwidth
        d_inv = 1 / self.antenna_elevation_width
        self.bw = np.rad2deg(d_inv * self.wavelength)

        # Calculate the start and end look angles
        self.start_look_angle = self.look_angle - self.bw/2
        self.end_look_angle = self.look_angle + self.bw/2 

        return

    # def get_instrument_noise(self, measured_los):
    #     """
    #     Some sample noise function
    #     """
    #
    #     # Sample a gaussian with mean measured_los and std of 10 meters
    #     noise = np.random.normal(loc=measured_los, scale=10, size=len(measured_los))
    #
    #     return noise

# end of file
