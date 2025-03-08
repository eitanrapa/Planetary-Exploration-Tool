#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import numpy as np
import pet


class ChirpChirp(pet.component, family="pet.instruments.inSAR.chirpChirp", implements=pet.protocols.instruments.inSAR):
    """
    Defines the Nightingale ChirpChirp instrument
    """

    wavelength = pet.properties.dimensional()
    wavelength.default = "13*cm"
    wavelength.doc = "radar instrument wavelength"

    look_angle = pet.properties.float()
    look_angle.default = 25
    look_angle.doc = "look angle of the radar instrument"

    antenna_elevation_width = pet.properties.dimensional()
    antenna_elevation_width.default = "0.5*m"
    antenna_elevation_width.doc = "elevation width of the antenna"

    look_direction = pet.properties.str()
    look_direction.default = "right"
    look_direction.doc = "look direction of the radar instrument"

    def __init__(self, **kwargs):
        """

        """

        super().__init__(**kwargs)

        d_inv = 1 / self.antenna_elevation_width

        self.bw = np.rad2deg(np.arcsin(d_inv * self.wavelength))

        self.start_look_angle = self.look_angle - self.bw
        self.end_look_angle = self.look_angle + self.bw

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
