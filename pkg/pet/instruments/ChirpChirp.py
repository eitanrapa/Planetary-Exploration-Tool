#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import numpy as np
import pet
import pandas as pd
import pathlib


class ChirpChirp(pet.component, family="pet.instruments.chirpchirp", implements=pet.protocols.instrument):
    """
    Defines the Nightingale ChirpChirp instrument
    """

    def __init__(self, instrument_parameters_path, **kwargs):
        self.instrument_parameters_path = pathlib.PosixPath(instrument_parameters_path)
        super().__init__(**kwargs)

        # Parameter table
        params = pd.read_excel(self.instrument_parameters_path, skiprows=1, usecols=[0, 1, 2, 3, 4], index_col=0)
        self.wavelength = params.loc["Radar Wavelength", "Value"]
        look_angle = params.loc["Antenna Off-nadir Boresight Angle", "Value"]
        antenna_elevation_width = params.loc["Antenna Elevation Width", "Value"]
        d_inv = 1 / antenna_elevation_width
        bw = d_inv*self.wavelength
        self.bw = np.degrees(bw)
        self.start_look_angle = look_angle - self.bw
        self.end_look_angle = look_angle + self.bw

# end of file
