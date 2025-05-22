#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import pet

# Make an instrument
instrument = pet.instruments.inSAR.chirpChirp(name="chirp chirp")

# Make a planet
planet = pet.planets.enceladus(name="enceladus")

# Print the bandwidth
sigma_phase, corr_tot, n_looks, nesn, sigma0 =\
    instrument.get_instrument_noise(planet, baseline=10, satellite_velocities=[100, 120, 140],
                                    look_angles=[20, 30, 40], incidence_angles=[10, 20 ,30],
                                    distances=[100e3, 100e3, 100e3],
                                    variable_backscatter=True)

print("Sigma phase: ", sigma_phase)
print("Correlation total: ", corr_tot)
print("Number of looks: ", n_looks)
print("NESN: ", nesn)
print("Sigma0: ", sigma0)

# end of file
