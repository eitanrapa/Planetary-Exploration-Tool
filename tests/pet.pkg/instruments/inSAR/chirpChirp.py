#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet

# Make an instrument
instrument = pet.instruments.inSAR.chirpChirp(name="chirp chirp")

# Print the bandwidth
print(instrument.bw)

# end of file
