#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet

# Create a file manager
fm = pet.spicetoolkit.fileManager(folder_path="/home/user/Documents/GitHub/Planetary-Exploration-Tool/input")

# Furnish some files
fm.furnsh(names_list=["cas_enceladus_ssd_spc_1024icq_v1.bds", "pck00011_n0066.tpc",
                      "insar_6stride_26d_v7_seo.bsp", "latest_leapseconds.tls"])

# Make a planet
planet = pet.planets.enceladus(name="enceladus")

# Make an instrument
instrument = pet.instruments.chirpchirp(
    instrument_parameters_path="/home/user/Documents/GitHub/Planetary-Exploration-Tool/input/Parameters Sband.xlsx")

# Make a con ops
conops = pet.conOps.nightingale5to1(name="nightingale",
                                    body_id=-303, start_time="2046 DEC 20 15:10:40.134", planet=planet)

# Get the times defining the first five tracks
times = conops.get_five_tracks()

track = pet.operations.track(start_time=times[0], end_time=times[1], planet=planet,
                             conops=conops, instrument=instrument, temporal_resolution=20, spatial_resolution=2000)

# Calculate the positions
# track.calculate_ground_swath()

# Save the track
# track.save()

track.load()

# Define a projection
projection = pet.projections.biaxialCylindrical(
    name="biaxial cylindrical", folder_path="/home/user/Documents/GitHub/Planetary-Exploration-Tool/figs")

# Plot the track
fig, ax, globe = planet.visualize_topography(projection=projection, return_fig=True)
track.visualize_swath(projection=projection, fig=fig, globe=globe, ax=ax)

fm.clear()

# end of file
