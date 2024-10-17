#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet

# Create a file manager
fm = pet.spicetoolkit.fileManager(folder_path="/home/user/Documents/other")

# Furnish some files
fm.furnsh(names_list=["cas_enceladus_ssd_spc_1024icq_v1.bds", "pck00011_n0066.tpc",
                      "insar_6stride_26d_v7_seo.bsp", "latest_leapseconds.tls"])

# Make a planet
planet = pet.planets.enceladus(name="enceladus")

# Make an instrument
instrument = pet.instruments.nightingale(name="nightingale", body_id=-303, start_look_angle=20, end_look_angle=30,
                                         start_time="2046 DEC 20 15:10:40.134", wavelength=0.13, planet=planet)

# Get the times defining the first five tracks
times = instrument.get_five_tracks()

# Get the orbit cycle time of the instrument
orbit_cycle_time = instrument.orbit_cycle

track = pet.mission.track(start_time=times[0], end_time=times[1], planet=planet, instrument=instrument,
                          temporal_resolution=20, spatial_resolution=2000)
track.calculate_ground_swath()

# Define a projection
projection = pet.projections.biaxialPlanar(name="biaxial planar", central_latitude=-90, north_extent=-30,
                                           folder_path="/home/user/Documents/GitHub/Planetary-Exploration-Tool/figs")

# Plot interferogram
track.visualize_swath(projection=projection)

fm.clear()

# end of file
