#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet

# Create a file manager
fm = pet.spicetoolkit.fileManager(folder_path="/home/eitanrapa/Documents/projects/other")

# Furnish some files
fm.furnsh(names_list=["cas_enceladus_ssd_spc_1024icq_v1.bds", "pck00011_n0066.tpc",
                      "insar_6stride_26d_v7_seo.bsp", "latest_leapseconds.tls"])

# Make an instrument
instrument = pet.instruments.nightingale(name="nightingale", body_id=-303, start_look_angle=20, end_look_angle=30)

# Make a planet
planet = pet.planets.enceladus(name="enceladus")

# Get times
times = instrument.get_five_tracks(planet)

# Make a ground swath
gs = pet.insar.groundSwath(name="1", start_time=times[1], end_time=times[2],
                           time_interval=10, ground_resolution=2000, planet=planet, instrument=instrument)

# Make a projection
projection = pet.projections.biaxialCylindrical(name="biaxial cylindrical")

# Make a visualization tool
visualization = pet.visualization.cartopyViz(name="cartopy tool",
                                             folder_path="/home/eitanrapa/Documents/projects/pet/figs")

# Plot
gs.visualize(visualization=visualization, planet=planet, instrument=instrument, projection=projection)

fm.clear()

# end of file
