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
instrument = pet.instruments.nightingale(name="nightingale", body_id=-303, start_look_angle=20, end_look_angle=30,
                                         orbit_cycle=110020, wavelength=0.13)

# Make a planet
planet = pet.planets.enceladus(name="enceladus")

# Make a displacement map
displacements = pet.insar.displacementMap(name="base",
                                          displacement_data_path=
                                          "/home/eitanrapa/Documents/projects/other/Simulation_Het_Results.hdf5")

# Make an interferogram
interferogram = pet.insar.simpleInterferogram(name="igram", pairing_one=0, pairing_two=1, track_number=2,
                                              instrument=instrument, planet=planet, displacements=displacements,
                                              time_interval=20, ground_resolution=500, baseline=0)

# Plot interferogram
# interferogram.visualize()

fm.clear()

# end of file
