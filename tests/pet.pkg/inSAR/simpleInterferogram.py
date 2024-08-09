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

# Make a planet
planet = pet.planets.enceladus(name="enceladus")

# Make an instrument
instrument = pet.instruments.nightingale(name="nightingale", body_id=-303, start_look_angle=20, end_look_angle=30,
                                         start_time="2046 DEC 20 15:10:40.134", wavelength=0.13, planet=planet)

# # Make a displacement map
# displacements = pet.insar.displacementMap(name="het",
#                                           displacement_data_path=
#                                           "/home/eitanrapa/Documents/projects/other/Simulation_Het_Results.hdf5",
#                                           planet=planet)

# Make a displacement map
displacements = pet.insar.displacementMap(name="Het",
                                          displacement_data_path=
                                          "/home/eitanrapa/Documents/projects/other/Simulation_Het_Results.hdf5",
                                          planet=planet)

# # Make an interferogram
# interferogram = pet.insar.simpleInterferogram(name="igram", pairing_one=0, pairing_two=1, track_number=2,
#                                               instrument=instrument, planet=planet, displacements=displacements,
#                                               time_interval=10, ground_resolution=200, baseline=0)
#
# # Save interferogram
# interferogram.save(path="/home/eitanrapa/Documents/projects/pet/files")

# Load interferogram
interferogram = pet.insar.simpleInterferogram(name="igram", instrument=instrument, planet=planet,
                                              displacements=displacements,
                                              load_name="/home/eitanrapa/Documents/GitHub/"
                                                        "Planetary-Exploration-Tool/files/interferogram.hdf5")

interferogram.pairing_one = 0
interferogram.pairing_two = 1

interferogram.recalculate_igram(baseline=0)

# # Save interferogram
# interferogram.save(path="/home/eitanrapa/Documents/GitHub/Planetary-Exploration-Tool/files", name="het_45")

# # Make a projection
# projection = pet.projections.biaxialCylindrical(name="biaxial cylindrical",
#                                                 folder_path="/home/eitanrapa/Documents/GitHub/"
#                                                             "Planetary-Exploration-Tool/figs")

# Define a projection
projection = pet.projections.biaxialPlanar(name="biaxial planar", central_latitude=-90, north_extent=-30,
                                           folder_path="/home/eitanrapa/Documents/"
                                                       "GitHub/Planetary-Exploration-Tool/figs")

# Plot interferogram
interferogram.visualize(projection=projection)

fm.clear()

# end of file
