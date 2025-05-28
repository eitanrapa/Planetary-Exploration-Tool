#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import pet

# Create a file manager
fm = pet.spiceTools.fileManager(folder_path="/home/user/Documents/GitHub/Planetary-Exploration-Tool/input")

# Furnish some files
fm.furnsh(names_list=["cas_enceladus_ssd_spc_1024icq_v1.bds", "pck00011_n0066.tpc",
                      "insar_6stride_26d_v7_seo.bsp", "latest_leapseconds.tls"])

# Make a planet
planet = pet.planets.enceladus(name="enceladus")

# Make an instrument
instrument = pet.instruments.inSAR.chirpChirp(name="chirp chirp")

# Make a con ops
campaign = pet.campaigns.orbiter.nightingale5to1(name="nightingale",
                                                 body_id=-303, start_time="2046 DEC 20 15:10:40.134", planet=planet)

# Get the times defining the first five tracks
times = campaign.get_five_tracks()

# Get the orbit cycle time of the instrument
orbit_cycle_time = campaign.orbit_cycle

# Make a displacement map
deformation_map = pet.natureSimulations.geophysicalModel.tidalDeformationMap(name="base",
                                                                             displacement_data_path="/home/user/"
                                                                                                    "Documents/GitHub/"
                                                                                                    "Planetary-"
                                                                                                    "Exploration-Tool/"
                                                                                                    "input/"
                                                                                                    "Simulation_Base"
                                                                                                    "_Results.hdf5",
                                                                             planet=planet)

# First track
track = pet.dataAcquisition.track.from_file(planet=planet, campaign=campaign, instrument=instrument,
                                            file_name="/home/user/Documents/GitHub/"
                                                       "Planetary-Exploration-Tool/files/track1")

interferogram = pet.dataAnalysis.surfaceDeformationAnalysis.simpleInterferogram(name="igram",
                                                                                instrument=instrument, planet=planet,
                                                                                deformation_map=deformation_map,
                                                                                track=track,
                                                                                campaign=campaign,
                                                                                time_offset_first_acquisition=0,
                                                                                time_offset_second_acquisition=
                                                                                orbit_cycle_time,
                                                                                perpendicular_baseline=10)
# Calculate interferogram
interferogram.calculate_igram()

# # Save interferogram
# interferogram.save(file_name="/home/user/Documents/GitHub/"
#                              "Planetary-Exploration-Tool/files/igram_base_1")
#
# # Load the interferogram
# interferogram = pet.dataAnalysis.surfaceDeformationAnalysis.simpleInterferogram.from_file(instrument=instrument,
#                                                                                           planet=planet,
#                                                                                           deformation_map=
#                                                                                           deformation_map,
#                                                                                           campaign=campaign,
#                                                                                           file_name="/home/user/"
#                                                                                                     "Documents/"
#                                                                                                     "GitHub/"
#                                                                                                     "Planetary"
#                                                                                                     "-Exploration-Tool/"
#                                                                                                     "files/igram_"
#                                                                                                     "base_1")

# Define a projection
projection = pet.projections.biaxialProjections.biaxialCylindrical(name="biaxial cylindrical",
                                                                   folder_path=
                                                                   "/home/user/Documents/GitHub/Planetary"
                                                                   "-Exploration-Tool/figs/")

# Plot interferogram
fig, ax, globe = planet.visualize_topography(projection=projection, return_fig=True)
interferogram.visualize_interferogram(projection=projection, fig=fig, globe=globe, ax=ax)

# Plot the displacements
fig, ax, globe = planet.visualize_topography(projection=projection, return_fig=True)
interferogram.visualize_displacements(projection=projection, fig=fig, globe=globe, ax=ax)

fm.clear()

# end of file
