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

# # Specify the baseline
# baseline = 4
# # Specify the basline uncertainty
# baseline_uncertainty = 0. #10
# # Specify the baseline orientation (roll=0 means horizontal)
# roll = 23.5
# # Specify the roll uncertainty
# roll_uncertainty = 0. #-8 / 3600
#
# # compute the interferogram
# interferogram = pet.dataAnalysis.surfaceDeformationAnalysis.simpleInterferogramRepeatPass(name="igram",
#                                                                                 instrument=instrument,
#                                                                                 planet=planet,
#                                                                                 track=track,
#                                                                                 campaign=campaign,
#                                                                                 deformation_map=deformation_map,
#                                                                                 baseline=baseline,
#                                                                                 baseline_uncertainty=baseline_uncertainty,
#                                                                                 roll=roll,
#                                                                                 roll_uncertainty=roll_uncertainty,
#                                                                                 time_offset_second_acquisition=
#                                                                                           orbit_cycle_time)
# # Calculate interferogram
# interferogram.calculate_igram()
#
# # Save interferogram
# interferogram.save_forward(file_name="/home/user/Documents/GitHub/"
#                              "Planetary-Exploration-Tool/files/igram_base_forward_1")
#
# # Calculate phase inverse
# interferogram.get_igram_inversion(use_dem="ellipsoid")
#
# # Save interferogram
# interferogram.save_inverse(file_name="/home/user/Documents/GitHub/"
#                              "Planetary-Exploration-Tool/files/igram_base_inverse_1")

# Load the interferogram
interferogram = pet.dataAnalysis.surfaceDeformationAnalysis.simpleInterferogramRepeatPass.from_file(instrument=
                                                                                                    instrument,
                                                                                          planet=planet,
                                                                                          campaign=campaign,
                                                                                          deformation_map=
                                                                                                    deformation_map,
                                                                                          forward_phase_file_name=
                                                                                                    "/home/user/"
                                                                                                    "Documents/"
                                                                                                    "GitHub/"
                                                                                                    "Planetary"
                                                                                                    "-Exploration-Tool/"
                                                                                                    "files/igram_"
                                                                                                    "base_forward_1",
                                                                                          inverse_phase_file_name=
                                                                                                    "/home/user/"
                                                                                                    "Documents/"
                                                                                                    "GitHub/"
                                                                                                    "Planetary"
                                                                                                    "-Exploration-Tool/"
                                                                                                    "files/igram_"
                                                                                                    "base_inverse_1")

# Define a projection
projection = pet.projections.biaxialProjections.biaxialCylindrical(name="biaxial cylindrical",
                                                                   folder_path=
                                                                   "/home/user/Documents/GitHub/Planetary"
                                                                   "-Exploration-Tool/figs/")

# Plot interferogram
interferogram.visualize_interferogram(projection=projection)

# Plot the displacements
interferogram.visualize_displacements(projection=projection)

fm.clear()

# end of file
