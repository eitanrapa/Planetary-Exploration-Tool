#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

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

# # # Get the times defining the first five tracks
# times = campaign.get_five_tracks()
#
# # First track
# track = pet.dataAcquisition.track(name="track1", start_time=times[0], end_time=times[1], planet=planet,
#                                   campaign=campaign, instrument=instrument, spatial_resolution="2*km",
#                                   temporal_resolution="1*minute")
#
# # Calculate the positions
# track.calculate_ground_swath()

# # Save the track
# track.save(file_name="/home/user/Documents/GitHub/Planetary-Exploration-Tool/files/track1")
#
# # Second track
# track = pet.dataAcquisition.track(start_time=times[1], end_time=times[2], planet=planet,
#                                   campaign=campaign, instrument=instrument)
#
# # Calculate the positions
# track.calculate_ground_swath()
#
# # Save the track
# track.save(file_name="/home/user/Documents/GitHub/Planetary-Exploration-Tool/files/track2")
#
# # Third track
# track = pet.dataAcquisition.track(start_time=times[2], end_time=times[3], planet=planet,
#                                   campaign=campaign, instrument=instrument)
#
# # Calculate the positions
# track.calculate_ground_swath()
#
# # Save the track
# track.save(file_name="/home/user/Documents/GitHub/Planetary-Exploration-Tool/files/track3")
#
# # Fourth track
# track = pet.dataAcquisition.track(start_time=times[3], end_time=times[4], planet=planet,
#                                   campaign=campaign, instrument=instrument)
#
# # Calculate the positions
# track.calculate_ground_swath()
#
# # Save the track
# track.save(file_name="/home/user/Documents/GitHub/Planetary-Exploration-Tool/files/track4")
#
# # Fifth track
# track = pet.dataAcquisition.track(start_time=times[4], end_time=times[5], planet=planet,
#                                   campaign=campaign, instrument=instrument)

# # Calculate the positions
# track.calculate_ground_swath()
#
# # Save the track
# track.save(file_name="/home/user/Documents/GitHub/Planetary-Exploration-Tool/files/track5")

track1 = pet.dataAcquisition.track.from_file(planet=planet, campaign=campaign, instrument=instrument,
                                             file_name="/home/user/Documents/GitHub/Planetary-Exploration-Tool/"
                                                       "files/track1")

track2 = pet.dataAcquisition.track.from_file(planet=planet, campaign=campaign, instrument=instrument,
                                             file_name="/home/user/Documents/GitHub/Planetary-Exploration-Tool/"
                                                       "files/track2")

track3 = pet.dataAcquisition.track.from_file(planet=planet, campaign=campaign, instrument=instrument,
                                             file_name="/home/user/Documents/GitHub/Planetary-Exploration-Tool/"
                                                       "files/track3")

track4 = pet.dataAcquisition.track.from_file(planet=planet, campaign=campaign, instrument=instrument,
                                             file_name="/home/user/Documents/GitHub/Planetary-Exploration-Tool/"
                                                       "files/track4")

track5 = pet.dataAcquisition.track.from_file(planet=planet, campaign=campaign, instrument=instrument,
                                             file_name="/home/user/Documents/GitHub/Planetary-Exploration-Tool/"
                                                       "files/track5")

# Define a projection
projection = pet.projections.biaxialProjections.biaxialCylindrical(name="biaxial cylindrical",
                                                                   folder_path="/home/user/Documents/GitHub/"
                                                                   "Planetary-Exploration-Tool/figs")
# Plot the track
fig, ax, globe = planet.visualize_topography(projection=projection, return_fig=True)

track1.visualize_swath(projection=projection, fig=fig, ax=ax, globe=globe)

fm.clear()

# end of file
