#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import numpy as np
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

# First track
track = pet.dataAcquisition.track.from_file(planet=planet, campaign=campaign, instrument=instrument,
                                            file_name="/home/user/Documents/GitHub/"
                                                       "Planetary-Exploration-Tool/files/track1")

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

# # Specify the baseline
# baseline = 4
# # Specify the basline uncertainty
# baseline_uncertainty = 0. #10
# # Specify the baseline orientation (roll=0 means horizontal)
# roll = 23.5
# # Specify the roll uncertainty
# roll_uncertainty = 0. #-8 / 3600
#
# for i in range(10):
#
#     # compute the interferogram
#     interferogram = pet.dataAnalysis.surfaceDeformationAnalysis.simpleInterferogramRepeatPass(name="igram",
#                                                                                     instrument=instrument,
#                                                                                     planet=planet,
#                                                                                     track=track,
#                                                                                     campaign=campaign,
#                                                                                     deformation_map=deformation_map,
#                                                                                     baseline=baseline,
#                                                                                     baseline_uncertainty=
#                                                                                               baseline_uncertainty,
#                                                                                     roll=roll,
#                                                                                     roll_uncertainty=roll_uncertainty,
#                                                                                     time_offset_first_acquisition=
#                                                                                                 orbit_cycle_time*i,
#                                                                                     time_offset_second_acquisition=
#                                                                                               orbit_cycle_time*(i + 1))
#     # Calculate interferogram
#     interferogram.calculate_igram()
#
#     # Save interferogram
#     interferogram.save_forward(file_name="/home/user/Documents/GitHub/"
#                                  "Planetary-Exploration-Tool/files/igram_base_forward_{}".format(i))
#
#     # Calculate phase inverse
#     interferogram.get_igram_inversion(use_dem="ellipsoid")
#
#     # Save interferogram
#     interferogram.save_inverse(file_name="/home/user/Documents/GitHub/"
#                                  "Planetary-Exploration-Tool/files/igram_base_inverse_{}".format(i))

path = "/home/user/Documents/GitHub/Planetary-Exploration-Tool/files/"
interferogram_forward_files = [path + 'igram_base_forward_{}'.format(i) for i in range(10)]
interferogram_inverse_files = [path + 'igram_base_inverse_{}'.format(i) for i in range(10)]
interferograms = pet.dataAnalysis.surfaceDeformationAnalysis.simpleInterferogramRepeatPass.from_files(
                                                                                            instrument=instrument,
                                                                                            planet=planet,
                                                                                            deformation_map=
                                                                                            deformation_map,
                                                                                            campaign=campaign,
                                                                                            forward_file_list=
                                                                                            interferogram_forward_files,
                                                                                            inverse_file_list=
                                                                                            interferogram_inverse_files)

# Create a grid of latitudes and longitudes
lats = np.linspace(89.5, -89.5, 359)
lons = np.linspace(179.5, -179.5, 719)
lats, lons = np.meshgrid(lats, lons)

# Make a time series object
time_series = pet.dataAnalysis.surfaceDeformationAnalysis.surfaceDeformation1DTimeSeries(name="1d time series",
                                                                                         planet=planet,
                                                                                         campaign=campaign,
                                                                                         instrument=instrument)

# Read the amplitudes, phases, and topographical errors
time_series.create_1d_time_series(interferograms=interferograms,
                                  spatial_points=np.asarray([lats.flatten(), lons.flatten()]).T, processors=3)

time_series.save(file_name="/home/user/Documents/GitHub/Planetary-Exploration-Tool/files/time_series_1d_base")

time_series = (
    pet.dataAnalysis.surfaceDeformationAnalysis.surfaceDeformation1DTimeSeries.from_file(
        planet=planet,
        campaign=campaign,
        instrument=
        instrument,
        file_name="/home/user/Documents/"
                  "GitHub/"
                  "Planetary-Exploration-Tool/"
                  "files/time_series_1d_base"))

# Define a projection
projection = pet.projections.biaxialProjections.biaxialCylindrical(name="biaxial cylindrical",
                                                                   folder_path="/home/user/Documents/GitHub/"
                                                                   "Planetary-Exploration-Tool/figs/")

time_series.visualize_time_series_amplitudes(projection=projection)
time_series.visualize_time_series_phases(projection=projection)

fm.clear()

# end of file
