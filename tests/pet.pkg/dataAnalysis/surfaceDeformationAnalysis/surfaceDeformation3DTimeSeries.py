#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import threading
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

#
# track1 = pet.dataAcquisition.track.from_file(planet=planet, campaign=campaign, instrument=instrument,
#                                              file_name="/home/user/Documents/GitHub/"
#                                                        "Planetary-Exploration-Tool/files/track1")
#
# track2 = pet.dataAcquisition.track.from_file(planet=planet, campaign=campaign, instrument=instrument,
#                                              file_name="/home/user/Documents/GitHub/"
#                                                        "Planetary-Exploration-Tool/files/track1")
#
# for i in range(10):
#     track1.modify_time(orbit_cycle_time * i)
#
#     track2.modify_time(orbit_cycle_time * (i + 1))
#
#     interferogram = pet.dataAnalysis.surfaceDeformationAnalysis.simpleInterferogram(name="igram".format(i),
#                                                                                     instrument=instrument,
#                                                                                     planet=planet,
#                                                                                     deformation_map=deformation_map,
#                                                                                     track1=track1, track2=track2,
#                                                                                     campaign=campaign, baseline=20)
#     interferogram.calculate_igram()
#
#     # Save interferogram
#     interferogram.save(file_name="/home/user/Documents/GitHub/"
#                                  "Planetary-Exploration-Tool/files/igram_base_{}".format(i))

# track1 = pet.dataAcquisition.track.from_file(planet=planet, campaign=campaign, instrument=instrument,
#                                              file_name="/home/user/Documents/GitHub/"
#                                                        "Planetary-Exploration-Tool/files/track2")
#
# track2 = pet.dataAcquisition.track.from_file(planet=planet, campaign=campaign, instrument=instrument,
#                                              file_name="/home/user/Documents/GitHub/"
#                                                        "Planetary-Exploration-Tool/files/track2")
#
# for i in range(10):
#     track1.modify_time(orbit_cycle_time * i)
#
#     track2.modify_time(orbit_cycle_time * (i + 1))
#
#     interferogram = pet.dataAnalysis.surfaceDeformationAnalysis.simpleInterferogram(name="igram".format(i),
#                                                                                     instrument=instrument,
#                                                                                     planet=planet,
#                                                                                     deformation_map=deformation_map,
#                                                                                     track1=track1, track2=track2,
#                                                                                     campaign=campaign, baseline=20)
#     interferogram.calculate_igram()
#
#     # Save interferogram
#     interferogram.save(file_name="/home/user/Documents/GitHub/"
#                                  "Planetary-Exploration-Tool/files/igram_base_1{}".format(i))

# track1 = pet.dataAcquisition.track.from_file(planet=planet, campaign=campaign, instrument=instrument,
#                                              file_name="/home/user/Documents/GitHub/"
#                                                        "Planetary-Exploration-Tool/files/track3")
#
# track2 = pet.dataAcquisition.track.from_file(planet=planet, campaign=campaign, instrument=instrument,
#                                              file_name="/home/user/Documents/GitHub/"
#                                                        "Planetary-Exploration-Tool/files/track3")
#
# for i in range(10):
#     track1.modify_time(orbit_cycle_time * i)
#
#     track2.modify_time(orbit_cycle_time * (i + 1))
#
#     interferogram = pet.dataAnalysis.surfaceDeformationAnalysis.simpleInterferogram(name="igram".format(i),
#                                                                                     instrument=instrument,
#                                                                                     planet=planet,
#                                                                                     deformation_map=deformation_map,
#                                                                                     track1=track1, track2=track2,
#                                                                                     campaign=campaign, baseline=20)
#     interferogram.calculate_igram()
#
#     # Save interferogram
#     interferogram.save(file_name="/home/user/Documents/GitHub/"
#                                  "Planetary-Exploration-Tool/files/igram_base_2{}".format(i))
#
# track1 = pet.dataAcquisition.track.from_file(planet=planet, campaign=campaign, instrument=instrument,
#                                              file_name="/home/user/Documents/GitHub/"
#                                                        "Planetary-Exploration-Tool/files/track4")
#
# track2 = pet.dataAcquisition.track.from_file(planet=planet, campaign=campaign, instrument=instrument,
#                                              file_name="/home/user/Documents/GitHub/"
#                                                        "Planetary-Exploration-Tool/files/track4")
#
# for i in range(10):
#     track1.modify_time(orbit_cycle_time * i)
#
#     track2.modify_time(orbit_cycle_time * (i + 1))
#
#     interferogram = pet.dataAnalysis.surfaceDeformationAnalysis.simpleInterferogram(name="igram".format(i),
#                                                                                     instrument=instrument,
#                                                                                     planet=planet,
#                                                                                     deformation_map=deformation_map,
#                                                                                     track1=track1, track2=track2,
#                                                                                     campaign=campaign, baseline=20)
#     interferogram.calculate_igram()
#
#     # Save interferogram
#     interferogram.save(file_name="/home/user/Documents/GitHub/"
#                                  "Planetary-Exploration-Tool/files/igram_base_3{}".format(i))

# track1 = pet.dataAcquisition.track.from_file(planet=planet, campaign=campaign, instrument=instrument,
#                                              file_name="/home/user/Documents/GitHub/"
#                                                        "Planetary-Exploration-Tool/files/track5")
#
# track2 = pet.dataAcquisition.track.from_file(planet=planet, campaign=campaign, instrument=instrument,
#                                              file_name="/home/user/Documents/GitHub/"
#                                                        "Planetary-Exploration-Tool/files/track5")
#
# for i in range(10):
#     track1.modify_time(orbit_cycle_time * i)
#
#     track2.modify_time(orbit_cycle_time * (i + 1))
#
#     interferogram = pet.dataAnalysis.surfaceDeformationAnalysis.simpleInterferogram(name="igram".format(i),
#                                                                                     instrument=instrument,
#                                                                                     planet=planet,
#                                                                                     deformation_map=deformation_map,
#                                                                                     track1=track1, track2=track2,
#                                                                                     campaign=campaign, baseline=20)
#     interferogram.calculate_igram()
#
#     # Save interferogram
#     interferogram.save(file_name="/home/user/Documents/GitHub/"
#                                  "Planetary-Exploration-Tool/files/igram_base_4{}".format(i))

path = "/home/user/Documents/GitHub/Planetary-Exploration-Tool/files/"
interferogram_files = [path + 'igram_base_{}'.format(i) for i in range(50)]
# interferogram_files = [path + 'igram_base_0', path + 'igram_base_1',
#                        path + 'igram_base_10', path + 'igram_base_11',
#                        path + 'igram_base_20', path + 'igram_base_21',
#                        path + 'igram_base_30', path + 'igram_base_31',
#                        path + 'igram_base_40', path + 'igram_base_41']
interferograms = pet.dataAnalysis.surfaceDeformationAnalysis.simpleInterferogram.from_files(instrument=instrument,
                                                                                            planet=planet,
                                                                                            deformation_map=
                                                                                            deformation_map,
                                                                                            campaign=campaign,
                                                                                            file_list=
                                                                                            interferogram_files)


lats = np.linspace(89.5, -89.5, 359)
lons = np.linspace(179.5, -179.5, 719)
lats, lons = np.meshgrid(lats, lons)

# # Make a time series object
# time_series = pet.dataAnalysis.surfaceDeformationAnalysis.surfaceDeformation3DTimeSeries(name="3d time series",
#                                                                                          planet=planet,
#                                                                                          campaign=campaign,
#                                                                                          instrument=instrument)
#
# # Read the amplitudes, phases, and topographical errors
# time_series.create_3d_time_series(interferograms=interferograms,
#                                   spatial_points=np.asarray([lats.flatten(), lons.flatten()]).T, processors=3)
#
# time_series.save(file_name="/home/user/Documents/GitHub/Planetary-Exploration-Tool/files/time_series_3d_base")

time_series = (
    pet.dataAnalysis.surfaceDeformationAnalysis.surfaceDeformation3DTimeSeries.from_file(
        planet=planet,
        campaign=campaign,
        instrument=
        instrument,
        file_name="/home/user/Documents/"
                  "GitHub/"
                  "Planetary-Exploration-Tool/"
                  "files/time_series_3d_base"))

# Define a projection
projection = pet.projections.biaxialProjections.biaxialPlanar(name="biaxial planar", north_extent=-30,
                                                              central_latitude=-90,
                                                              folder_path="/home/user/Documents/GitHub/"
                                                                   "Planetary-Exploration-Tool/figs/")

time_series.visualize_time_series_amplitudes(projection=projection, direction="east")
time_series.visualize_time_series_phases(projection=projection, direction="east")

fm.clear()

# end of file
