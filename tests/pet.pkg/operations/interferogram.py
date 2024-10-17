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

# Make a displacement map
deformation_map = pet.models.deformationMap(name="het",
                                            displacement_data_path=
                                            "/home/user/Documents/other/Simulation_Het_Results.hdf5",
                                            planet=planet)
# Get the times defining the first five tracks
times = instrument.get_five_tracks()

# Get the orbit cycle time of the instrument
orbit_cycle_time = instrument.orbit_cycle

track1 = pet.mission.track(start_time=times[0], end_time=times[1], planet=planet, instrument=instrument,
                           temporal_resolution=20, spatial_resolution=2000)
time_space1, swath1 = track1.calculate_ground_swath()
time_space2 = time_space1 + orbit_cycle_time

# Make an interferogram
interferogram = pet.operations.interferogram(instrument=instrument, planet=planet, deformation_map=deformation_map,
                                             swath1=swath1, swath2=swath1, time_space1=time_space1,
                                             time_space2=time_space2, baseline=10)
interferogram.calculate_igram()

# # Save interferogram
# interferogram.save(path="/home/user/Documents/GitHub/Planetary-Exploration-Tool/files")

# # Load interferogram
# interferogram = pet.insar.simpleInterferogram(instrument=instrument, planet=planet,
#                                               deformation_map=deformation_map,
#                                               load_path="/files/base_0_0_1_interferogram.hdf5")

# Define a projection
projection = pet.projections.biaxialPlanar(name="biaxial planar", central_latitude=-90, north_extent=-30,
                                           folder_path="/home/user/Documents/GitHub/Planetary-Exploration-Tool/figs")

# Plot interferogram
interferogram.visualize_interferogram(projection=projection)

# Plot the displacements

interferogram.visualize_displacements(projection=projection)

fm.clear()

# end of file
