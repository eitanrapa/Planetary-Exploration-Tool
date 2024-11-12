#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet

# Create a file manager
fm = pet.spicetoolkit.fileManager(folder_path="/home/user/Documents/GitHub/Planetary-Exploration-Tool/input")

# Furnish some files
fm.furnsh(names_list=["cas_enceladus_ssd_spc_1024icq_v1.bds", "pck00011_n0066.tpc",
                      "insar_6stride_26d_v7_seo.bsp", "latest_leapseconds.tls"])

# Make a planet
planet = pet.planets.enceladus(name="enceladus")

# Make an instrument
instrument = pet.instruments.chirpchirp(
    instrument_parameters_path="/home/user/Documents/GitHub/Planetary-Exploration-Tool/input/Parameters Sband.xlsx")

# Make a con ops
conops = pet.conOps.nightingale5to1(name="nightingale",
                                    body_id=-303, start_time="2046 DEC 20 15:10:40.134", planet=planet)

# Get the times defining the first five tracks
times = conops.get_five_tracks()

# Get the orbit cycle time of the instrument
orbit_cycle_time = conops.orbit_cycle

# Make a displacement map
deformation_map = pet.geophysical.deformationMap(name="het",
                                                 displacement_data_path=
                                                 "/home/user/Documents/GitHub/"
                                                 "Planetary-Exploration-Tool/input/Simulation_Het_Results.hdf5",
                                                 planet=planet)

track1 = pet.mission.track(start_time=times[0], end_time=times[1], planet=planet, instrument=instrument,
                           conops=conops, temporal_resolution=20, spatial_resolution=2000)
track1.load()

track2 = pet.mission.track(start_time=times[0], end_time=times[1], planet=planet, instrument=instrument,
                           conops=conops, temporal_resolution=20, spatial_resolution=2000)
track2.load()
track2.start_time = track1.start_time + orbit_cycle_time
track2.end_time = track1.end_time + orbit_cycle_time
track2.data["time"].values = track1.data["time"].values + orbit_cycle_time

# Make an interferogram
interferogram = pet.operations.interferogram(instrument=instrument, planet=planet, deformation_map=deformation_map,
                                             track1=track1, track2=track2, conops=conops, baseline=10)
interferogram.calculate_igram()
# interferogram.load()

# Save interferogram
# interferogram.save()

# Define a projection
projection = pet.projections.biaxialPlanar(name="biaxial planar", central_latitude=-90, north_extent=-30,
                                           folder_path="/home/user/Documents/GitHub/Planetary-Exploration-Tool/figs")

# Plot interferogram
fig, ax, globe = planet.visualize_topography(projection=projection, return_fig=True)
interferogram.visualize_interferogram(projection=projection, fig=fig, globe=globe, ax=ax)

# Plot the displacements
# fig, ax, globe = planet.visualize_topography(projection=projection, return_fig=True)
# interferogram.visualize_displacements(projection=projection, fig=fig, globe=globe, ax=ax)

fm.clear()

# end of file
