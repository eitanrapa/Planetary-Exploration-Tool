#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved
import numpy as np
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
deformation_map = pet.geophysical.deformationMap(name="base",
                                                 displacement_data_path=
                                                 "/home/user/Documents/GitHub/"
                                                 "Planetary-Exploration-Tool/input/Simulation_Base_Results.hdf5",
                                                 planet=planet)

# Make a bunch of interferograms
interferograms = []

track1 = pet.operations.track(start_time=times[0], end_time=times[1], planet=planet, instrument=instrument,
                              conops=conops, temporal_resolution=20, spatial_resolution=2000)
track1.load()

start_time = track1.start_time
end_time = track1.end_time
times_values = track1.data["time"].values

track2 = pet.operations.track(start_time=times[0], end_time=times[1], planet=planet, instrument=instrument,
                              conops=conops, temporal_resolution=20, spatial_resolution=2000)
track2.load()

for i in range(10):
    track1.start_time = start_time + orbit_cycle_time * i
    track1.end_time = end_time + orbit_cycle_time * i
    track1.data["time"].values = times_values + orbit_cycle_time * i

    track2.start_time = start_time + orbit_cycle_time * (i + 1)
    track2.end_time = end_time + orbit_cycle_time * (i + 1)
    track2.data["time"].values = times_values + orbit_cycle_time * (i + 1)

    interferogram = pet.operations.interferogram(instrument=instrument, planet=planet,
                                                 deformation_map=deformation_map,
                                                 track1=track1, track2=track2, conops=conops, baseline=10)
    # interferogram.calculate_igram()
    # interferogram.save()
    interferogram.load()
    interferograms.append(interferogram)

# Make a time series object
time_series = pet.operations.timeSeries(planet=planet, conops=conops, instrument=instrument)

# # Read the amplitudes, phases, and topographical errors
# amplitudes, amplitude_uncertainties, phases, phases_uncertainties, zs = time_series.create_1d_time_series(
#     interferograms=interferograms, spatial_points=np.asarray([[0, 90.0]]))
#
# print(amplitudes, amplitude_uncertainties, phases, phases_uncertainties, zs)

track1 = pet.operations.track(start_time=times[1], end_time=times[2], planet=planet, instrument=instrument,
                              conops=conops, temporal_resolution=20, spatial_resolution=2000)
track1.load()

start_time = track1.start_time
end_time = track1.end_time
times_values = track1.data["time"].values

track2 = pet.operations.track(start_time=times[1], end_time=times[2], planet=planet, instrument=instrument,
                              conops=conops, temporal_resolution=20, spatial_resolution=2000)
track2.load()

for i in range(10):
    track1.start_time = start_time + orbit_cycle_time * i
    track1.end_time = end_time + orbit_cycle_time * i
    track1.data["time"].values = times_values + orbit_cycle_time * i

    track2.start_time = start_time + orbit_cycle_time * (i + 1)
    track2.end_time = end_time + orbit_cycle_time * (i + 1)
    track2.data["time"].values = times_values + orbit_cycle_time * (i + 1)

    interferogram = pet.operations.interferogram(instrument=instrument, planet=planet, deformation_map=deformation_map,
                                                 track1=track1, track2=track2, conops=conops, baseline=10)
    # interferogram.calculate_igram()
    # interferogram.save()
    interferogram.load()
    interferograms.append(interferogram)

track1 = pet.operations.track(start_time=times[2], end_time=times[3], planet=planet, instrument=instrument,
                              conops=conops, temporal_resolution=20, spatial_resolution=2000)
track1.load()

start_time = track1.start_time
end_time = track1.end_time
times_values = track1.data["time"].values

track2 = pet.operations.track(start_time=times[2], end_time=times[3], planet=planet, instrument=instrument,
                              conops=conops, temporal_resolution=20, spatial_resolution=2000)
track2.load()

for i in range(10):
    track1.start_time = start_time + orbit_cycle_time * i
    track1.end_time = end_time + orbit_cycle_time * i
    track1.data["time"].values = times_values + orbit_cycle_time * i

    track2.start_time = start_time + orbit_cycle_time * (i + 1)
    track2.end_time = end_time + orbit_cycle_time * (i + 1)
    track2.data["time"].values = times_values + orbit_cycle_time * (i + 1)

    interferogram = pet.operations.interferogram(instrument=instrument, planet=planet, deformation_map=deformation_map,
                                                 track1=track1, track2=track2, conops=conops, baseline=10)
    # interferogram.calculate_igram()
    # interferogram.save()
    interferogram.load()
    interferograms.append(interferogram)

track1 = pet.operations.track(start_time=times[3], end_time=times[4], planet=planet, instrument=instrument,
                              conops=conops, temporal_resolution=20, spatial_resolution=2000)
track1.load()

start_time = track1.start_time
end_time = track1.end_time
times_values = track1.data["time"].values

track2 = pet.operations.track(start_time=times[3], end_time=times[4], planet=planet, instrument=instrument,
                              conops=conops, temporal_resolution=20, spatial_resolution=2000)
track2.load()

for i in range(10):
    track1.start_time = start_time + orbit_cycle_time * i
    track1.end_time = end_time + orbit_cycle_time * i
    track1.data["time"].values = times_values + orbit_cycle_time * i

    track2.start_time = start_time + orbit_cycle_time * (i + 1)
    track2.end_time = end_time + orbit_cycle_time * (i + 1)
    track2.data["time"].values = times_values + orbit_cycle_time * (i + 1)

    interferogram = pet.operations.interferogram(instrument=instrument, planet=planet, deformation_map=deformation_map,
                                                 track1=track1, track2=track2, conops=conops, baseline=10)
    # interferogram.calculate_igram()
    # interferogram.save()
    interferogram.load()
    interferograms.append(interferogram)

track1 = pet.operations.track(start_time=times[4], end_time=times[5], planet=planet, instrument=instrument,
                              conops=conops, temporal_resolution=20, spatial_resolution=2000)
track1.load()

start_time = track1.start_time
end_time = track1.end_time
times_values = track1.data["time"].values

track2 = pet.operations.track(start_time=times[4], end_time=times[5], planet=planet, instrument=instrument,
                              conops=conops, temporal_resolution=20, spatial_resolution=2000)
track2.load()

for i in range(10):
    track1.start_time = start_time + orbit_cycle_time * i
    track1.end_time = end_time + orbit_cycle_time * i
    track1.data["time"].values = times_values + orbit_cycle_time * i

    track2.start_time = start_time + orbit_cycle_time * (i + 1)
    track2.end_time = end_time + orbit_cycle_time * (i + 1)
    track2.data["time"].values = times_values + orbit_cycle_time * (i + 1)

    interferogram = pet.operations.interferogram(instrument=instrument, planet=planet, deformation_map=deformation_map,
                                                 track1=track1, track2=track2, conops=conops, baseline=10)
    # interferogram.calculate_igram()
    # interferogram.save()
    interferogram.load()
    interferograms.append(interferogram)

lats = np.linspace(89.5, -89.5, 359)
lons = np.linspace(179.5, -179.5, 719)
lats, lons = np.meshgrid(lats, lons)

# Read the amplitudes, phases, and topographical errors
amplitudes, amplitude_uncertainties, phases, phases_uncertainties, zs = time_series.create_3d_time_series(
    interferograms=interferograms, spatial_points=np.asarray([lats.flatten(), lons.flatten()]).T)

# Save the amplitudes and phases
np.save("amplitudes.npy", amplitudes)
np.save("phase.npy", phases)

fm.clear()

# end of file
