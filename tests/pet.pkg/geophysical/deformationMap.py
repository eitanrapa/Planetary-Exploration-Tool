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

# Make a displacement map
deformation_map = pet.geophysical.deformationMap(name="het",
                                                 displacement_data_path=
                                                 "/home/user/Documents/GitHub/Planetary-Exploration-Tool/"
                                                 "input/Simulation_Het_Results.hdf5",
                                                 planet=planet)
# Make a con ops
conops = pet.conOps.nightingale5to1(name="nightingale",
                                    body_id=-303, start_time="2046 DEC 20 15:10:40.134", planet=planet)

# Define a projection
projection = pet.projections.biaxialPlanar(name="biaxial planar", central_latitude=-90, north_extent=-30,
                                           folder_path="/home/user/Documents/GitHub/Planetary-Exploration-Tool/figs")

# Visualize displacements
deformation_map.visualize(projection=projection, time_point=0, direction="east")

# Visualize displacements
deformation_map.visualize(projection=projection, time_point=0, direction="north")

# Visualize displacements
deformation_map.visualize(projection=projection, time_point=0 , direction="up")

# # Visualize displacements
# displacements.visualize(projection=projection, time_point=1482219658.317582 + 110020, direction="east")
#
# # Visualize displacements
# displacements.visualize(projection=projection, time_point=1482219658.317582 + 110020, direction="north")
#
# # Visualize displacements
# displacements.visualize(projection=projection, time_point=1482219658.317582 + 110020, direction="up")

# end of file
