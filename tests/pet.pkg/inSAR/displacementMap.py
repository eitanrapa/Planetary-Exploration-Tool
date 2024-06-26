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

# Make a displacement map
displacements = pet.insar.displacementMap(
    name="het", displacement_data_path="/home/eitanrapa/Documents/projects/other/Simulation_Het_Results.hdf5",
    planet=planet)

# Define a projection
projection = pet.projections.biaxialPlanar(name="biaxial planar", central_latitude=-90, north_extent=-30,
                                             folder_path="/home/eitanrapa/Documents/projects/pet/figs")

# Visualize displacements
displacements.visualize(projection=projection, time_point=0, direction="east")

# Visualize displacements
displacements.visualize(projection=projection, time_point=0, direction="north")

# Visualize displacements
displacements.visualize(projection=projection, time_point=0, direction="up")

# end of file
