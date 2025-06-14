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

print([axis for axis in planet.get_axes()])

# Make a projection
projection = pet.projections.biaxialProjections.biaxialPlanar(name="biaxial planar",
                                                              central_latitude=-90, north_extent=-30,
                                                              folder_path="/home/user/Documents/"
                                                                          "GitHub/Planetary-Exploration-Tool/figs")

# Visualize
planet.visualize_topography(projection=projection)

fm.clear()

# end of file
