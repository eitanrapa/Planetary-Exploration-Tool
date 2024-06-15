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

# Make an instrument
ins = pet.instruments.nightingale(name="nightingale", body_id=-303, start_look_angle=25, end_look_angle=35)

# Make a planet
planet = pet.planets.enceladus(name="enceladus")

# Make a projection
projection = pet.projections.biaxialPlanar(name="biaxial planar", central_latitude=-90)

# Make a visualization tool
visualization = pet.visualization.cartopyViz(name="cartopy tool", north_extent=-30,
                                             folder_path="/home/eitanrapa/Documents/projects/pet/figs")

# Plot the orbit
ins.plot_orbit(visualization=visualization, start_time="2046 DEC 20 15:10:40.134", end_time="2046 DEC 22 03:10:40.134",
               planet=planet, projection=projection)

fm.clear()

# end of file
