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

# Make a con ops
conops = pet.conOps.nightingale5to1(name="nightingale",
                                    body_id=-303, start_time="2046 DEC 20 15:10:40.134", planet=planet)

# Make a projection
projection = pet.projections.biaxialPlanar(name="biaxial planar", central_latitude=-90, north_extent=-30,
                                           folder_path="/home/user/Documents/GitHub/Planetary-Exploration-Tool/figs")

# Create a time conversion instance
time_conversion = pet.spicetoolkit.timeConversion()

times = conops.get_five_tracks()

# convert to utcs
utcs = time_conversion.convert_utcs(ets=times)

# Plot the orbits
conops.plot_orbit(start_time=utcs[0], end_time=utcs[5], projection=projection)

# utc = "2046 DEC 20 15:10:40.134"
#
# et = time_conversion.convert_ets(utcs=utc)
#
# position, velocity = conops.get_states(times=et)
#
# # Plot the orbit
# conops.plot_orbit(start_time="2046 DEC 20 15:10:40.134", end_time="2046 DEC 22 03:10:40.134",
#                   projection=projection)

fm.clear()

# end of file
