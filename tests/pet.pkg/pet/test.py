#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet

# Create a file manager
fm = pet.spicetoolkit.filemanager(folder_path="/home/eitanrapa/Documents/projects/other")

# Furnish some files
fm.furnsh(names_list=["cas_enceladus_ssd_spc_1024icq_v1.bds", "pck00011_n0066.tpc",
                      "insar_6stride_26d_v7_seo.bsp", "latest_leapseconds.tls"])

instrument = pet.instruments.nightingale(name="nightingale", body_id=-303, start_look_angle=25, end_look_angle=35)

planet = pet.planets.enceladus(name="enceladus")

gs = pet.insar.groundSwath(name="1", start_time="2046 DEC 20 15:10:40.134", end_time="2046 DEC 20 15:10:50.134",
                           temporal_resolution=10, spatial_resolution=2000, planet=planet, instrument=instrument)

times = instrument.get_five_tracks(planet)

fm.clear()

# end of file
