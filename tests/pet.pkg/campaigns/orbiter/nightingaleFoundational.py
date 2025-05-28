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

# Make a con ops
campaign = pet.campaigns.orbiter.nightingaleFoundational(name="nightingale",
                                                 body_id=-303, start_time="2046 DEC 20 15:10:40.134", planet=planet)

# Create a time conversion instance
time_conversion = pet.spiceTools.timeConversion()

time = campaign.get_et_start_time()

print(time)

fm.clear()

# end of file
