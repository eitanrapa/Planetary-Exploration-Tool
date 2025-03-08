#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet

# Create a file manager
fm = pet.spiceTools.fileManager(folder_path="/home/user/Documents/GitHub/Planetary-Exploration-Tool/input")

# Furnish some files
fm.furnsh(names_list=["cas_enceladus_ssd_spc_1024icq_v1.bds", "pck00011_n0066.tpc",
                      "insar_6stride_26d_v7_seo.bsp", "latest_leapseconds.tls"])

# Create a time conversion object
time_conversion = pet.spiceTools.timeConversion()

# Define a time
utc = "2046 DEC 20 15:10:40.134"

# Convert to utc
et = time_conversion.convert_ets(utcs=utc)

# Convert back
utc_converted = time_conversion.convert_utcs(ets=et)

print(utc, utc_converted)

# Define a time
utcs = ["2046 DEC 20 15:10:40.134", "2046 DEC 21 15:10:40.134"]

# Convert to utc
ets = time_conversion.convert_ets(utcs=utcs)

# Convert back
times_converted = time_conversion.convert_utcs(ets=ets)

print(utcs, times_converted)

# end of file
