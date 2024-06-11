#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pathlib
import cspyce as spice
import pet

path = pathlib.PosixPath("/home/eitanrapa/Documents/projects/other")
spice.furnsh(path / "cas_enceladus_ssd_spc_1024icq_v1.bds")  # Topography of Enceladus
spice.furnsh(path / "pck00011_n0066.tpc")  # Reference frames
spice.furnsh(path / "insar_6stride_26d_v7_seo.bsp")  # Ephemeris data
spice.furnsh(path / "latest_leapseconds.tls")  # Leap seconds file

ins = pet.instruments.nightingale(name="nightingale", body_id=-303, start_look_angle=25, end_look_angle=35)

planet = pet.planets.enceladus(name="enceladus")

projection = pet.projections.biaxialPlanar(name="biaxial planar", central_latitude=-90)

ins.plot_orbit(start_time="2046 DEC 20 15:10:40.134", end_time="2046 DEC 22 03:10:40.134", planet=planet,
               projection=projection, north_extent=-30)

# end of file
