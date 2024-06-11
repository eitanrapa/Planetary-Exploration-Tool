#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import cspyce as spice
import pet
import pathlib

path = pathlib.PosixPath("/home/eitanrapa/Documents/projects/other")
spice.furnsh(path / "cas_enceladus_ssd_spc_1024icq_v1.bds")  # Topography of Enceladus
spice.furnsh(path / "pck00011_n0066.tpc")  # Reference frames

planet = pet.planets.enceladus(name="enceladus")

projection = pet.projections.biaxialConic(name="biaxial conic")

planet.visualize_topography(projection=projection)

spice.kclear()

# end of file