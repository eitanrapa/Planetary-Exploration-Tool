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

# gs = pet.insar.groundSwath(name="1", start_time="2046 DEC 20 15:10:40.134", end_time="2046 DEC 20 15:10:50.134",
#                              time_interval=10, ground_resolution=0.2)

ins = pet.instruments.nightingale(name="nightingale", body_id=-303, start_look_angle=25, end_look_angle=35)

planet = pet.planets.enceladus(name="enceladus")

times = ins.get_five_tracks(planet)

# planet.visualize_topography()

# gs.calculate_swath(instrument=ins, planet=planet)

displacements = pet.insar.displacementMap(name="base", displacement_data_path=path / "Simulation_Base_Results.hdf5",
                                          planet=planet)

# displacements.attach(swath=gs)

# gs.visualize()

spice.kclear()

# import matplotlib.pyplot as plt
#
# # Plot X,Y,Z
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# x = []
# y = []
# z = []
# for beam in swath:
#     for point in beam:
#         x.append(point.x)
#         y.append(point.y)
#         z.append(point.z)
# ax.plot_trisurf(x, y, z, color='white', edgecolor='none', linewidth=0, antialiased=False, alpha=0.5)
# ax.scatter(x, y, z, c='red')
# plt.show()

# end of file
