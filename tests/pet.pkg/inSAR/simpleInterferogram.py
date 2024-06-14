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

ins = pet.instruments.nightingale(name="nightingale", body_id=-303, start_look_angle=25, end_look_angle=35,
                                  orbit_cycle=110020, wavelength=0.13)

planet = pet.planets.enceladus(name="enceladus")

displacements = pet.insar.displacementMap(name="base",
                                          displacement_data_path=path / "Simulation_Base_Results.hdf5")

interferogram = pet.insar.simpleInterferogram(name="igram", pairing_one=0, pairing_two=1, track_number=2)

interferogram.calculate_igram(instrument=ins, planet=planet, displacements=displacements, time_interval=60,
                              ground_resoluton=2, baseline=0.01)

# interferogram.visualize()

fm.clear()

# end of file
