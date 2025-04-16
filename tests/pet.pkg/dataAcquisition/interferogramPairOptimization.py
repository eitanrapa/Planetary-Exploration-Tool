#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import pet
import json

# Create a file manager
fm = pet.spiceTools.fileManager(folder_path="/input")

# Furnish some files
fm.furnsh(names_list=["cas_enceladus_ssd_spc_1024icq_v1.bds", "pck00011_n0066.tpc",
                      "insar_6stride_26d_v7_seo.bsp", "latest_leapseconds.tls"])

# Make a planet
planet = pet.planets.enceladus(name="enceladus")

# Make a con ops
campaign = pet.campaign.nightingale5to1(name="nightingale",
                                    body_id=-303, start_time="2046 DEC 20 15:10:40.134", planet=planet)

# Load the acquisition cadence file
with open('/input/acquisition_cadence.json', 'r') as file:
    acquisition_cadence = json.load(file)

# Make an onboard processing object
onboard_processing = pet.simulations.onboardProcessing(planet=planet, campaign=campaign, acquisition_radius=400/3,
                                                       acquisition_cadence=acquisition_cadence)

# Get the acquisition times and baselines from the simulation
acquisition_times, baselines = onboard_processing.simulate_acquisitions(n_acquisitions_stored=5,
                                                                        downlink_pair_capacity=1,
                                                                        total_time=2.246e6, distribution="gaussian")

print(acquisition_times, baselines)

# end of file
