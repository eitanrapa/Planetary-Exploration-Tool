#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np


# Create a file manager
fm = pet.spicetoolkit.fileManager(folder_path="/home/user/Documents/GitHub/Planetary-Exploration-Tool/input")

# Furnish some files
fm.furnsh(names_list=["cas_enceladus_ssd_spc_1024icq_v1.bds", "pck00011_n0066.tpc",
                      "insar_6stride_26d_v7_seo.bsp", "latest_leapseconds.tls"])
# Make a planet
planet = pet.planets.enceladus(name="enceladus")

# Make a displacement map
deformation_map = pet.geophysical.deformationMap(name="base",
                                                 displacement_data_path=
                                                 "/home/user/Documents/GitHub/Planetary-Exploration-Tool/"
                                                 "input/Simulation_Base_Results.hdf5",
                                                 planet=planet)
# Make a con ops
conops = pet.conOps.nightingale5to1(name="nightingale",
                                    body_id=-303, start_time="2046 DEC 20 15:10:40.134", planet=planet)

directions = ["east", "north", "up"]

for i in range(len(directions)):

    # Define a projection
    projection = pet.projections.biaxialPlanar(name="biaxial planar", central_latitude=-90, north_extent=30,
                                               folder_path="/home/user/Documents/GitHub/Planetary-Exploration-Tool/figs")

    fig, ax, globe = projection.proj(planet=planet)

    # Create a dense mesh grid of latitude, longitude points
    lats = np.linspace(89.5, -89.5, 359)
    lons = np.linspace(179.5, -179.5, 719)
    lats, lons = np.meshgrid(lats, lons)

    # a_fits, phi_fits = deformation_map.time_series(spatial_points=np.asarray([lats.flatten(), lons.flatten()]).T,
    #                                                direction=directions[i])
    # a_fits = np.asarray(a_fits)
    # phi_fits = np.asarray(phi_fits)

    # Load .npy file
    a_fits_recovered = np.load("/home/user/Documents/GitHub/Planetary-Exploration-Tool/tests/pet.pkg/operations/amplitudes.npy")
    phi_fits_recovered = np.load("/home/user/Documents/GitHub/Planetary-Exploration-Tool/tests/pet.pkg/operations/phase.npy")

    # Remove any element that corresponds to a latitude greater than -40 degrees
    # a_fits_recovered = a_fits_recovered[lats.flatten() < -40]
    # a_fits = a_fits[lats.flatten() < -40]

    # Get the first element of each row
    a_fits_recovered = a_fits_recovered[:, i]

    # # Plot the a_fits_recovered
    # im = ax.scatter(lons.flatten()[lats.flatten() < -40], lats.flatten()[lats.flatten() < -40],
    #                 c=a_fits_recovered - a_fits, cmap='viridis', transform=ccrs.PlateCarree(globe=globe))

    # Plot the a_fits_recovered
    im = ax.scatter(lons.flatten(), lats.flatten(),
                    c=a_fits_recovered, cmap='viridis', transform=ccrs.PlateCarree(globe=globe))

    # Add a colorbar
    plt.colorbar(im, ax=ax, orientation='vertical', shrink=0.25)

    # Add labels and legend
    plt.title('Displacements ' + directions[i], pad=20)

    # Save the plot
    plt.show()


# # Define a projection
# projection = pet.projections.biaxialPlanar(name="biaxial planar", central_latitude=-90, north_extent=-30,
#                                            folder_path="/home/user/Documents/GitHub/Planetary-Exploration-Tool/figs")

# # Visualize displacements
# deformation_map.visualize(projection=projection, time_point=1482224658.317582, direction="east")
#
# # Visualize displacements
# deformation_map.visualize(projection=projection, time_point=1482224658.317582, direction="north")
#
# # Visualize displacements
# deformation_map.visualize(projection=projection, time_point=1482224658.317582, direction="up")
#
# # Visualize displacements
# deformation_map.visualize(projection=projection, time_point=1482334678.317582, direction="east")
#
# # Visualize displacements
# deformation_map.visualize(projection=projection, time_point=1482334678.317582, direction="north")
#
# # Visualize displacements
# deformation_map.visualize(projection=projection, time_point=1482334678.317582, direction="up")

# # Visualize displacements
# displacements.visualize(projection=projection, time_point=1482219658.317582 + 110020, direction="east")
#
# # Visualize displacements
# displacements.visualize(projection=projection, time_point=1482219658.317582 + 110020, direction="north")
#
# # Visualize displacements
# displacements.visualize(projection=projection, time_point=1482219658.317582 + 110020, direction="up")

fm.clear()

# end of file
