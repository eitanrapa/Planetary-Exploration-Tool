#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pet


def convert_positions(cartesian_positions, planet):
    """
    Convert flattened points to geodetic coordinates
    :param cartesian_positions: Points to convert
    :return: Points in geodetic coordinates
    """
    # Get planet axes
    a, b, c = planet.get_axes()

    # Create a coordinate conversion object
    coordinate_conversions = pet.conversions.coordinateConversions(name="conversions", a=a, b=b, c=c)

    # Get the corresponding latitude and longitude with a triaxial ellipsoid conversion
    flattened_swath_coordinates = np.asarray(coordinate_conversions.geodetic(
        cartesian_coordinates=cartesian_positions))

    return flattened_swath_coordinates





# Create a file manager
fm = pet.spiceTools.fileManager(folder_path="/data/largeHome/bene_an/Projects/nightingale/Simulation/Planetary-Exploration-Tool/data/")

# Furnish some files
fm.furnsh(names_list=["cas_enceladus_ssd_spc_1024icq_v1.bds", "pck00011_n0066.tpc",
                      "insar_6stride_26d_v7_seo.bsp", "latest_leapseconds.tls"])

# Make a planet
planet = pet.planets.enceladus(name="enceladus")

# Make an instrument
instrument = pet.instruments.inSAR.chirpChirp(name="chirp chirp")

# Make a con ops
campaign = pet.campaigns.orbiter.nightingale5to1(name="nightingale",
                                                 body_id=-303, start_time="2046 DEC 20 15:10:40.134", planet=planet)

# Get the times defining the first five tracks
times = campaign.get_five_tracks()

# Get the orbit cycle time of the instrument
orbit_cycle_time = campaign.orbit_cycle

# Get the track
track = pet.dataAcquisition.track(name="track", start_time=times[0], end_time=times[0] + 60*90, planet=planet,
                                  campaign=campaign, instrument=instrument, spatial_resolution=700,
                                  temporal_resolution=20, interpol2FinalRes=1, interferogram_resolution_az=700,
                                  interferogram_resolution_rg=700)
track.calculate_ground_swath_new()



# Specify the baseline
baseline = 200.
# Specify the basline uncertainty
baseline_uncertainty = 0.
# Specify the baseline orientation (roll=0 means horizontal)
roll = np.deg2rad(0)
# Specify the roll uncertainty
roll_uncertainty = 0.#np.deg2rad(-8 / 3600)
# compute the interferogram
interferogram = pet.dataAnalysis.surfaceDeformationAnalysis.simpleInterferogramTopo(name="igram",
                                                                                instrument=instrument,
                                                                                planet=planet,
                                                                                track=track,
                                                                                campaign=campaign,
                                                                                baseline=baseline,
                                                                                baseline_uncertainty=baseline_uncertainty,
                                                                                roll=roll,
                                                                                roll_uncertainty=roll_uncertainty)
# Calculate interferogram
interferogram.calculate_igram()
# Save interferogram
interferogram.save(file_name="/data/largeHome/bene_an/Projects/nightingale/Simulation/Planetary-Exploration-Tool/output/files/iinterferogram")

# read interferometric results----------------
xyz_meas = interferogram.xyz_meas
llh = convert_positions(xyz_meas,planet)
longitudes = llh[:,1]
latitudes = llh[:,0]
heights = llh[:,2]
#get height error
heightError_corr = -1*planet.get_height_above_surface(xyz_meas)
#---------------------------------------------

#plot height-------------------------------------------------------------------
projection = pet.projections.biaxialProjections.biaxialCylindrical(name="biaxial cylindrical",
                                                                        folder_path="/data/largeHome/bene_an/Projects/nightingale/Simulation/Planetary-Exploration-Tool/output/plots/")
fig1, ax1, globe1 = projection.proj(planet=planet)
# Make the colormap cyclical
cm = plt.cm.get_cmap('terrain')
im = ax1.scatter(longitudes, latitudes, cmap=cm,
                transform=ccrs.PlateCarree(globe=globe1), c=heightError_corr, s=0.01,marker=',')
plt.colorbar(im, fraction=0.02, pad=0.1)
plt.title('height error [m]', pad=12)
plt.savefig('heightError.png', format='png',dpi=300,bbox_inches="tight")

fig2, ax2, globe2 = projection.proj(planet=planet)
# Make the colormap cyclical
cm = plt.cm.get_cmap('terrain')
im = ax2.scatter(longitudes, latitudes, cmap=cm,
                transform=ccrs.PlateCarree(globe=globe2), c=heights, s=0.01,marker=',')
plt.colorbar(im, fraction=0.02, pad=0.1)
plt.title('measured height [m]', pad=12)
plt.savefig('height_meas.png', format='png',dpi=300,bbox_inches="tight")

fig3, ax3, globe3 = projection.proj(planet=planet)
# Make the colormap cyclical
cm = plt.cm.get_cmap('Spectral')
im = ax3.scatter(longitudes, latitudes, cmap=cm,
                transform=ccrs.PlateCarree(globe=globe3), c=10*np.log10(interferogram.NESN), s=0.01,marker=',')
plt.colorbar(im, fraction=0.02, pad=0.1)
plt.title('NESN [dB]', pad=12)
plt.savefig('NESN.png', format='png',dpi=300,bbox_inches="tight")

fig4, ax4, globe4 = projection.proj(planet=planet)
# Make the colormap cyclical
cm = plt.cm.get_cmap('Spectral')
im = ax4.scatter(longitudes, latitudes, cmap=cm,
                transform=ccrs.PlateCarree(globe=globe4), c=interferogram.corr_tot, s=0.01,marker=',')
plt.colorbar(im, fraction=0.02, pad=0.1)
plt.title('total decorrelation []', pad=12)
plt.savefig('coherence.png', format='png',dpi=300,bbox_inches="tight")

fig5, ax5, globe5 = projection.proj(planet=planet)
# Make the colormap cyclical
cm = plt.cm.get_cmap('Spectral')
im = ax5.scatter(longitudes, latitudes, cmap=cm,
                transform=ccrs.PlateCarree(globe=globe5), c=interferogram.sigma_phase, s=0.01,marker=',')
plt.colorbar(im, fraction=0.02, pad=0.1)
plt.title('phase standard deviation [rad]', pad=12)
plt.savefig('phase_std.png', format='png',dpi=300,bbox_inches="tight")
#------------------------------------------------------------------------------


#plot some interferometric products--------------------------------------------
kz = interferogram.data['kz'].kz.values
HoA = 2 * np.pi / kz
incAng = interferogram.data['incAng'].incAng.values
phase_absolut = interferogram.data['phase_absolut'].phase_absolut.values
phase_ref = interferogram.data['phase_ref'].phase_ref.values
phase_flat = phase_absolut - phase_ref
phase_flat_wrapped = phase_flat%(2*np.pi)
height_true = track.data["height"].values

interferogram.visualize_interferogram(phase_flat, projection, "hsv","flattened phase [rad]", "phase_flat.png")
interferogram.visualize_interferogram(phase_absolut, projection, "hsv", "absolut phase [rad]", "phase_absolut.png")
interferogram.visualize_interferogram(phase_ref, projection, "hsv","reference phase [rad]", "phase_reference.png")
interferogram.visualize_interferogram(phase_flat_wrapped, projection, "hsv","wrapped phase [rad]", "phase_wrapped.png")
interferogram.visualize_interferogram(np.rad2deg(incAng), projection, "Spectral","incident angle knowledge [deg]", "incAng_know.png")
interferogram.visualize_interferogram(np.rad2deg(track.data.values), projection, "Spectral", "incident angle [deg]", "incAng.png")
interferogram.visualize_interferogram(HoA, projection, "Spectral","height of ambiguity knowledge [m]", "HoA_know.png")
interferogram.visualize_interferogram(height_true, projection, "terrain","DEM [m]", "DEM.png")
#------------------------------------------------------------------------------

fm.clear()

# end of file
