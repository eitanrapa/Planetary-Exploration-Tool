#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from scipy.interpolate import griddata
import cspyce as spice
import pdb
import os
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

# # Furnish some files
# fm.furnsh(names_list=["cas_enceladus_ssd_spc_1024icq_v1.bds", "pck00011_n0066.tpc",
#                       "insar_6stride_26d_v7_seo.bsp", "latest_leapseconds.tls"])

# # Furnish some files
# fm.furnsh(names_list=["cas_enceladus_ssd_spc_1024icq_v1.bds", "pck00011_n0066.tpc",
#                       "orb61refv5_foundation_eastward_6mo_scpse.bsp", "latest_leapseconds.tls"])

# Furnish some files
fm.furnsh(names_list=["cas_enceladus_ssd_spc_1024icq_v1.bds", "pck00011_n0066.tpc",
                      "insar_5stride_130d_v1_scpse (1).bsp", "latest_leapseconds.tls"])


# Make a planet
planet = pet.planets.enceladus(name="enceladus")

# Make an instrument
instrument = pet.instruments.inSAR.chirpChirp(name="chirp chirp")

# Make a con ops
# campaign = pet.campaigns.orbiter.nightingale5to1(name="nightingale",
#                                                  body_id=-303, start_time="2046 DEC 20 15:10:40.134", planet=planet)
campaign = pet.campaigns.orbiter.nightingaleFoundational(name="nightingale",
                                                 body_id=-303, start_time="2046 DEC 20 15:10:40.134", planet=planet)

#for foundational (start time of kernel is Start Time: 2046 JUN 14 01:33:22)
campaign = pet.campaigns.orbiter.nightingaleFoundational(name="nightingale",
                                                 body_id=-303, start_time="2046 JUN 15 00:00:00.001", planet=planet)

# Get the start time
startTime = campaign.start_time_et

# # Get the orbit cycle time of the instrument
# orbit_cycle_time = campaign.orbit_cycle

# Get the track
start = startTime + 3600*6 - 3600*2
end = start + 3600*6 * 2
track = pet.dataAcquisition.track(name="track", start_time=start, end_time=end, planet=planet,
                                  campaign=campaign, instrument=instrument, spatial_resolution=200,
                                  temporal_resolution=7, interpol2FinalRes=0, interferogram_resolution_az=400,
                                  interferogram_resolution_rg=400)
track.calculate_ground_swath_new()

#100 and 5

# Specify the baseline
baseline = 52
# Specify the basline uncertainty
baseline_uncertainty = 0.
# Specify the baseline orientation (roll=0 means horizontal)
roll = np.deg2rad(23.5)
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
heightError_corr = planet.get_height_above_surface(xyz_meas)
#---------------------------------------------

kz = interferogram.data['kz'].kz.values
HoA = 2 * np.pi / kz
incAng = track.data.values
phase_absolut = interferogram.data['phase_absolut'].phase_absolut.values
phase_ref = interferogram.data['phase_ref'].phase_ref.values
phase_flat = phase_absolut - phase_ref
phase_flat_wrapped = phase_flat%(2*np.pi)
height_true = track.data["height"].values


def plotData(projection, longitudes, latitudes, data, cmap, s, marker, vmin, vmax, extent,filename,title, isSouth=False):
    if isSouth:
        os.environ["PROJ_IGNORE_CELESTIAL_BODY"] = "YES"
        custom_globe = ccrs.Globe(semimajor_axis=256600.,
                            semiminor_axis=251400.,
                            ellipse=None)
        projectionSouth = ccrs.SouthPolarStereo(globe=custom_globe)
        fig = plt.figure('south')
        ax = plt.axes(projection=projectionSouth)
        # Add features
        gl = ax.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
        gl.xlocator = plt.MultipleLocator(30)
        gl.ylocator = plt.MultipleLocator(10)
        im = ax.scatter(longitudes, latitudes, transform=ccrs.PlateCarree(), c=data,s=s, cmap=cmap, marker=marker,vmin=vmin,vmax=vmax)
        # Set extent to show the Antarctic region
        ax.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
        plt.colorbar(im, fraction=0.022, pad=0.1)
        plt.title(title, pad=12)
        plt.savefig(filename, format='png',dpi=300,bbox_inches="tight")
        fig.clf()
    else:
        os.environ["PROJ_IGNORE_CELESTIAL_BODY"] = "YES"
        custom_globe = ccrs.Globe(semimajor_axis=256600.,
                            semiminor_axis=251400.,
                            ellipse=None)
        projectionGlob = ccrs.PlateCarree(globe=custom_globe)
        fig = plt.figure('global')
        ax = plt.axes(projection=projectionGlob)
        # Add features
        gl = ax.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
        gl.xlocator = plt.MultipleLocator(30)
        gl.ylocator = plt.MultipleLocator(10)
        im = ax.scatter(longitudes, latitudes, transform=ccrs.PlateCarree(), c=data,s=s, cmap=cmap, marker=marker,vmin=vmin,vmax=vmax)
        # Set extent to show the Antarctic region
        ax.set_extent(extent, ccrs.PlateCarree())
        plt.colorbar(im, fraction=0.022, pad=0.1)
        plt.title(title, pad=12)
        plt.savefig(filename, format='png',dpi=300,bbox_inches="tight")
        fig.clf()





#plot height-------------------------------------------------------------------
# extent = [np.min(longitudes), np.max(longitudes), np.min(latitudes), np.max(latitudes)]
extent = [-180, 180, -90, 90]
s=0.2
marker='.'
isSouth=True
projection = pet.projections.biaxialProjections.biaxialCylindrical(name="biaxial cylindrical",
                                                                        folder_path="/data/largeHome/bene_an/Projects/nightingale/Simulation/Planetary-Exploration-Tool/output/plots/")

plotTime = (track.data['time'].values - np.min(track.data['time'].values)) / 3600 / 24
plotData(projection, longitudes, latitudes, plotTime, cmap='jet', s=s,
marker=marker, vmin=np.min(plotTime), vmax=np.max(plotTime),
extent=extent,filename='output/plots/time.png',
title='time [day]', isSouth=isSouth)

plotData(projection, longitudes, latitudes, heightError_corr, cmap='jet', s=s,
marker=marker, vmin=np.mean(heightError_corr) - np.std(heightError_corr)*3,
vmax=np.mean(heightError_corr) + np.std(heightError_corr)*3,
extent=extent,filename='output/plots/heightError.png',
title='height error [m]', isSouth=isSouth)

plotData(projection, longitudes, latitudes, heights, cmap='terrain', s=s,
marker=marker, vmin=np.min(heights), vmax=np.max(heights),
extent=extent,filename='output/plots/height_meas.png',
title='measured height [m]', isSouth=isSouth)

plotData(projection, longitudes, latitudes, 10*np.log10(interferogram.NESN), cmap='Spectral', s=s,
marker=marker, vmin=np.min(10*np.log10(interferogram.NESN))+3, vmax=np.max(10*np.log10(interferogram.NESN))-2,
extent=extent,filename='output/plots/NESN.png',
title='NESN [dB]', isSouth=isSouth)

plotData(projection, longitudes, latitudes, interferogram.corr_tot, cmap='gray', s=s,
marker=marker, vmin=0., vmax=0.5,
extent=extent,filename='output/plots/coherence.png',
title='coherence []', isSouth=isSouth)

sigma0Plot = 10*np.log10(interferogram.sigma0)
plotData(projection, longitudes, latitudes, sigma0Plot, cmap='gray', s=s,
marker=marker, vmin=np.min(sigma0Plot)+2+4+5, vmax=np.max(sigma0Plot),
extent=extent,filename='output/plots/sigma0.png',
title='backscatter [dB]', isSouth=isSouth)

plotData(projection, longitudes, latitudes, np.rad2deg(incAng), cmap='Spectral', s=s,
marker=marker, vmin=np.min(np.rad2deg(incAng)), vmax=np.max(np.rad2deg(incAng)),
extent=extent,filename='output/plots/incAng.png',
title='incidence angle [deg]', isSouth=isSouth)

plotData(projection, longitudes, latitudes, phase_flat_wrapped, cmap='hsv', s=s,
marker=marker, vmin=np.min(phase_flat_wrapped), vmax=np.max(phase_flat_wrapped),
extent=extent,filename='output/plots/phase_wrapped.png',
title='wrapped phase [rad]', isSouth=isSouth)

plotData(projection, longitudes, latitudes, interferogram.sigma_height, cmap='hsv', s=s,
marker=marker, vmin=1, vmax=11,
extent=extent,filename='output/plots/sigma_height.png',
title='height standard deviation (model) [m]', isSouth=isSouth)

plotData(projection, longitudes, latitudes, interferogram.sigma_disp*1e3, cmap='hsv', s=s,
marker=marker, vmin=np.min(interferogram.sigma_disp*1e3), vmax=np.max(interferogram.sigma_disp*1e3)-14-25,
extent=extent,filename='output/plots/sigma_disp.png',
title='displacement standard deviation (model) [mm]', isSouth=isSouth)
# #------------------------------------------------------------------------------

# campaign.plot_orbit(projection=projection, start_time=campaign.start_time,
#                     end_time="2046 JUN 29 00:00:00.001", temporal_resolution=10)

# a, b, c = planet.get_axes()




#plot some interferometric products--------------------------------------------
# interferogram.visualize_interferogram(phase_flat, projection, "hsv","flattened phase [rad]", "output/plots/phase_flat.png")
# interferogram.visualize_interferogram(phase_absolut, projection, "hsv", "absolut phase [rad]", "output/plots/phase_absolut.png")
# interferogram.visualize_interferogram(phase_ref, projection, "hsv","reference phase [rad]", "output/plots/phase_reference.png")
# interferogram.visualize_interferogram(phase_flat_wrapped, projection, "hsv","wrapped phase [rad]", "output/plots/phase_wrapped.png")
# interferogram.visualize_interferogram(np.rad2deg(incAng), projection, "Spectral","incident angle knowledge [deg]", "output/plots/incAng_know.png")
# interferogram.visualize_interferogram(np.rad2deg(track.data.values), projection, "Spectral", "incident angle [deg]", "output/plots/incAng.png")
# interferogram.visualize_interferogram(HoA, projection, "Spectral","height of ambiguity knowledge [m]", "output/plots/HoA_know.png")
# interferogram.visualize_interferogram(height_true, projection, "terrain","DEM [m]", "output/plots/DEM.png")
#------------------------------------------------------------------------------

fm.clear()

# end of file
