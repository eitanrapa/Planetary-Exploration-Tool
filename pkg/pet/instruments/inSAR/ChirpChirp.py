#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import numpy as np
import pdb
import matplotlib.pyplot as plt
import sys
import pet
from scipy.interpolate import RegularGridInterpolator
import scipy.constants as const


class ChirpChirp(pet.component, family="pet.instruments.inSAR.chirpChirp", implements=pet.protocols.instruments.inSAR):
    """
    Defines the Nightingale ChirpChirp instrument
    """

    wavelength = pet.properties.float()
    wavelength.default = 0.13
    wavelength.doc = "radar instrument wavelength [m]"

    look_angle = pet.properties.float()
    look_angle.default = 23.5
    look_angle.doc = "look angle of the radar instrument"

    antenna_length = pet.properties.float()
    antenna_length.default = 2.
    antenna_length.doc = "antenna length (azimuth) [m]"

    antenna_height = pet.properties.float()
    antenna_height.default = 0.5
    antenna_height.doc = "antenna height (elevation) [m]"

    range_bandwidth = pet.properties.float()
    range_bandwidth.default = 30 * 1e6  # 30 MHz
    range_bandwidth.doc = "range bandwidth of the radar instrument [Hz]"

    antenna_efficiency = pet.properties.float()
    antenna_efficiency.default = -3
    antenna_efficiency.doc = "antenna efficiency [dB]"

    peak_power = pet.properties.float()
    peak_power.default = 200
    peak_power.doc = "peak power of the radar instrument [W]"

    pulse_duration = pet.properties.float()
    pulse_duration.default = 111 * 1e-6  # 60 microseconds
    pulse_duration.doc = "pulse duration of the radar instrument [s]"

    pulse_repetition_frequency = pet.properties.float()
    pulse_repetition_frequency.default = 450
    pulse_repetition_frequency.doc = "pulse repetition frequency of the radar instrument [Hz]"

    az_bandwidth_proc = pet.properties.float()
    az_bandwidth_proc.default = 64
    az_bandwidth_proc.doc = "processed azimuth bandwidth [Hz]"

    az_resolution_proc = pet.properties.float()
    az_resolution_proc.default = 1
    az_resolution_proc.doc = "processed azimuth resolution [m]"

    noise_figure = pet.properties.float()
    noise_figure.default = 3
    noise_figure.doc = "noise figure of the radar instrument [dB]"

    power_margin = pet.properties.float()
    power_margin.default = 5
    power_margin.doc = "power margin of the radar instrument [dB]"

    receiver_noise_temperature = pet.properties.float()
    receiver_noise_temperature.default = 240
    receiver_noise_temperature.doc = "receiver noise temperature of the radar instrument [K]"

    processing_azimuth_resolution = pet.properties.float()
    processing_azimuth_resolution.default = 20
    processing_azimuth_resolution.doc = "processing azimuth resolution of the radar instrument [m]"

    processing_ground_range_resolution = pet.properties.float()
    processing_ground_range_resolution.default = 20
    processing_ground_range_resolution.doc = "processing ground resolution of the radar instrument [m]"

    look_direction = pet.properties.str()
    look_direction.default = "right"
    look_direction.doc = "look direction of the radar instrument (right or left)"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Calculate the beamwidth
        d_inv = 1 / self.antenna_height
        self.bw = np.rad2deg(d_inv * self.wavelength)

        # Calculate the start and end look angles
        self.start_look_angle = self.look_angle - self.bw/2
        self.end_look_angle = self.look_angle + self.bw/2 

        return

    def get_instrument_noise(self, planet, baseline, satellite_velocity, look_angles, incidence_angles, distances):
        """
        Get instrument noise for the radar instrument
        :param planet: Planet object
        :param baseline: Baseline of the radar instrument [m]
        :param satellite_velocity: Satellite velocity [m/s]
        :param look_angles: Look angles of the radar instrument [rad]
        :param incidence_angles: Incidence angles of the radar instrument [rad]
        :param distances: Distances of the radar instrument [m]
        :return: Standard deviation of the phase of the interferogram [rad]
        """

        # Constants
        c0 = const.c  # speed of light
        kB = const.k  # Boltzman constant

        # Derived quantities
        SLC_res = [self.az_resolution_proc, c0 / (2 * self.range_bandwidth)]  # resolution of SAR image [azimuth, range] [m]
        # azBw_beam = float(2 * satellite_velocity / self.wavelength *
        #                   np.sin(self.wavelength / self.antenna_length))  # Doppler bandwidth
        ground_range_res = SLC_res[1] / np.sin(incidence_angles)  # ground range resolution

        # get 2-D antenna pattern and corresponding azimuth and elevation angles
        az_axis, el_axis, pattern2D = self.getAntennaPattern(wvlength=self.wavelength, lenAz=self.antenna_length,
                                                             lenEl=self.antenna_height,
                                                             efficiency=10**(self.antenna_efficiency/10))

        print("        Max. antenna gain is " + str(20*np.log10(np.max(pattern2D))) + "...", file=sys.stderr)
        losses_dB = 1.5 + self.power_margin  # losses [dB]
        losses = 10 ** (losses_dB / 10)

        vel = np.mean(np.linalg.norm(satellite_velocity, axis=-1))
        # compute the noise equivalent sigma naught
        NESN = self.computeNESN(incAng=incidence_angles, slant_range=distances, elevAng=look_angles-np.deg2rad(self.look_angle),
                                wvlength=self.wavelength, peakPower=self.peak_power, chirpDuration=self.pulse_duration,
                                PRF=self.pulse_repetition_frequency, rangeBw=self.range_bandwidth, doppBw=self.az_bandwidth_proc,
                                vSat=vel, noiseFigure=10**(self.noise_figure/10), losses=losses,
                                patTx2D=pattern2D, patRx2D=pattern2D, az_axis=az_axis, el_axis=el_axis, k_Boltz=kB,
                                noiseTemp=self.receiver_noise_temperature, Naux=256, dopplerCentroid=0,  c0=c0)
        print("        Mean NESN is " + str(10*np.log10(np.mean(NESN))) + "...", file=sys.stderr)

        # compute the SNR
        sigma0 = 3.51*np.cos(incidence_angles)**1.23#10**(planet.surface_backscatter/10)
        SNR = self.nesn2snr(nesn=NESN, sigma0=sigma0)
        print("        Mean SNR is " + str(10*np.log10(np.mean(SNR))) + "...", file=sys.stderr)

        # compute all relevant decorrelation sources
        corr_thermal = self.thermalCorrelation(SNR=SNR)
        corr_basline = self.baselineCorrelation(B=self.range_bandwidth, incAng=incidence_angles, baseline=baseline,
                                                wl=self.wavelength, r0=distances, c0=c0)
        corr_volume = self.volumeCorrelation(eps=planet.surface_permittivity, penDepth=planet.radar_penetration_depth,
                                             baseline=baseline, r0=distances, wl=self.wavelength,
                                             incAng=incidence_angles)

        # compute the total correlation
        corr_tot = corr_thermal * corr_basline * corr_volume
        print("        Mean thermal decorrelation is " + str(np.mean(corr_thermal)) + "...", file=sys.stderr)
        print("        Mean baseline decorrelation is " + str(np.mean(corr_basline)) + "...", file=sys.stderr)
        print("        Mean volume decorrelation is " + str(np.mean(corr_volume)) + "...", file=sys.stderr)
        print("        Mean total decorrelation is " + str(np.mean(corr_tot)) + "...", file=sys.stderr)

        
        # compute the available number of looks to average the interferogram
        Nlooks = ((self.processing_ground_range_resolution / ground_range_res) *
                  (self.processing_azimuth_resolution / SLC_res[0]))

        print("        Mean number of looks is " + str(np.mean(Nlooks)) + "...", file=sys.stderr)

        # compute the standard deviation of the phase of the interferogram
        sigma_phase = self.coherence2phasenoise(coherence=corr_tot, effectiveLooks=Nlooks)

        print("        Mean phase standard deviation is " + str(np.mean(sigma_phase)) + "...", file=sys.stderr)

        return sigma_phase, corr_tot, Nlooks, NESN, sigma0

    def thermalCorrelation(self, SNR):
        thermal_corr = (1 / (1 + 1 / SNR))
        return thermal_corr

    def volumeCorrelation(self, eps, penDepth, baseline, r0, wl, incAng):
        refAng = np.arcsin(np.sin(incAng)/np.sqrt(eps))
        HoA = wl * r0 * np.sin(incAng) / (2*baseline)
        HoA_vol = HoA/np.sqrt(eps) * np.cos(refAng) / np.cos(incAng)
        dp2 = penDepth/2 * np.cos(refAng)

        corr = abs(1/(1 + 1j * 2 * np.pi * dp2/HoA_vol))
        # corr = 1 / np.sqrt(1 + (2 * np.pi * np.sqrt(eps) * penDepth * baseline / (r0 * wl * np.tan(incAng))) ** 2)
        return corr

    def criticalBaseline(self, wl, r0, incAng, B, c0):
        critical_baseline = B * wl * r0 * np.tan(incAng) / c0
        return critical_baseline

    def baselineCorrelation(self, B, incAng, baseline, wl, r0, c0):
        baseCrit = self.criticalBaseline(wl=wl, r0=r0, incAng=incAng, B=B, c0=c0)
        corr = 1 - baseline / baseCrit
        return corr

    def nesn2snr(self, nesn, sigma0):
        nesn2snr = sigma0 / nesn
        return nesn2snr

    def coherence2phasenoise(self, coherence, effectiveLooks):
        """ Estimates the Standard Deviation of the Interferometric Phase Error from coherence, based on Cramer-Rao  """
        phasenoise = np.sqrt((1. - coherence ** 2.) / (2. * effectiveLooks * coherence ** 2.))
        return phasenoise

    def getOutputDoppler(self, prf, Na):
        doppler = -prf / 2 + np.arange(Na) / Na * prf
        return doppler

    def getPowerInBand(self, pat, Nabw, idxAz):
        powAbs = abs(pat)
        powInBand = np.roll((np.cumsum(powAbs) - np.roll(np.cumsum(powAbs), Nabw)),
                            -int(0.5 * Nabw))
        return powInBand[idxAz]

    def Regular2DInterp(self, data, naInt, nrInt):
        Na = len(data[:, 0])
        Nr = len(data[0, :])
        fint = RegularGridInterpolator((np.arange(Na), np.arange(Nr)), data,
                                       bounds_error=False, fill_value=0)
        return fint((naInt, nrInt))

    def evalPattern(self, patIn, azAng, elAng, azOut, elOut):
        try:
            # assumption is elevation is the fast dimension!!!
            idxAz = (azOut - azAng[0]) / (azAng[1] - azAng[0])
            idxEl = (elOut - elAng[0]) / (elAng[1] - elAng[0])
            patOut = self.Regular2DInterp(data=patIn, naInt=idxAz, nrInt=idxEl)
        except:
            idxAz = (azOut - azAng[0]) / (azAng[1] - azAng[0])
            idxEl = (elOut - elAng[0]) / (elAng[1] - elAng[0])
            patOut = self.Regular2DInterp(data=patIn, naInt=idxAz, nrInt=idxEl)
        return patOut

    def computeNESN(self, incAng, slant_range, elevAng, wvlength, peakPower, chirpDuration,
                    PRF, rangeBw, doppBw, vSat, noiseFigure, losses, patTx2D, patRx2D,
                    az_axis, el_axis, k_Boltz, noiseTemp, Naux, dopplerCentroid, c0):
        """ Computes NESN of a single channel system
            :param patTx2D: 2D antenna pattern with [az, el]
        """

        faOut = self.getOutputDoppler(prf=PRF, Na=Naux)
        Nabw = int(np.round(doppBw / (faOut[1] - faOut[0])))  # number of bins corresponding to Doppler bandwidth
        idxAz = np.round((dopplerCentroid - faOut[0]) / (faOut[1] - faOut[0])).astype(int)  # index of doppler centroid
        thOut = np.arcsin(wvlength * faOut / 2 / vSat)
        angAp = abs(thOut[-1] - thOut[
            0])  # angular extent of the PRF
        powerDensity = noiseFigure * k_Boltz * noiseTemp
        noisePow = powerDensity * rangeBw

        Nr_out = len(slant_range)
        nesn = np.zeros([Nr_out], np.float64)
        for ii in range(Nr_out):
            patSignal = abs(self.evalPattern(patIn=patTx2D, azAng=az_axis, elAng=el_axis, azOut=thOut,
                                             elOut=elevAng[ii]) *
                            self.evalPattern(patIn=patRx2D, azAng=az_axis, elAng=el_axis, azOut=thOut,
                                             elOut=elevAng[ii]))
            # Signal power at antenna surface
            signalPow = peakPower * wvlength ** 2 * c0 * chirpDuration * angAp / 2 / (
                    4 * np.pi * slant_range[ii]) ** 3 / losses / np.sin(incAng[ii])
            # Power spectral densities [dbW]
            psdSignal = signalPow * (patSignal) ** 2 / PRF
            psdNoise = noisePow * (np.ones(Naux)) ** 2 / PRF
            # Integrate power over processed band (Nabw)
            powIntSignal = self.getPowerInBand(pat=psdSignal, Nabw=Nabw, idxAz=idxAz)
            powIntNoise = self.getPowerInBand(pat=psdNoise, Nabw=Nabw, idxAz=idxAz)
            nesn[ii] = powIntNoise / powIntSignal
        return nesn

    def computeAntennaGain(self, length, height, wvlength, efficiency):
        factor = 4 * np.pi / wvlength ** 2 * efficiency
        gain = factor * length * height
        return gain

    def getElementPattern(self, thOut, wvlength, L):
        th3dB = wvlength / L
        return np.sinc(np.sin(thOut) / th3dB) + 0j

    def get2DPattern(self, azAng, elAng, wvlength, lenAz, lenEl, efficiency):
        antennaGain = self.computeAntennaGain(length=lenAz, height=lenEl, wvlength=wvlength, efficiency=efficiency)
        pattern2D = (self.getElementPattern(thOut=azAng, wvlength=wvlength, L=lenAz).reshape(len(azAng), 1) *
                     self.getElementPattern(thOut=elAng, wvlength=wvlength, L=lenEl).reshape(1, len(elAng)))
        return pattern2D / np.max(abs(pattern2D)) * np.sqrt(antennaGain)

    def getAntennaPattern(self, wvlength, lenAz, lenEl, efficiency):
        elExt = 3 * wvlength / lenEl
        azExt = 3 * wvlength / lenAz
        az_axis = np.linspace(-azExt / 2, azExt / 2, 501)
        el_axis = np.linspace(-elExt / 2, elExt / 2, 501)
        pattern = abs(self.get2DPattern(azAng=az_axis, elAng=el_axis, wvlength=wvlength, lenAz=lenAz,
                                        lenEl=lenEl, efficiency=efficiency))
        return az_axis, el_axis, pattern

# end of file
