#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import numpy as np
import sys
import pet
from scipy.interpolate import RegularGridInterpolator
import scipy.constants as const

# Constants
C = const.c  # speed of light
K_B = const.k  # Boltzman constant


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
    antenna_length.default = 2.0
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
        self.start_look_angle = self.look_angle - self.bw / 2
        self.end_look_angle = self.look_angle + self.bw / 2

        return

    @pet.export
    def get_instrument_noise(self, planet, baseline, satellite_velocities, look_angles, incidence_angles, distances,
                             variable_backscatter=False):
        """
        Helper function to get the instrument noise
        """

        # Make sure the inputs are numpy arrays
        satellite_velocities = np.asarray(satellite_velocities)
        look_angles = np.asarray(look_angles)
        incidence_angles = np.asarray(incidence_angles)
        distances = np.asarray(distances)

        # Make look angles and incidence angles radians
        look_angles = np.deg2rad(look_angles)
        incidence_angles = np.deg2rad(incidence_angles)

        return self._get_instrument_noise(planet, baseline, satellite_velocities, look_angles,
                                          incidence_angles, distances, variable_backscatter)

    def _get_instrument_noise(self, planet, baseline, satellite_velocities, look_angles, incidence_angles, distances,
                             variable_backscatter=False):
        """
        Get instrument noise for the radar instrument
        :param planet: Planet object
        :param baseline: Perpendicular baseline of the radar instrument [m]
        :param satellite_velocities: Satellite velocities [m/s]
        :param look_angles: Look angles of the radar instrument [rad]
        :param incidence_angles: Incidence angles of the radar instrument [rad]
        :param distances: Distances of the radar instrument [m]
        :param variable_backscatter: Flag for variable backscatter
        :return: Standard deviation of the phase of the interferogram [rad]
        """

        # Derived quantities
        slc_res = [self.az_resolution_proc,
                   C / (2 * self.range_bandwidth)]  # resolution of SAR image [azimuth, range] [m]
        ground_range_res = slc_res[1] / np.sin(incidence_angles)  # ground range resolution

        # Get satellite speeds
        satellite_speeds = np.linalg.norm(satellite_velocities, axis=-1)

        # get 2-D antenna pattern and corresponding azimuth and elevation angles
        az_axis, el_axis, pattern_2d = self.get_antenna_pattern(wavelength=self.wavelength, len_az=self.antenna_length,
                                                                len_el=self.antenna_height,
                                                                efficiency=10 ** (self.antenna_efficiency / 10))

        print("        Max. antenna gain is " + str(20 * np.log10(np.max(pattern_2d))) + "...", file=sys.stderr)

        # Calculate losses
        losses_db = 1.5 + self.power_margin  # losses [dB]
        losses = 10 ** (losses_db / 10)

        # compute the noise equivalent sigma naught, take the mean of the speeds
        nesn = self.compute_nesn(incidence_angle=incidence_angles, slant_range=distances,
                                elevation_angle=look_angles - np.deg2rad(self.look_angle),
                                wavelength=self.wavelength, peak_power=self.peak_power,
                                chirp_duration=self.pulse_duration,
                                prf=self.pulse_repetition_frequency, range_bw=self.range_bandwidth,
                                dopp_bw=self.az_bandwidth_proc,
                                v_sat=np.mean(satellite_speeds), noise_figure=10 ** (self.noise_figure / 10),
                                losses=losses,
                                pat_tx_2d=pattern_2d, pat_rx_2d=pattern_2d, az_axis=az_axis, el_axis=el_axis,
                                noise_temp=self.receiver_noise_temperature, n_aux=256, doppler_centroid=0)

        print("        Mean NESN is " + str(10 * np.log10(np.mean(nesn))) + "...", file=sys.stderr)

        # compute the SNR
        if variable_backscatter:
            sigma0 = 3.51 * np.cos(incidence_angles) ** 1.23
        else:
            sigma0 = 10**(planet.surface_backscatter/10)

        snr = self.nesn2snr(nesn=nesn, sigma0=sigma0)

        print("        Mean SNR is " + str(10 * np.log10(np.mean(snr))) + "...", file=sys.stderr)

        # compute all relevant decorrelation sources
        corr_thermal = self.thermal_correlation(snr=snr)
        corr_baseline = self.baseline_correlation(bandwidth=self.range_bandwidth,
                                                  incidence_angle=incidence_angles, baseline=baseline,
                                                  wavelength=self.wavelength, slant_range=distances)
        corr_volume = self.volume_correlation(eps=planet.surface_permittivity,
                                              penetration_depth=planet.radar_penetration_depth,
                                              baseline=baseline, slant_range=distances, wavelength=self.wavelength,
                                              incidence_angle=incidence_angles)

        # compute the total correlation
        corr_tot = corr_thermal * corr_baseline * corr_volume
        print("        Mean thermal decorrelation is " + str(np.mean(corr_thermal)) + "...", file=sys.stderr)
        print("        Mean baseline decorrelation is " + str(np.mean(corr_baseline)) + "...", file=sys.stderr)
        print("        Mean volume decorrelation is " + str(np.mean(corr_volume)) + "...", file=sys.stderr)
        print("        Mean total decorrelation is " + str(np.mean(corr_tot)) + "...", file=sys.stderr)

        # compute the available number of looks to average the interferogram
        n_looks = ((self.processing_ground_range_resolution / ground_range_res) *
                  (self.processing_azimuth_resolution / slc_res[0]))

        print("        Mean number of looks is " + str(np.mean(n_looks)) + "...", file=sys.stderr)

        # compute the standard deviation of the phase of the interferogram
        sigma_phase = self.coherence2phasenoise(coherence=corr_tot, effective_looks=n_looks)

        print("        Mean phase standard deviation is " + str(np.mean(sigma_phase)) + "...", file=sys.stderr)

        # Return sigma phase, correlation, number of looks and NESN
        return sigma_phase, corr_tot, n_looks, nesn, sigma0

    def get_2d_pattern(self, az_ang, elevation_angle, wavelength, len_az, len_el, efficiency):
        """
        Computes the 2D antenna pattern
        :param az_ang: azimuth angles [rad]
        :param elevation_angle: elevation angles [rad]
        :param wavelength: wavelength [m]
        :param len_az: length of the antenna in azimuth [m]
        :param len_el: length of the antenna in elevation [m]
        :param efficiency: antenna efficiency [dB]
        :return: 2D antenna pattern
        """

        antenna_gain = self.compute_antenna_gain(length=len_az, height=len_el, wavelength=wavelength,
                                                efficiency=efficiency)
        pattern_2d = (self.get_element_pattern(th_out=az_ang, wavelength=wavelength,
                                               length=len_az).reshape(len(az_ang), 1) *
                      self.get_element_pattern(th_out=elevation_angle, wavelength=wavelength,
                                               length=len_el).reshape(1, len(elevation_angle)))
        return pattern_2d / np.max(abs(pattern_2d)) * np.sqrt(antenna_gain)

    def get_antenna_pattern(self, wavelength, len_az, len_el, efficiency):
        """
        Get the 2D antenna pattern
        :param wavelength: wavelength [m]
        :param len_az: length of the antenna in azimuth [m]
        :param len_el: length of the antenna in elevation [m]
        :param efficiency: antenna efficiency [dB]
        :return: azimuth axis, elevation axis, 2D antenna pattern
        """

        el_ext = 3 * wavelength / len_el
        az_ext = 3 * wavelength / len_az
        az_axis = np.linspace(-az_ext / 2, az_ext / 2, 501)
        el_axis = np.linspace(-el_ext / 2, el_ext / 2, 501)
        pattern = abs(self.get_2d_pattern(az_ang=az_axis, elevation_angle=el_axis, wavelength=wavelength, len_az=len_az,
                                        len_el=len_el, efficiency=efficiency))
        return az_axis, el_axis, pattern

    def eval_pattern(self, pat_in, az_ang, elevation_ang, az_out, el_out):
        """
        Evaluate the antenna pattern
        :param pat_in: 2D antenna pattern
        :param az_ang: azimuth angles [rad]
        :param elevation_ang: elevation angles [rad]
        :param az_out: output azimuth angles [rad]
        :param el_out: output elevation angles [rad]
        """

        # assumption is elevation is the fast dimension!!!
        idx_az = (az_out - az_ang[0]) / (az_ang[1] - az_ang[0])
        idx_el = (el_out - elevation_ang[0]) / (elevation_ang[1] - elevation_ang[0])
        pat_out = self.regular_2d_interp(data=pat_in, na_int=idx_az, nr_int=idx_el)

        return pat_out

    def compute_nesn(self, incidence_angle, slant_range, elevation_angle, wavelength, peak_power, chirp_duration,
                    prf, range_bw, dopp_bw, v_sat, noise_figure, losses, pat_tx_2d, pat_rx_2d,
                    az_axis, el_axis, noise_temp, n_aux, doppler_centroid):
        """
        Computes NESN of a single channel system
        :param incidence_angle: incidence angle [rad]
        :param slant_range: slant range [m]
        :param elevation_angle: elevation angle [rad]
        :param wavelength: wavelength [m]
        :param peak_power: peak power [W]
        :param chirp_duration: chirp duration [s]
        :param prf: pulse repetition frequency [Hz]
        :param range_bw: range bandwidth [Hz]
        :param dopp_bw: Doppler bandwidth [Hz]
        :param v_sat: satellite velocity [m/s]
        :param noise_figure: noise figure [dB]
        :param losses: losses [dB]
        :param pat_tx_2d: 2D antenna pattern with [az, el]
        :param pat_rx_2d: 2D antenna pattern with [az, el]
        :param az_axis: azimuth axis [rad]
        :param el_axis: elevation axis [rad]
        :param noise_temp: noise temperature [K]
        :param n_aux: number of auxiliary samples
        :param doppler_centroid: Doppler centroid [Hz]
        :return: NESN [dB]
        """

        fa_out = self.get_output_doppler(prf=prf, na=n_aux)
        na_bw = int(np.round(dopp_bw / (fa_out[1] - fa_out[0])))  # number of bins corresponding to Doppler bandwidth
        idx_az = np.round((doppler_centroid - fa_out[0]) / (fa_out[1] - fa_out[0])).astype(int)  # index of doppler centroid

        th_out = np.arcsin(wavelength * fa_out / 2 / v_sat)
        ang_ap = abs(th_out[-1] - th_out[0])  # angular extent of the PRF
        power_density = noise_figure * K_B * noise_temp
        noise_pow = power_density * range_bw

        nr_out = len(slant_range)
        nesn = np.zeros([nr_out], np.float64)
        for ii in range(nr_out):
            pat_signal = abs(self.eval_pattern(pat_in=pat_tx_2d, az_ang=az_axis, elevation_ang=el_axis, az_out=th_out,
                                               el_out=elevation_angle[ii]) *
                             self.eval_pattern(pat_in=pat_rx_2d, az_ang=az_axis, elevation_ang=el_axis, az_out=th_out,
                                               el_out=elevation_angle[ii]))
            # Signal power at antenna surface
            signal_pow = peak_power * wavelength ** 2 * C * chirp_duration * ang_ap / 2 / (
                    4 * np.pi * slant_range[ii]) ** 3 / losses / np.sin(incidence_angle[ii])
            # Power spectral densities [dbW]
            psd_signal = signal_pow * pat_signal ** 2 / prf
            psd_noise = noise_pow * (np.ones(n_aux)) ** 2 / prf
            # Integrate power over processed band (Nabw)
            pow_int_signal = self.get_power_in_band(pat=psd_signal, na_bw=na_bw, idx_az=idx_az)
            pow_int_noise = self.get_power_in_band(pat=psd_noise, na_bw=na_bw, idx_az=idx_az)
            nesn[ii] = pow_int_noise / pow_int_signal
        return nesn

    def baseline_correlation(self, bandwidth, incidence_angle, baseline, wavelength, slant_range):
        """
        Computes the baseline decorrelation based on the critical baseline
        :param bandwidth: range bandwidth [Hz]
        :param incidence_angle: incidence angle [rad]
        :param baseline: baseline [m]
        :param wavelength: wavelength [m]
        :param slant_range: slant range [m]
        :return: Baseline decorrelation
        """

        base_crit = self.critical_baseline(wavelength=wavelength, slant_range=slant_range,
                                           incidence_angle=incidence_angle, bandwidth=bandwidth)
        corr = 1 - baseline / base_crit
        return corr

    @staticmethod
    def thermal_correlation(snr):
        """
        Computes the thermal decorrelation based on the SNR
        :param snr: Signal to noise ratio
        :return: Thermal decorrelation
        """

        thermal_corr = (1 / (1 + 1 / snr))
        return thermal_corr

    @staticmethod
    def volume_correlation(eps, penetration_depth, baseline, slant_range, wavelength, incidence_angle):
        """
        Computes the volume decorrelation based on the penetration depth
        :param eps: permittivity
        :param penetration_depth: penetration depth [m]
        :param baseline: baseline [m]
        :param slant_range: slant range [m]
        :param wavelength: wavelength [m]
        :param incidence_angle: incidence angle [rad]
        :return: Volume decorrelation
        """

        corr = 1 / np.sqrt(1 + (2 * np.pi * np.sqrt(eps) * penetration_depth * baseline /
                                (slant_range * wavelength * np.tan(incidence_angle))) ** 2)
        return corr

    @staticmethod
    def critical_baseline(wavelength, slant_range, incidence_angle, bandwidth):
        """
        Calculate the critical baseline
        :param wavelength: wavelength [m]
        :param slant_range: slant range [m]
        :param incidence_angle: incidence angle [rad]
        :param bandwidth: range bandwidth [Hz]
        :return: critical baseline [m]
        """

        critical_baseline = bandwidth * wavelength * slant_range * np.tan(incidence_angle) / C
        return critical_baseline

    @staticmethod
    def nesn2snr(nesn, sigma0):
        nesn2snr = sigma0 / nesn
        return nesn2snr

    @staticmethod
    def coherence2phasenoise(coherence, effective_looks):
        """
        Estimates the Standard Deviation of the Interferometric Phase Error from coherence, based on Cramer-Rao
        :param coherence: coherence of the interferogram
        :param effective_looks: number of effective looks
        :return: standard deviation of the phase of the interferogram [rad]
        """
        phase_noise = np.sqrt((1. - coherence ** 2.) / (2. * effective_looks * coherence ** 2.))
        return phase_noise

    @staticmethod
    def get_output_doppler(prf, na):
        """
        Get the output Doppler frequency
        :param prf: pulse repetition frequency [Hz]
        :param na: number of azimuth samples
        :return: output Doppler frequency [Hz]
        """

        doppler = -prf / 2 + np.arange(na) / na * prf
        return doppler

    @staticmethod
    def get_power_in_band(pat, na_bw, idx_az):
        """
        Get the power in the band
        :param pat: 2D antenna pattern
        :param na_bw: number of azimuth bandwidths
        :param idx_az: azimuth index
        :return: Power in the band
        """

        pow_abs = abs(pat)
        pow_in_band = np.roll((np.cumsum(pow_abs) - np.roll(np.cumsum(pow_abs), na_bw)),
                            -int(0.5 * na_bw))
        return pow_in_band[idx_az]

    @staticmethod
    def regular_2d_interp(data, na_int, nr_int):
        """
        Regular 2D interpolation
        :param data: 2D data array
        :param na_int: azimuth index
        :param nr_int: range index
        :return: Interpolated data
        """

        n_azimuth = len(data[:, 0])
        n_range = len(data[0, :])
        fint = RegularGridInterpolator((np.arange(n_azimuth), np.arange(n_range)), data,
                                       bounds_error=False, fill_value=0)
        return fint((na_int, nr_int))

    @staticmethod
    def compute_antenna_gain(length, height, wavelength, efficiency):
        """
        Compute the antenna gain
        :param length: length of the antenna [m]
        :param height: height of the antenna [m]
        :param wavelength: wavelength [m]
        :param efficiency: antenna efficiency [dB]
        :return: antenna gain [dB]
        """

        factor = 4 * np.pi / wavelength ** 2 * efficiency
        gain = factor * length * height
        return gain

    @staticmethod
    def get_element_pattern(th_out, wavelength, length):
        """
        Get the element pattern of the antenna
        :param th_out: output angle [rad]
        :param wavelength: wavelength [m]
        :param length: length of the antenna [m]
        :return: element pattern
        """

        th3db = wavelength / length
        return np.sinc(np.sin(th_out) / th3db) + 0j

# end of file
