#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import numpy as np
from itertools import combinations


class OnboardProcessing(pet.component):
    """
    Class that encapsulates a standalone simulation to
    """

    def __init__(self, planet, conops, acquisition_cadence, acquisition_radius, **kwargs):
        super().__init__(**kwargs)
        self.planet = planet
        self.conops = conops
        self.acquisition_cadence = acquisition_cadence
        self.acquisition_radius = acquisition_radius

    def extremal_distribution(self, n):
        """

        """

        r = self.acquisition_radius
        theta = np.random.uniform(low=0, high=2 * np.pi, size=n)

        px = r * np.cos(theta)
        py = r * np.sin(theta)
        p = np.asarray([px, py]).T

        return p

    def gaussian_distribution(self, n):
        """

        """

        r = self.acquisition_radius
        px = np.random.normal(0, r, n)
        py = np.random.normal(0, r, n)
        p = np.asarray([px, py]).T
        return p

    def uniform_distribution(self, n):
        """

        """

        r = self.acquisition_radius
        r = r ** 2
        length = np.sqrt(np.random.uniform(low=0, high=r, size=n))
        angle = np.pi * np.random.uniform(low=0, high=2, size=n)

        px = length * np.cos(angle)
        py = length * np.sin(angle)
        p = np.asarray([px, py]).T

        return p

    def get_baseline(self, position1, position2):
        """

        """

        return np.linalg.norm([position1[0] - position2[0], position1[1] - position2[1]])

    def find_intersection(self, p1, p2, p3, slope_perpendicular):
        """

        """

        # Line between p1 and p2
        slope_line = (p2[1] - p1[1]) / (p2[0] - p1[0])
        intercept_line = p1[1] - slope_line * p1[0]

        # Line from p3 with slope perpendicular
        intercept_perpendicular = p3[1] - slope_perpendicular * p3[0]

        # Solve for x, y intersection
        x_intersect = (intercept_perpendicular - intercept_line) / (slope_line - slope_perpendicular)
        y_intersect = slope_line * x_intersect + intercept_line

        return (x_intersect, y_intersect)

    def get_baseline_perp(self, position1, position2, surface_point):
        """

        """

        p1 = [position1, position2][
            np.argmax([np.linalg.norm(point - surface_point) for point in [position1, position2]])]
        p2 = surface_point
        p3 = [position1, position2][
            np.argmin([np.linalg.norm(point - surface_point) for point in [position1, position2]])]

        # Calculate the slope of the original line
        slope_line = (p2[1] - p1[1]) / (p2[0] - p1[0])

        # Perpendicular slope is the negative reciprocal of the original slope
        slope_perpendicular = -1 / slope_line

        # Find the intersection point
        intersection = self.find_intersection(p1, p2, p3, slope_perpendicular)

        return np.linalg.norm(np.asarray(intersection) - np.asarray(p3))

    def classic_stack(self, current_acquisitions, current_times, new_point, new_time, n_acquisitions_stored):
        """
        Stack new acquisitions with timestamps, maintaining only the last `n_acquisitions_stored` acquisitions.
        """

        # Remove the oldest acquisition and timestamp to make space for the new one
        current_acquisitions = np.delete(current_acquisitions, 0, axis=0)
        current_acquisitions = np.append(current_acquisitions, [new_point], axis=0)
        current_times = np.delete(current_times, 0, axis=0)
        current_times = np.append(current_times, [new_time], axis=0)

        return current_acquisitions, current_times

    def lowest_baseline_pairings(self, current_acquisitions, current_times, downlink_pair_capacity, surface_point,
                                 pairs_list):
        """
        Select pairs with the lowest baseline values, associating them with acquisition times.
        """

        possible_pairs = list(combinations(current_acquisitions, 2))
        current_times_pairs = [(current_times[i], current_times[j]) for i, j in
                               combinations(range(len(current_times)), 2)]

        # Compute baseline values for each pair
        pair_baselines = [(pair, self.get_baseline_perp(pair[0], pair[1], surface_point)) for pair in possible_pairs]

        # Sort pairs by their baseline values
        sorted_pairs = sorted(pair_baselines, key=lambda x: x[1])
        sorted_times = [x for _, x in sorted(zip(pair_baselines, current_times_pairs), key=lambda pair: pair[1])]

        sorted_current_times_pairs_isolated = [element[0] for element in sorted_times]
        sorted_pairs_isolated = [element[0] for element in sorted_pairs]

        sorted_pairs_return = []
        sorted_pairs_times_return = []

        i = 0

        while len(sorted_pairs_return) < downlink_pair_capacity:
            if len(pairs_list) == 0:
                sorted_pairs_return.append(sorted_pairs_isolated[i])
                sorted_pairs_times_return.append(sorted_current_times_pairs_isolated[i])
            elif not (np.any(np.all(pairs_list == np.array(sorted_pairs_isolated[i]), axis=(1, 2)))):
                sorted_pairs_return.append(sorted_pairs_isolated[i])
                sorted_pairs_times_return.append(sorted_current_times_pairs_isolated[i])
            i += 1
        return sorted_pairs_return, sorted_pairs_times_return

    def simulate_acquisitions(self, track_number, n_acquisitions_stored, downlink_pair_capacity,
                              total_time, distribution):
        """
        Simulate acquisitions based on a distribution, tracking acquisition and pairing times.
        """

        # Select distribution
        if distribution == "extremal":
            p = self.extremal_distribution(n=int(total_time / self.conops.orbit_cycle))
        elif distribution == "uniform":
            p = self.uniform_distribution(n=int(total_time / self.conops.orbit_cycle))
        elif distribution == "gaussian":
            p = self.gaussian_distribution(n=int(total_time / self.conops.orbit_cycle))
        else:
            raise Exception("input correct distribution")

        # Get track times and initializations
        times = self.conops.get_five_tracks()
        pairings = []
        pairings_times = []
        baselines_perp = []
        current_acquisitions = p[:n_acquisitions_stored]
        current_times = [times[track_number] * self.conops.orbit_cycle * i for i in range(n_acquisitions_stored)]
        surface_point = [10000, -10000]

        for i in range(len(p[n_acquisitions_stored:])):
            # Calculate acquisition time
            new_time = times[track_number] * self.conops.orbit_cycle * (i + n_acquisitions_stored)

            # Stack acquisition and update times
            current_acquisitions, current_times = self.classic_stack(
                current_acquisitions, current_times, p[n_acquisitions_stored + i], new_time, n_acquisitions_stored
            )

            # Get new pairings with baseline and associated times
            current_pairings, current_pairing_times = self.lowest_baseline_pairings(
                current_acquisitions=current_acquisitions,
                current_times=current_times,
                downlink_pair_capacity=downlink_pair_capacity,
                surface_point=surface_point,
                pairs_list=pairings
            )

            # Store pairings and times
            pairings.extend(current_pairings)
            pairings_times.extend(current_pairing_times)
            baselines_perp.extend(
                [self.get_baseline_perp(pairing[0], pairing[1], surface_point) for pairing in current_pairings]
            )

        return pairings_times, baselines_perp

# end of file
