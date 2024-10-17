#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
import numpy as np
from itertools import combinations


class AcquisitionCadence(pet.component):
    """

    """

    def __init__(self, planet, instrument, deformation_map, stored_acquisition_number, acquisition_radius,
                 acquisitions, **kwargs):
        super().__init__(**kwargs)
        self.planet = planet
        self.instrument = instrument
        self.deformation_map = deformation_map
        self.stored_acquisition_number = stored_acquisition_number
        self.acquisition_radius = acquisition_radius
        self.acquisitions = acquisitions

    def extremal_distribution(self, n):
        """

        """

        r = self.acquisition_radius
        theta = np.random.uniform(low=0, high=2 * np.pi, size=n)

        px = r * np.cos(theta)
        py = r * np.sin(theta)
        p = np.asanyarray([px, py]).T

        return p

    def gaussian_distribution(self, n):
        """

        """

        r = self.acquisition_radius
        px = np.random.normal(0, r, n)
        py = np.random.normal(0, r, n)
        p = np.asanyarray([px, py]).T
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
        p = np.asanyarray([px, py]).T

        return p

    def classic_stack(self, current_acquisitions, new_point):
        """

        """

        temp_set = current_acquisitions
        temp_set = np.delete(temp_set, 0, axis=0)
        temp_set = np.append(temp_set, [new_point], axis=0)
        current_acquisitions = temp_set

        return current_acquisitions

    def lowest_baseline_pairings(self, current_acquisitions, n_pairs_per_acquisition, surface_point, pairs_list):
        """

        """

        possible_pairs = list(combinations(current_acquisitions, 2))

        # Compute baseline values for each pair
        pair_baselines = [(pair, self.get_baseline_perp(pair[0], pair[1], surface_point)) for pair in possible_pairs]

        # Sort pairs by their baseline values
        sorted_pairs = sorted(pair_baselines, key=lambda x: x[1])
        sorted_pairs_isolated = [element[0] for element in sorted_pairs]

        sorted_pairs_return = []
        i = 0

        while len(sorted_pairs_return) < n_pairs_per_acquisition:
            if len(pairs_list) == 0:
                sorted_pairs_return.append(sorted_pairs_isolated[i])
            elif (np.any(np.all(pairs_list == np.array(sorted_pairs_isolated[i]), axis=(1, 2)))) == False:
                sorted_pairs_return.append(sorted_pairs_isolated[i])
            i += 1
        return sorted_pairs_return

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

        return np.linalg.norm(np.asanyarray(intersection) - np.asanyarray(p3))

    def get_baselines(self, track_number, n_passes, distribution):
        """

        """

        if distribution == "extremal":
            p = self.extremal_distribution(n=n_passes)
        elif distribution == "uniform":
            p = self.extremal_distribution(n=n_passes)
        elif distribution == "gaussian":
            p = self.extremal_distribution(n=n_passes)
        else:
            raise Exception("input correct distribution")

        # Get the times defining the first five tracks
        times = self.instrument.get_five_tracks()

        pairings = []
        baselines = []
        baselines_perp = []
        current_acquisitions = p[:n_acquisitions]
        for i in range(len(p[n_acquisitions:])):
            current_acquisitions = classic_stack(current_acquisitions, p[n_acquisitions + i])
            current_pairings = lowest_baseline_pairings(current_acquisitions, n_pairs_per_acquisition,
                                                        surface_point, pairings)
            pairings.extend(current_pairings)
            baselines.extend([get_baseline(pairing[0], pairing[1]) for pairing in current_pairings])
            baselines_perp.extend(
                [get_baseline_perp(pairing[0], pairing[1], surface_point) for pairing in current_pairings])
        return current_acquisitions, pairings[burn_in:], baselines[burn_in:], baselines_perp[burn_in:]

# end of file
