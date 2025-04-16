#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

import pet
import numpy as np
from itertools import combinations


class InterferogramPairOptimization(pet.component, family="pet.dataAcquisition.interferogramPairOptimization",
                                    implements=pet.protocols.dataAcquisition):
    """
    Class that encapsulates the onboard interferogram selecting algorithm
    """

    planet = pet.protocols.planet()
    planet.doc = "target planet"

    campaign = pet.protocols.campaigns.orbiter()
    campaign.doc = "campaign (must be an oribiter)"

    acquisition_cadence = pet.properties.path()
    acquisition_cadence.doc = "acquisition cadence path"

    acquisition_radius = pet.properties.float()
    acquisition_radius.doc = "radius of acquisition tube"

    def get_baseline(self, position1, position2):
        """
        Get the baseline (distance between two points)
        :param position1: Position of first point in x, y, z coordinates
        :param position2: Positions of second point in x, y, z coordinates
        :return: The baseline
        """

        return np.linalg.norm([position1[0] - position2[0], position1[1] - position2[1]])

    def find_intersection(self, p1, p2, p3):
        """
        Find the intersection between a line originating from p2 towards the line between p1 and p3
        :param p1: Position of first point in x, y, z coordinates
        :param p2: Position of second point in x, y, z coordinates
        :param p3: Position of third point in x, y, z coordinates
        :return: x and y of intersect
        """

        # Line between p1 and p2
        slope_line = (p2[1] - p1[1]) / (p2[0] - p1[0])
        intercept_line = p1[1] - slope_line * p1[0]

        # Perpendicular slope is the negative reciprocal of the original slope
        slope_perpendicular = -1 / slope_line

        # Line from p3 with slope perpendicular
        intercept_perpendicular = p3[1] - slope_perpendicular * p3[0]

        # Solve for x, y intersection
        x_intersect = (intercept_perpendicular - intercept_line) / (slope_line - slope_perpendicular)
        y_intersect = slope_line * x_intersect + intercept_line

        return x_intersect, y_intersect

    def get_baseline_perp(self, position1, position2, surface_point):
        """
        Get the perpendicular component of the baseline between two points and a point on the surface
        :param position1: Position of first point in x, y, z coordinates
        :param position2: Position of second point in x, y, z coordinates
        :param surface_point: Position of surface point in x, y, z coordinates
        :return: Perpendicular baseline component
        """

        # Define the points
        p1 = [position1, position2][
            np.argmax([np.linalg.norm(point - surface_point) for point in [position1, position2]])]
        p2 = surface_point
        p3 = [position1, position2][
            np.argmin([np.linalg.norm(point - surface_point) for point in [position1, position2]])]

        # Find the intersection point
        intersection = self.find_intersection(p1=p1, p2=p2, p3=p3)

        # Return the distance
        return np.linalg.norm(np.asarray(intersection) - np.asarray(p3))

    def classic_stack(self, current_acquisitions, current_times, new_point, new_time):
        """
        Stack new acquisitions with timestamps
        :param current_acquisitions: Current acquisitions
        :param current_times: Current times of acquisitions
        :param new_point: New point to be added
        :param new_time: Time of new point
        :return: Current acquisitions and current times after replacement
        """

        # Remove the oldest acquisition and timestamp to make space for the new one
        current_acquisitions = np.delete(current_acquisitions, 0, axis=0)
        current_acquisitions = np.append(current_acquisitions, [new_point], axis=0)
        current_times = np.delete(current_times, 0, axis=0)
        current_times = np.append(current_times, [new_time], axis=0)

        # Return the current acquisitions and times
        return current_acquisitions, current_times

    def lowest_baseline_pairings(self, current_acquisitions, current_times, downlink_pair_capacity, surface_point,
                                 pairs_list):
        """
        Select pairs with the lowest baseline values, associating them with acquisition times.
        :param current_acquisitions: Current acquisitions
        :param current_times: Current times of acquisitions
        :param downlink_pair_capacity: How many interferograms can be sent at a time
        :param surface_point: Arbitrary point on a surface
        :param pairs_list: Pairs that have already been made
        :return: Pairs to be made, and times at which they are made
        """

        # Make the possible list of pairs
        possible_pairs = list(combinations(current_acquisitions, 2))

        # Get the same pairs but times at which the acquisitions are made
        current_times_pairs = [(current_times[i], current_times[j]) for i, j in
                               combinations(range(len(current_times)), 2)]

        # Compute baseline values for each pair
        pair_baselines = [(pair, self.get_baseline_perp(pair[0], pair[1], surface_point)) for pair in possible_pairs]

        # Sort pairs by their baseline values
        sorted_pairs = sorted(pair_baselines, key=lambda x: x[1])

        # Sort times in the same way
        sorted_times = [x for _, x in sorted(zip(pair_baselines, current_times_pairs), key=lambda pair: pair[1])]

        # Isolate the current times and current pairs
        sorted_current_times_pairs_isolated = [[element[0], element[1]] for element in sorted_times]
        sorted_pairs_isolated = [element[0] for element in sorted_pairs]

        # Make a list of pairs to return
        sorted_pairs_return = []
        sorted_pairs_times_return = []

        i = 0

        # Add pairs while the returned list is less than the downlink capacity
        while len(sorted_pairs_return) < downlink_pair_capacity:
            if len(pairs_list) == 0:
                sorted_pairs_return.append(sorted_pairs_isolated[i])
                sorted_pairs_times_return.append(sorted_current_times_pairs_isolated[i])
            elif not (np.any(np.all(pairs_list == np.array(sorted_pairs_isolated[i]), axis=(1, 2)))):
                sorted_pairs_return.append(sorted_pairs_isolated[i])
                sorted_pairs_times_return.append(sorted_current_times_pairs_isolated[i])
            i += 1

        # Return the pairs and times
        return sorted_pairs_return, sorted_pairs_times_return

    def simulate_acquisitions_per_track(self, track_number, n_acquisitions_stored, downlink_pair_capacity,
                                        total_time, distribution):
        """
        Simulate acquisitions based on a distribution, tracking acquisition and pairing times, for a specific track
        :param track_number: The number (from 1 to 5) of the track to simulate
        :param n_acquisitions_stored: The number of acquisitions that can be stored on-board the spacecraft
        :param downlink_pair_capacity: How many interferograms can be sent at a time
        :param total_time: How much time to run simulation for
        :param distribution: Which distribution to choose from for positions of the satellite
        :return: Pairing times and baselines for the pairs made
        """

        # Select distribution
        if distribution == "extremal":
            p = self.extremal_distribution(n=int(total_time / self.campaign.orbit_cycle))
        elif distribution == "uniform":
            p = self.uniform_distribution(n=int(total_time / self.campaign.orbit_cycle))
        elif distribution == "gaussian":
            p = self.gaussian_distribution(n=int(total_time / self.campaign.orbit_cycle))
        else:
            raise Exception("input correct distribution")

        # Get track times and initializations
        times = self.campaign.get_five_tracks()
        pairings = []
        pairings_times = []
        baselines_perp = []
        current_acquisitions = p[:n_acquisitions_stored]
        current_times = [times[track_number] * self.campaign.orbit_cycle * i for i in range(n_acquisitions_stored)]

        # Define an arbitrary surface point
        surface_point = [10000, -10000]

        # Iterate through the samples from the distribution
        for i in range(len(p[n_acquisitions_stored:])):
            # Calculate acquisition time
            new_time = times[track_number] * self.campaign.orbit_cycle * (i + n_acquisitions_stored)

            # Stack acquisition and update times
            current_acquisitions, current_times = self.classic_stack(
                current_acquisitions=current_acquisitions, current_times=current_times,
                new_point=p[n_acquisitions_stored + i], new_time=new_time)

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

        # Return pairing times and baselines
        return pairings_times, baselines_perp

    def simulate_acquisitions(self, n_acquisitions_stored, downlink_pair_capacity,
                              total_time, distribution):
        """
        Simulate acquisitions based on a distribution, tracking acquisition and pairing times
        :param n_acquisitions_stored: The number of acquisitions that can be stored on-board the spacecraft
        :param downlink_pair_capacity: How many interferograms can be sent at a time
        :param total_time: How much time to run simulation for
        :param distribution: Which distribution to choose from for positions of the satellite
        :return: Pairing times and baselines for the pairs made
        """

        # Make a list of total pairing times and baselines for each track
        pairing_times_total = []
        baselines_perp_total = []

        # Iterate through the acquisition cadence file
        for key, value in self.acquisition_cadence.items():
            if value is not None:  # Check if the value is not null
                pairing_times, baselines_perp = self.simulate_acquisitions_per_track(
                    track_number=int(key), n_acquisitions_stored=n_acquisitions_stored,
                    downlink_pair_capacity=downlink_pair_capacity,
                    total_time=total_time, distribution=distribution)
                pairing_times_total.append(pairing_times)
                baselines_perp_total.append(baselines_perp)

        # Return the pairing times and baselines
        return np.asarray(pairing_times_total), np.asarray(baselines_perp_total)


# end of file
