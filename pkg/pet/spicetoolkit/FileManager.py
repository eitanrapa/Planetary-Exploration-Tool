#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved
import cspyce as spice
import pathlib
import numpy as np


class FileManager:
    """
    Class that manages SPICE kernel loading
    """

    def __init__(self, folder_path):
        self.folder_path = pathlib.PosixPath(folder_path)

    def _furnsh(self, names_list):
        """
        Activate a set of spice kernel files
        :param names_list: List of paths for kernels
        :return: Nothing returned
        """

        # Loop over paths
        for name in names_list:

            # Use spice to furnish
            spice.furnsh(self.folder_path / name)

    def furnsh(self, names_list):
        """
        Helper function for _furnsh
        """

        # Make sure names_list is a numpy array
        names_list = np.asanyarray(names_list)

        # Make sure dimensions is 1
        if names_list.ndim == 0:
            names_list = np.asanyarray([names_list])

        # Check list is made of strings
        for name in names_list:
            print(name)
            if not isinstance(name, str):
                raise Exception("One or more paths in list is not a string")

        self._furnsh(names_list)

        # Catch file not found error
        try:
            return self._furnsh(names_list)
        except OSError:
            print("One or more paths could not be found")

    def clear(self):
        """
        Clear the SPICE kernels
        :return: Nothing returned
        """

        # Clear the kernels
        spice.kclear()

# end of file
