#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet
from pyre.units.SI import meter, second
from pyre.units.angle import degree
import pandas as pd


class Nightingale(pet.component, family="pet.instruments.nightingale", implements=pet.protocols.instrument):

    # Parameter table
    paramsDF = pd.read_excel("/home/erapapor/kraken-bak/Enceladus-Exploration-Twin-files/" +
                             "parameter_tables/Parameters Sband - 200 km - 0.5mx3.5m hga antenna worst.xlsx",
                             skiprows=1, usecols=[0, 1, 2, 3, 4], index_col=0)

# end of file
