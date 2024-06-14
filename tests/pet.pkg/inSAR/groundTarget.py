#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet

# Make a ground target
ground = pet.insar.groundTarget(x=10, y=10, z=10)

# Get the positions
print(ground.get_position())

# end of file
