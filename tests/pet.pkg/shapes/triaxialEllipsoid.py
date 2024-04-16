#!/usr/bin/env python3
# -*- Python -*-
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

import pet


def test():
    # access the {pet.shapes} sub-package
    import pet.shapes
    # make a triaxial ellipsoid
    tri = pet.shapes.triaxialEllipsoid(name="sphere", a=10, b=10, c=10)
    # verify
    assert tri.pyre_name == "sphere"
    assert tri.a == 10
    assert tri.b == 10
    assert tri.c == 10
    # make another one
    # from configuration file
    tri = pet.shapes.triaxialEllipsoid(name="tri")
    # verify
    assert tri.pyre_name == "tri"
    assert tri.a == 4
    assert tri.b == 5
    assert tri.c == 6
    # all done
    return


# main
if __name__ == "__main__":
    # do...
    test()

# end of file
