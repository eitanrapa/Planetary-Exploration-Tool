#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# the pet development team
# california institute of technology
# (c) 2023-2025 all rights reserved
#


"""
Version check
"""


def test():
    # access the {pet} extension
    from pet import libpet
    # verify that the static and dynamic versions match
    assert libpet.version.static() == libpet.version.dynamic()
    # all done
    return


# main
if __name__ == "__main__":
    # do...
    test()


# end of file
