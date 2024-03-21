#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved


"""
Version check
"""


def test():
    # access the {pet} package
    import pet
    # access the {pet} extension
    from pet import libpet
    # verify that the static and current versions match
    assert libpet.version.static() == pet.version()
    # verify that the dynamic and current versions match
    assert libpet.version.dynamic() == pet.version()
    # all done
    return



# main
if __name__ == "__main__":
    # do...
    test()


# end of file

