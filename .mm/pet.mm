# -*- makefile -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved


# pet consists of a python package
pet.packages := pet.pkg
# libraries
pet.libraries := pet.lib
# python extensions
pet.extensions := pet.ext
# a ux bundle
pet.webpack :=
# and some tests
pet.tests := pet.lib.tests pet.ext.tests pet.pkg.tests


# load the packages
include $(pet.packages)
# the libraries
include $(pet.libraries)
# the extensions
include $(pet.extensions)
# the ux
include $(pet.webpack)
# and the test suites
include $(pet.tests)


# end of file
