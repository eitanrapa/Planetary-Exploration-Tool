# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

from .Conversions import Conversions as conversions

# attempt to
try:
    # pull the extension module
    from . import pet as libpet
# if this fails
except ImportError:
    # indicate the bindings are not accessible
    libpet = None


# end of file 
