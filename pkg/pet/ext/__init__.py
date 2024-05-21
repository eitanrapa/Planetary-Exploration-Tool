# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2024 all rights reserved

# attempt to
try:
    # pull the extension module
    from . import pet as libpet
    from .Conversions import Conversions as conversions
# if this fails
except ImportError:
    # indicate the bindings are not accessible
    libpet = None


# end of file 
