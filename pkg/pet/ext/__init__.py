# -*- coding: utf-8 -*-
#
# the pet development team
# (c) 2023-2025 all rights reserved

# attempt to
try:
    # pull the extension module
    from . import pet as libpet
# if this fails
except ImportError:
    # indicate the bindings are not accessible
    libpet = None


# end of file 
