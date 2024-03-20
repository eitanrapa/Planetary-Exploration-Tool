// -*- C++ -*-
//
// the pet development team
// (c) 2023-2024 all rights reserved


// external
#include "external.h"
// namespace setup
#include "forward.h"


// the module entry point
PYBIND11_MODULE(pet, m)
{
    // the doc string
    m.doc() = "the libpet bindings";

    // bind the opaque types
    pet::py::opaque(m);
    // register the exception types
    pet::py::exceptions(m);
    // version info
    pet::py::version(m);
}


// end of file
