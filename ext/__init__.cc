// -*- C++ -*-
//
// the pet development team
// (c) 2023-2024 all rights reserved


// external
#include "external.h"
// namespace setup
#include "forward.h"

#include "pet/geodesy.h"
#include "pet/metadata.h"

namespace py = pybind11;

// the module entry point
PYBIND11_MODULE(pet, m)
{
    // the doc string
    m.doc() = "the libpet bindings";

    // Bind the geodesy functions
    m.def("cartesian", &pet::extension::cartesian, py::return_value_policy::reference);
    m.def("geodetic", &pet::extension::geodetic, py::return_value_policy::reference);

    // bind the opaque types
    pet::py::opaque(m);
    // register the exception types
    pet::py::exceptions(m);
    // version info
    pet::py::version(m);
}


// end of file
