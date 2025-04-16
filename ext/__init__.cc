// -*- C++ -*-
//
// the pet development team
// (c) 2023-2025 all rights reserved


// external
#include "external.h"
// namespace setup
#include "forward.h"


// the module entry point
PYBIND11_MODULE(pet, m)
{

    // the doc string
    m.doc() = "the libpet bindings";

    // Geodetic to Cartesian function
    m.def("cartesian", &pet::cartesian, "Compute Cartesian coordinates",
          py::arg("te"), py::arg("gp"), py::arg("cp"));

    // Cartesian to Geodetic function
    m.def("geodetic", &pet::geodetic, "Compute Geodetic coordinate",
        py::arg("te"), py::arg("cp"), py::arg("gp"), py::arg("tol"));

    // bind the opaque types
    pet::py::opaque(m);
    // register the exception types
    pet::py::exceptions(m);
    // version info
    pet::py::version(m);
    //bind CartesianPoint
    pet::py::cartesianPoint(m);
    //bind GeodeticPoint
    pet::py::geodeticPoint(m);
    //bind TriaxialEllipsoid
    pet::py::triaxialEllipsoid(m);

    return;
}


// end of file
