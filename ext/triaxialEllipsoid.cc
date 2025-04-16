// -*- C++ -*-
//
// the pet development team
// (c) 2023-2025 all rights reserved


// external
#include "external.h"
// namespace setup
#include "forward.h"


void pet::py::triaxialEllipsoid(py::module &m) {

    // type aliases
    using triaxialEllipsoid_t = pet::TriaxialEllipsoid;

    // Triaxial Ellipsoid struct
    py::class_<triaxialEllipsoid_t>(m, "TriaxialEllipsoid")
        .def(py::init<double, double, double>(), py::arg("a"), py::arg("b"), py::arg("c"))
        .def_readwrite("a", &triaxialEllipsoid_t::a)
        .def_readwrite("b", &triaxialEllipsoid_t::b)
        .def_readwrite("c", &triaxialEllipsoid_t::c);

return;
}


// end of file
