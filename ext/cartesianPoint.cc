// -*- C++ -*-
//
// the pet development team
// (c) 2023-2025 all rights reserved


// external
#include "external.h"
// namespace setup
#include "forward.h"


void pet::py::cartesianPoint(py::module &m) {

    // type aliases
    using cartesianPoint_t = pet::CartesianPoint;

    // CartesianPoint struct
    py::class_<cartesianPoint_t>(m, "CartesianPoint")
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"))
        .def_readwrite("x", &cartesianPoint_t::x)
        .def_readwrite("y", &cartesianPoint_t::y)
        .def_readwrite("z", &cartesianPoint_t::z);

return;
}


// end of file
