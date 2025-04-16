// -*- C++ -*-
//
// the pet development team
// (c) 2023-2025 all rights reserved


// external
#include "external.h"
// namespace setup
#include "forward.h"


void pet::py::geodeticPoint(py::module &m) {

    // type aliases
    using geodeticPoint_t = pet::GeodeticPoint;

    // GeodeticPoint struct
    py::class_<geodeticPoint_t>(m, "GeodeticPoint")
        .def(py::init<double, double, double>(), py::arg("latitude"), py::arg("longitude"), py::arg("height"))
        .def_readwrite("latitude", &geodeticPoint_t::latitude)
        .def_readwrite("longitude", &geodeticPoint_t::longitude)
        .def_readwrite("height", &geodeticPoint_t::height);

return;
}


// end of file
