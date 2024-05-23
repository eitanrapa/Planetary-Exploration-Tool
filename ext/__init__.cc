// -*- C++ -*-
//
// the pet development team
// (c) 2023-2024 all rights reserved


// external
#include "external.h"
// namespace setup
#include "forward.h"
// conversion code
#include <pet/pet.h>

namespace py = pybind11;

// the module entry point
PYBIND11_MODULE(pet, m)
{
    // the doc string
    m.doc() = "the libpet bindings";

    // CartesianPoint struct
    py::class_<CartesianPoint>(m, "CartesianPoint")
        .def(py::init<>())
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"))
        .def_readwrite("x", &CartesianPoint::x)
        .def_readwrite("y", &CartesianPoint::y)
        .def_readwrite("z", &CartesianPoint::z);

    // GeodeticPoint struct
    py::class_<GeodeticPoint>(m, "GeodeticPoint")
        .def(py::init<>())
        .def(py::init<double, double, double>(), py::arg("latitude"), py::arg("longitude"), py::arg("height"))
        .def_readwrite("latitude", &GeodeticPoint::latitude)
        .def_readwrite("longitude", &GeodeticPoint::longitude)
        .def_readwrite("height", &GeodeticPoint::height);

    // Triaxial Ellipsoid struct
    py::class_<TriaxialEllipsoid>(m, "TriaxialEllipsoid")
        .def(py::init<>())
        .def(py::init<double, double, double>(), py::arg("a"), py::arg("b"), py::arg("c"))
        .def_readwrite("a", &TriaxialEllipsoid::a)
        .def_readwrite("b", &TriaxialEllipsoid::b)
        .def_readwrite("c", &TriaxialEllipsoid::c);

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
}


// end of file
