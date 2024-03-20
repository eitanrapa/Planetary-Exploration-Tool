// -*- c++ -*-
//
// the pet development team
// (c) 2023-2024 all rights reserved

// code guard
#if !defined(pet_py_external_h)
#define pet_py_external_h


// STL
#include <string>

// journal
#include <pyre/journal.h>


// pybind support
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


// type aliases
namespace pet::py {
    // import {pybind11}
    namespace py = pybind11;
    // get the special {pybind11} literals
    using namespace py::literals;

    // sizes of things
    using size_t = std::size_t;
    // strings
    using string_t = std::string;
}


#endif

// end of file
