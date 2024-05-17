// -*- c++ -*-
//
// the pet development team
// (c) 2023-2024 all rights reserved

// code guard
#if !defined(pet_py_forward_h)
#define pet_py_forward_h


// the {project.name} namespace
namespace pet::py {

    // bindings of opaque types
    void opaque(py::module &);
    // exceptions
    void exceptions(py::module &);
    // version info
    void version(py::module &);
}


#endif

// end of file
