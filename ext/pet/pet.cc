// -*- C++ -*-
//
// the pet development team
// (c) 2023-2024 all rights reserved
//

// #include <portinfo>
#include <Python.h>

// the module method declarations
#include "geodesy.h"
#include "metadata.h"

// put everything in my private namespace
namespace pet {
    namespace extension {

        // the module method table
        PyMethodDef module_methods[] = {
            // admin
            { copyright__name__, copyright, METH_VARARGS, copyright__doc__ },
            { version__name__, version, METH_VARARGS, version__doc__ },

            // geodesy
            { cartesian__name__, cartesian, METH_VARARGS, cartesian__doc__ },
            { geodetic__name__, geodetic, METH_VARARGS, geodetic__doc__ },

            // sentinel
            { 0, 0, 0, 0 }
        };

        // the module documentation string
        const char * const __doc__ = "an example of a python extension";

        // the module definition structure
        PyModuleDef module_definition = {
            // header
            PyModuleDef_HEAD_INIT,
            // the name of the module
            "pet",
            // the module documentation string
            __doc__,
            // size of the per-interpreter state of the module; -1 if this state is global
            -1,
            // the methods defined in this module
            module_methods
        };

    } // of namespace extension
} // of namespace pet


// initialization function for the module
// *must* be called PyInit_hello
PyMODINIT_FUNC
PyInit_pet()
{
    // create the module
    PyObject * module = PyModule_Create(&pet::extension::module_definition);
    // and return it
    return module;
}

// end of file