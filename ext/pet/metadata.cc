// -*- C++ -*-
//
// the pet development team
// (c) 2023-2024 all rights reserved
//


// externals
#include <Python.h>
// my declarations
#include "metadata.h"


// copyright
const char * const
pet::extension::
copyright__name__ = "copyright";

const char * const
pet::extension::
copyright__doc__ = "the project copyright string";

PyObject *
pet::extension::
copyright(PyObject *, PyObject *)
{
    // the value
    const char * const copyright_note =
        "pet: (c) 2023-2024 the pet development team";

    // build a python string and return
    return Py_BuildValue("s", copyright_note);
}


// version
const char * const
pet::extension::
version__name__ = "version";

const char * const
pet::extension::
version__doc__ = "the project version string";

PyObject *
pet::extension::
version(PyObject *, PyObject *)
{
        return Py_BuildValue("s", "0.0.1");
}


// end of file