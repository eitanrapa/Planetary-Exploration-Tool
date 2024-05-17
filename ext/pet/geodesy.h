// -*- C++ -*-
//
// the pet development team
// (c) 2023-2024 all rights reserved
//

#if !defined(pet_extension_geodesy_h)
#define pet_extension_geodesy_h

// place everything in my private namespace
namespace pet {
    namespace extension {

        // cartesian: get cartesian coordinates
        extern const char * const cartesian__name__;
        extern const char * const cartesian__doc__;
        PyObject * cartesian(PyObject *, PyObject *);

        // geodetic: get geodetic coordinates
        extern const char * const geodetic__name__;
        extern const char * const geodetic__doc__;
        PyObject * geodetic(PyObject *, PyObject *);

    } // of namespace extension`
} // of namespace pet

#endif

// end of file