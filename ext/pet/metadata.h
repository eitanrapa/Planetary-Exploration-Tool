// -*- C++ -*-
//
// the pet development team
// (c) 2023-2024 all rights reserved
//


#if !defined(pet_extension_metadata_h)
#define pet_extension_metadata_h


// place everything in my private namespace
namespace pet {
    namespace extension {

        // copyright note
        extern const char * const copyright__name__;
        extern const char * const copyright__doc__;
        PyObject * copyright(PyObject *, PyObject *);

        // version
        extern const char * const version__name__;
        extern const char * const version__doc__;
        PyObject * version(PyObject *, PyObject *);

    } // of namespace extension`
} // of namespace pet

#endif

// end of file