// -*- C++ -*-
//
// the pet development team
// (c) 2023-2024 all rights reserved
//

// externals
#include <Python.h>
#include <pyre/journal.h>
// my api
#include <pet/pet.h>
// my declarations
#include "geodesy.h"


// cartesian
const char * const
pet::extension::
cartesian__name__ = "cartesian";

const char * const
pet::extension::
cartesian__doc__ = "cartesian coordinates";

PyObject *
pet::extension::
cartesian(PyObject *, PyObject * args)
{
    CartesianPoint first_point = {0.0, 0.0, 0.0};
    GeodeticPoint geo_point;
    TriaxialEllipsoid earth;

    // Parse the arguments from Python
    if (!PyArg_ParseTuple(args, "dddddd",
                          &earth.a, &earth.b, &earth.c,
                          &geo_point.latitude, &geo_point.longitude, &geo_point.height)) {
        return nullptr;  // Return null if parsing failed
    }

    // all done
    pet::cartesian(earth, geo_point, first_point);
    return Py_BuildValue("items", first_point.x, first_point.y, first_point.z);
}

// geodetic
const char * const
pet::extension::
geodetic__name__ = "geodetic";

const char * const
pet::extension::
geodetic__doc__ = "geodetic coordinates";

PyObject *
pet::extension::
geodetic(PyObject *, PyObject * args)
{
 	int iter = 1;		//  iterations counter
    CartesianPoint first_point;
    GeodeticPoint geo_point = {0.0, 0.0, 0.0};
    TriaxialEllipsoid earth;
    double error; 	//  convergence tolerance

    // Parse the arguments from Python
    if (!PyArg_ParseTuple(args, "ddddddd",
                          &earth.a, &earth.b, &earth.c,
                          &first_point.x, &first_point.y, &first_point.z,
                          &error)) {
        return nullptr;  // Return null if parsing failed
    }

    // all done
    iter =  pet::geodetic(earth, first_point, geo_point, error);
    return Py_BuildValue("items", iter, geo_point.latitude, geo_point.longitude, geo_point.height);
}

// end of file