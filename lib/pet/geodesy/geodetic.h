// -*- C++ -*-
//
// the pet development team
// (c) 2023-2025 all rights reserved
//


// code guard
#pragma once

#include "CartesianPoint.h"
#include "GeodeticPoint.h"
#include "TriaxialEllipsoid.h"

// Function declarations
namespace pet {
    // geodetic
    int geodetic(const TriaxialEllipsoid& te, const CartesianPoint& cp, GeodeticPoint& gp, double tol);
    int bisect2(double x0,double y0, double cz, double tol,double &m,double &Gm);
    int bisect3(double x0,double y0,double z0, double cx, double cy, double tol,double &m,double &Hm);
    void xyz2fl(double x,double y,double z, double mee, double mex, double &pr,double &lr);
    void fl_octal(double x,double y,double z,double &pr,double &lr);
    int sign(double t);
}


// end of file
