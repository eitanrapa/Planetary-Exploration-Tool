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
    // cartesian
    void cartesian(const TriaxialEllipsoid& te, const GeodeticPoint& gp, CartesianPoint& cp);
}


// end of file
