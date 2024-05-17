#include <iostream>
#include <cmath>
#include <iomanip>

#define zero 0.0
#define one 1.0

#ifndef CARTESIAN_H
#define CARTESIAN_H

#ifndef MYSTRUCT_DEFINED
#define MYSTRUCT_DEFINED

struct CartesianPoint{
double x;
double y;
double z;
};

struct GeodeticPoint{
double latitude;
double longitude;
double height;
};

struct TriaxialEllipsoid{
double a;
double b;
double c;
};

#endif // MYSTRUCT_DEFINED

// Function declarations
namespace pet {
    // cartesian
    void cartesian(const TriaxialEllipsoid& te, const GeodeticPoint& gp, CartesianPoint& cp);
}

#endif