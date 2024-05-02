#include <iostream>
#include <cmath>
#include <iomanip>

#define zero 0.0
#define one 1.0

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

// Function declarations

void direct(const TriaxialEllipsoid& te, const GeodeticPoint& gp, CartesianPoint& cp);
