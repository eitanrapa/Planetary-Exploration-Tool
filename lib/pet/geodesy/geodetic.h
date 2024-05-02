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

int cart2geod(const TriaxialEllipsoid& te, const CartesianPoint& cp, GeodeticPoint& gp, double tol);
int bisect2(double x0,double y0, double cz, double tol,double &m,double &Gm);
int bisect3(double x0,double y0,double z0, double cx, double cy, double tol,double &m,double &Hm);
void xyz2fl(double x,double y,double z, double mee, double mex, double &pr,double &lr);
void fl_octal(double x,double y,double z,double &pr,double &lr);
int sign(double t);