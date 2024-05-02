#include "cartesian.h"
#include <iostream>

int main()
{
    CartesianPoint first_point = {0.0, 0.0, 0.0};
    GeodeticPoint geo_point = {55.75, 46.5, -1589188.43460915};
    TriaxialEllipsoid earth = {6378172.0000000000, 6378103.0000000000, 6356753.0000000000};

//  === Call the principal function of the method ===

    direct(earth, geo_point, first_point);

//  =================================================

    std::cout << first_point.x << ", " << first_point.y << ", " << first_point.z << std::endl;
}