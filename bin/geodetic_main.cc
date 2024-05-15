#include "geodetic.h"
#include <iostream>
#include <iomanip>

int main()
{
 	int iter = 1;		//  iterations counter
    CartesianPoint first_point = {1860981.1161016132, 1961008.4466714826, 3935221.2408977931};
    GeodeticPoint geo_point = {0.0, 0.0, 0.0};
    TriaxialEllipsoid earth = {6378172.0000000000, 6378103.0000000000, 6356753.0000000000};
    double error = 1.0e-7; 	//  convergence tolerance

//  === Call the principal function of the method ===

    iter=cart2geod(earth, first_point, geo_point, error);  // returns number of bisection iterations
    
//  =================================================

    std::cout << geo_point.latitude << ", " << geo_point.longitude << ", " << std::setprecision (15) << geo_point.height << std::endl;
    std::cout << iter << std::endl;
}