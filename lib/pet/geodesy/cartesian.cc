#include "cartesian.h"

static auto pih = M_PI_2;
static auto pi = M_PI;

void pet::cartesian(const TriaxialEllipsoid& te, const GeodeticPoint& gp, CartesianPoint& cp)
//  Computes the direct transformation of geodetic to Cartesian coordinates
//  Angular coordinates in radians
{
    double Ex2,Ee2,lex2,lee2,ax2,ay2,b2,kx,ky,kx2,ky2,cx,cy,cz,mex,mee;
    double ps,pc,ls,lc, D, N;
    double longitude_rad, latitude_rad;

	ax2=te.a*te.a;
    ay2=te.b*te.b;
    b2=te.c*te.c;
    Ex2=ax2-b2;
    Ee2=ax2-ay2;
    kx=Ex2/te.a;
    ky=(ay2-b2)/te.b;
    kx2=kx*kx;
    ky2=ky*ky;
    cx=ax2/b2;
    cy=ay2/b2;
    cz=ax2/ay2;
    lex2=Ex2/ax2;
    lee2=Ee2/ax2;
    mex=one-lex2;
    mee=one-lee2;

    latitude_rad = gp.latitude*(pi/180.0);
	longitude_rad = gp.longitude*(pi/180.0);

	ps=std::sin(latitude_rad);
	pc=std::cos(latitude_rad);
	ls=std::sin(longitude_rad);
	lc=std::cos(longitude_rad);

	D=one-lex2*ps*ps-lee2*pc*pc*ls*ls;
	N=te.a/std::sqrt(D);

    cp.x=(N+gp.height)*pc*lc;
    cp.y=(N*mee+gp.height)*pc*ls;
    cp.z=(N*mex+gp.height)*ps;
}
