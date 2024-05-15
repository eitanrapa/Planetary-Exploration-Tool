#include "geodetic.h"

static auto pih = M_PI_2;
static auto pi = M_PI;

// ---------------------------------------
//  FUNCTIONS
// ---------------------------------------

int cart2geod(const TriaxialEllipsoid& te, const CartesianPoint& cp, GeodeticPoint& gp, double tol)
//  Principal function for the transformation of Cartesian to geodetic coordinates
//  Angular coordinates in radians
{
    int n;
    double Ex2,Ee2,lex2,lee2,ax2,ay2,b2,kx,ky,kx2,ky2,cx,cy,cz,mex,mee;
    double x,y,z,x2,y2,x0,y0,z0,dx,dy,dz;
    double sx,sy,sz,m,Mm;

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

	x=std::abs(cp.x);
	y=std::abs(cp.y);
	z=std::abs(cp.z);
	x2=x*x;
	y2=y*y;
	
//----- Determination of foot point -----	
//       (using bisection method)

	if (z==zero)  {
		if ((x==zero) && (y==zero)) {
			sx=zero;
			sy=zero;
			sz=te.c;
			n=0;
			m=zero;
			Mm=zero;  }
		else if ((ky2*x2+kx2*y2)<(kx2*ky2))  {
			sx=te.a*x/kx;
			sy=te.b*y/ky;
			sz=te.c*std::sqrt(one-(x2/kx2)-(y2/ky2));
			n=0;
			m=zero;
			Mm=zero;  }
		else if (x==zero)  {
			sx=zero;
			sy=te.b;
			sz=zero;
			n=0;
			m=zero;
			Mm=zero;  }
		else if (y==zero)  {
			sx=te.a;
			sy=zero;
			sz=zero;
			n=0;
			m=zero;
			Mm=zero;  }
		else  {
			x0=x/te.a;
			y0=y/te.b;
			n=bisect2(x0,y0,cz,tol,m,Mm);
			sx=cz*x/(cz+m);
			sz=zero;
			sy=y/(one+m);
		}
	}
	else  {
		if ((x==zero) && (y==zero)) {
			sx=zero;
			sy=zero;
			sz=te.c;
			n=0;
			m=zero;
			Mm=zero;  }
		else  {
			x0=x/te.a;
			y0=y/te.b;
			z0=z/te.c;
			n=bisect3(x0,y0,z0,cx,cy,tol,m,Mm);
			sx=cx*x/(cx+m);
			sy=cy*y/(cy+m);
			if ((m<zero) && ((ky2*x2 + kx2*y2) < (kx2*ky2)))  
				sz=te.c*std::sqrt(one-((sx*sx)/ax2)-((sy*sy)/ay2));
			else
				sz=z/(one+m);
		}
	}
	
//- Determination of latitude & longitude of foot point -----	

	xyz2fl(sx,sy,sz,mee,mex,gp.latitude,gp.longitude);

//- Determination of geodetic height -----	

	dx=x-sx;
	dy=y-sy;
	dz=z-sz;
	gp.height=std::sqrt(dx*dx+dy*dy+dz*dz);

	if ((x+y+z)<(sx+sy+sz))   { gp.height=-gp.height; }

//- Adjust latitude & longitude according to signs of (cp.x,cp.y,cp.z)

	fl_octal(cp.x,cp.y,cp.z,gp.latitude,gp.longitude);

//- Convert to degrees

    gp.latitude = gp.latitude*(180.0/pi);
    gp.longitude = gp.longitude*(180.0/pi);

	return n;
}

void xyz2fl(double x,double y,double z, double mee, double mex, double &latitude,double &longitude)
//  Computes the transformation of Cartesian to geodetic coordinates on the surface of the ellipsoid
//  assuming x,y,z are all non-negative
//  Angular coordinates in radians
{
	double nom,den,dex,xme,rot;
	
	nom=mee*z;
	xme=mee*x;
	dex=xme*xme+y*y;
	den=mex*std::sqrt(dex);
	rot=std::sqrt(dex);
	
	if (den==zero)  {
		latitude=pih;
		longitude=zero;
	}
	else  {
		if (nom<=den)  
			{ latitude=std::atan(nom/den); }
		else
			{ latitude=pih-std::atan(den/nom); }
	
	    if (y<=xme)  {
			den=xme+rot;
			longitude=2.*std::atan(y/den);
		}
		else   {
			den=y+rot;
			longitude=pih-2.0*std::atan(xme/den);
		}
	}
}

void fl_octal(double x,double y,double z,double &latitude,double &longitude)
//  Adjusts latitude & longitude according to signs of x,y,z  (angles in radians)
{
	if (z<zero)   {latitude=-latitude;}
	
	if (x>=zero)  {
		if (y<zero)  {longitude=-longitude;}  }
	else  {
		if (y>=zero)  {longitude=pi-longitude;}
		else  {longitude=longitude-pi;}
	}	
}

int bisect2(double x0,double y0,double cz, double tol,double &m,double &Gm)
//  Implements the bisection method on the X-Y plane
{
	int n;
	double d,d1,d2,Gd,g2,MC;

	n=0;
	m=-2.;
	d1=y0-one;
	g2=cz*cz*x0*x0;
	d2=std::sqrt(g2+y0*y0)-one;
	Gd=g2/((cz+d1)*(cz+d1));
	d=0.50*(d2-d1);

	while (d>tol) {
		n++;
		MC=m;
		m=d1+d;
		Gm=(g2/((cz+m)*(cz+m))) + ((y0*y0)/((one+m)*(one+m))) - one;
		
		if  (MC==(m+tol))   return n;
		
		if  (Gm==zero)  { return n; }
		else  {
			if (sign(Gm)==sign(Gd))  {
				d1=m;
				Gd=Gm;
			}
			else
			    { d2=m; }
		}
		d=0.50*(d2-d1);
	}
	n++;
	m=d1+d;
	Gm=(g2/((cz+m)*(cz+m))) + ((y0*y0)/((one+m)*(one+m))) - one;
	
	return n;	
}

int bisect3(double x0,double y0,double z0, double cx, double cy, double tol,double &m,double &Hm)
//  Implements the bisection method in 3D space
{
	int n;
	double d,d1,d2,Hd,g2,g3,MC;

	n=0;
	m=-2.;
	d1=z0-one;
	g2=cx*cx*x0*x0;
	g3=cy*cy*y0*y0;
	d2=std::sqrt(g2+g3+z0*z0)-one;
	Hd=g2/((cx+d1)*(cx+d1)) + g3/((cy+d1)*(cy+d1));
	d=0.50*(d2-d1);

	while (d>tol) {
		n++;
		MC=m;
		m=d1+d;
		Hm=(g2/((cx+m)*(cx+m))) + (g3/((cy+m)*(cy+m))) + ((z0*z0)/((one+m)*(one+m))) - one;
		
		if  (MC==(m+tol))   return n;
		
		if  (Hm==zero)  { return n; }
		else  {
			if (sign(Hm)==sign(Hd))  {
				d1=m;
				Hd=Hm;
			}
			else
			   { d2=m; }
		}
		d=0.50*(d2-d1);
	}
	n++;
	m=d1+d;
	Hm=(g2/((cx+m)*(cx+m))) + (g3/((cy+m)*(cy+m))) + ((z0*z0)/((one+m)*(one+m))) - one;
	
	return n;	
}

int sign(double t)
//  Returns the sign of its argument
{
	int n;
	if (t>zero) n=1;
	else if (t<zero) n=-1;
	else  n=0;

	return n;
}