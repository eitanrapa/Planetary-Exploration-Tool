//
// Code for Cartesian to geodetic coordinate conversion
// Authors: Panou & Korakitis, 2019
// Manuscript: "Cartesian to geodetic coordinates conversion
// on a triaxial ellipsoid using the bisection method"
//

#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

// ---------------------------------------
//  Global variables
// ---------------------------------------
long double const zero = 0.0;
long double const one = 1.0;

long double const pi=4.0*atan(one);
long double const pih=pi/2;
long double const d2r=pi/180.0;				// Conversion from deg to rad

long double Ex2,Ee2,lex2,lee2,ax,ay,b,ax2,ay2,b2;
long double kx,ky,kx2,ky2,cx,cy,cz,mex,mee,D,N;

//-----------------------------------------------------------------
// Function prototypes
// ----------------------------------------------------------------
void ellipse();
void direct(long double &xf,long double &yf,long double &zf,long double pr,long double lr,long double hm);
int cart2geod(long double tol,long double xi,long double yi,long double zi,long double &pr,long double &lr,long double &hm);
int bisect2(long double x0,long double y0,long double tol,long double &m,long double &Gm);
int bisect3(long double x0,long double y0,long double z0,long double tol,long double &m,long double &Hm);
void xyz2fl(long double x,long double y,long double z,long double &pr,long double &lr);
void fl_octal(long double x,long double y,long double z,long double &pr,long double &lr);
int sign(long double t);

//------------------------
// Main program function
// -----------------------
int main()
{
 	int iter = 1;		//  iterations counter
 	long double x0,y0,z0,latitude,longitude,altitude;
 	long double phi,lamda,x,y,z,dx,dy,dz,dS,logds;
    long double error = 1.0e-19; 	//  convergence tolerance

/*  === Use appropriate values for semi-axes of ellipsoid === 

	(Virtual body)
	ax=1000000.0000000000;
	ay= 700000.0000000000;
	 b= 300000.0000000000;

    (Oblate Earth)
	ax=6378137.0000000000;
	ay=6378137.0000000000;
	 b=6356752.3142451795;

    (Triaxial Earth)
	ax=6378172.0000000000;
	ay=6378103.0000000000;
	 b=6356753.0000000000;

    (Moon)
	ax=1735550.0000000000;
	ay=1735324.0000000000;
	 b=1734989.0000000000;

*///========================================================
//  Example
	ax = 6378172.0000000000;
	ay = 6378103.0000000000;
	 b = 6356753.0000000000;

//  === Compute secondary ellipsoid parameters ===

	ellipse();
	
//  ==============================================

//  === Enter values of Cartesian coordinates  x0, y0, z0 ===
//  (from keyboard, data file etc.)
//  Example
	x0 =  1860981.1161016132;                                                                                 
	y0 =  1961008.4466714826; 
	z0 =  3935221.2408977931;
	
//  =========================================================

//  === Call the principal function of the method ===

    iter=cart2geod(error,x0,y0,z0,latitude,longitude,altitude);  // returns number of bisection iterations
    
//  =================================================

//  === Optionally, express latitude and longitude in degrees ===

	  phi = latitude/d2r;
	lamda = longitude/d2r;  

//  ================================================================	

//  === Compute transformation error dS =========

	direct(x,y,z,latitude,longitude,altitude);
	
	dx = x0-x;
	dy = y0-y;
	dz = z0-z;
	
	dS = sqrt(dx*dx+dy*dy+dz*dz);
	logds = log10(dS);  // Optionally, in logarithmic scale
	
//  =============================================

//  === Output results as appropriate ===
//  (to screen, data file etc.)
//  Example
	cout<<fixed<<setprecision(9);
	cout<<"\n  Latitude = "<<phi;
	cout<<"\n Longitude = "<<lamda;
	cout<<"\n  Altitude = "<<altitude;
	cout<<setprecision(2)<<endl;
	cout<<"\n log10(ds) = "<<logds;
	cout<<"\nIterations = "<<iter<<endl;
	
//  =========================================================

    return 0;
}


// ---------------------------------------
//  FUNCTIONS
// ---------------------------------------

void ellipse()
//  Computes several ellipsoid parameters
{
    ax2=ax*ax;
    ay2=ay*ay;
    b2=b*b;
    Ex2=ax2-b2;
    Ee2=ax2-ay2;
    kx=Ex2/ax;
    ky=(ay2-b2)/ay;
    kx2=kx*kx;
    ky2=ky*ky;
    cx=ax2/b2;
    cy=ay2/b2;
    cz=ax2/ay2;
    lex2=Ex2/ax2;
    lee2=Ee2/ax2;
    mex=one-lex2;
    mee=one-lee2;
}

void direct(long double &xf,long double &yf,long double &zf,long double latitude,long double longitude,long double height)
//  Computes the direct transformation of geodetic to Cartesian coordinates
//  Angular coordinates in radians
{
    long double ps,pc,ls,lc;
	
	ps=sin(latitude);
	pc=cos(latitude);
	ls=sin(longitude);
	lc=cos(longitude);
	
	D=one-lex2*ps*ps-lee2*pc*pc*ls*ls;
	N=ax/sqrt(D);

    xf=(N+height)*pc*lc;
    yf=(N*mee+height)*pc*ls;
    zf=(N*mex+height)*ps;
}

int cart2geod(long double tol,long double xi,long double yi,long double zi,long double &latitude,long double &longitude,long double &height)
//  Principal function for the transformation of Cartesian to geodetic coordinates
//  Angular coordinates in radians
{
    int n;
    long double x,y,z,x2,y2,x0,y0,z0,dx,dy,dz; 
    long double sx,sy,sz,m,Mm;
	
	x=abs(xi);
	y=abs(yi);
	z=abs(zi);
	x2=x*x;
	y2=y*y;
	
//----- Determination of foot point -----	
//       (using bisection method)

	if (z==zero)  {
		if ((x==zero) && (y==zero)) {
			sx=zero;
			sy=zero;
			sz=b;
			n=0;
			m=zero;
			Mm=zero;  }
		else if ((ky2*x2+kx2*y2)<(kx2*ky2))  {
			sx=ax*x/kx;
			sy=ay*y/ky;
			sz=b*sqrt(one-(x2/kx2)-(y2/ky2));
			n=0;
			m=zero;
			Mm=zero;  }
		else if (x==zero)  {
			sx=zero;
			sy=ay;
			sz=zero;
			n=0;
			m=zero;
			Mm=zero;  }
		else if (y==zero)  {
			sx=ax;
			sy=zero;
			sz=zero;
			n=0;
			m=zero;
			Mm=zero;  }
		else  {
			x0=x/ax;
			y0=y/ay;
			n=bisect2(x0,y0,tol,m,Mm);
			sx=cz*x/(cz+m);
			sz=zero;
			sy=y/(one+m);
		}
	}
	else  {
		if ((x==zero) && (y==zero)) {
			sx=zero;
			sy=zero;
			sz=b;
			n=0;
			m=zero;
			Mm=zero;  }
		else  {
			x0=x/ax;
			y0=y/ay;
			z0=z/b;
			n=bisect3(x0,y0,z0,tol,m,Mm);
			sx=cx*x/(cx+m);
			sy=cy*y/(cy+m);
			if ((m<zero) && ((ky2*x2 + kx2*y2) < (kx2*ky2)))  
				sz=b*sqrt(one-((sx*sx)/ax2)-((sy*sy)/ay2));
			else
				sz=z/(one+m);
		}
	}
	
//- Determination of latitude & longitude of foot point -----	

	xyz2fl(sx,sy,sz,latitude,longitude);

//- Determination of geodetic height -----	

	dx=x-sx;
	dy=y-sy;
	dz=z-sz;
	height=sqrt(dx*dx+dy*dy+dz*dz);

	if ((x+y+z)<(sx+sy+sz))   { height=-height; }

//- Adjust latitude & longitude according to signs of (xi,yi,zi)

	fl_octal(xi,yi,zi,latitude,longitude);
	
	return n;
}

void xyz2fl(long double x,long double y,long double z,long double &latitude,long double &longitude)
//  Computes the transformation of Cartesian to geodetic coordinates on the surface of the ellipsoid
//  assuming x,y,z are all non-negative
//  Angular coordinates in radians
{
	long double nom,den,dex,xme,rot;
	
	nom=mee*z;
	xme=mee*x;
	dex=xme*xme+y*y;
	den=mex*sqrt(dex);
	rot=sqrt(dex);
	
	if (den==zero)  {
		latitude=pih;
		longitude=zero;
	}
	else  {
		if (nom<=den)  
			{ latitude=atan(nom/den); }
		else
			{ latitude=pih-atan(den/nom); }
	
	    if (y<=xme)  {
			den=xme+rot;
			longitude=2.*atan(y/den);
		}
		else   {
			den=y+rot;
			longitude=pih-2.0*atan(xme/den);
		}
	}
}

void fl_octal(long double x,long double y,long double z,long double &latitude,long double &longitude)
//  Adjusts latitude & longitude according to signs of x,y,z  (angles in radians)
{
	if (z<zero)   latitude=-latitude;
	
	if (x>=zero)  {
		if (y<zero)  {longitude=-longitude;}  }
	else  {
		if (y>=zero)  {longitude=pi-longitude;}
		else  {longitude=longitude-pi;}
	}	
}

int bisect2(long double x0,long double y0,long double tol,long double &m,long double &Gm)
//  Implements the bisection method on the X-Y plane
{
	int n;
	long double d,d1,d2,Gd,g2,MC;

	n=0;
	m=-2.;
	d1=y0-one;
	g2=cz*cz*x0*x0;
	d2=sqrt(g2+y0*y0)-one;
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

int bisect3(long double x0,long double y0,long double z0,long double tol,long double &m,long double &Hm)
//  Implements the bisection method in 3D space
{
	int n;
	long double d,d1,d2,Hd,g2,g3,MC;

	n=0;
	m=-2.;
	d1=z0-one;
	g2=cx*cx*x0*x0;
	g3=cy*cy*y0*y0;
	d2=sqrt(g2+g3+z0*z0)-one;
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

int sign(long double t)
//  Returns the sign of its argument
{
	int n;
	if (t>zero) n=1;
	else if (t<zero) n=-1;
	else  n=0;
	
	return n;
}
