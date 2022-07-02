// file: src/conformal_maps.c

// simple conformal maps used to generate SLE
// Presently there are two: one introduces a titled slit 
// starting from the origin. The other introduces a vertical 
// slit starting at a point offset from the origin.
// Future: map that introduces a circular arc orthogonal to axis.

/*
    This file is part of fast_sle (version 1.0), which implements 
    a fast algorithm for simulating the Schramm-Loewner evolution (SLE)  

    Copyright (C) 2005 Tom Kennedy

    fast_sle is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    fast_sle is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with fast_sle; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    E-mail contact: tgk@math.arizona.edu
    Ordinary mail: 
        Tom Kennedy
        Mathematics Department
        University of Arizona
        Tucson, ZA 85721, USA
*/

// system includes
#include <cstdio> 
using namespace std; 
#include "math.h"
#include <iostream>
using namespace std; 
#include <cstdlib>
#include "sys/time.h"
#include "time.h"
#include "limits.h"
using namespace std;
#include <complex>
#define I complex<double>(0.,1.);

////////////////////////////////////////////////////////////////////////
// Maps that take half plane to half plane minus a curve 
// These are zipping maps
////////////////////////////////////////////////////////////////////////

complex<double> vertical_slit_ch_zip(double dt,double dx,complex<double> z)
// map_choice=1
// applies the conformal map z -> i*sqrt(-z^2+4*dt)+dx
{
    if (z.imag()<-1.e-10) {
      printf("err vertical_slit_ch_zip: z is below half plane %le\n",z.imag());
      exit(1);
    }
    complex<double> w,zz;
    w=-z*z+4*dt;
    zz=sqrt(w);
    zz*=I;
    zz+=dx;
    return(zz);
} 


complex<double> tilted_slit_ch_zip(
    double alpha, double xl, double xr,complex<double> z)
// map_choice=2
// Applies the conformal map that produces a slit starting at the origin
// with an angle of alpha*pi from the positive real axis:
//    f(z)=(z+xl)^(1-alpha)*(z-xr)^alpha
//  
// The interval [-xl,xr] gets mapped onto the slit. The length of this
// interval determines the length of the slit. Shifting this interval 
// (relative to 0) does not change this length. To map 0 to the tip of the 
// slit, we should take xr/xl=alpha/(1-alpha).
// To get a slit with capacity 2*dt we take
//   xl=2*sqrt(dt*(1-alpha)/alpha);
//   xr=2*sqrt(dt*alpha/(1-alpha));
{
    double x,y,r1,r2,theta1,theta2,r,theta;
    complex<double> zz;

    if (z.imag()<-1.e-10) {
        printf("err in tilted_slit_cmap: argument %le is below half plane",
            z.imag());
        exit(1);
    }

    // z+xl=R1*exp(i*theta1), with r1=R1^(1-alpha)
    x=z.real()+xl;
    y=z.imag();
    if (y<1.e-13) y=0.;  // kludge to avoid roundoff error
    r1=exp(log(x*x+y*y)*(1-alpha)/2.); // faster than pow()
    theta1=atan2(y,x); 

    // z-xr=R2*exp(i*theta2), with r2=R2^alpha
    x=z.real()-xr;
    y=z.imag();
    if (y<1.e-13) y=0.;  // kludge to avoid roundoff error
    r2=exp(log(x*x+y*y)*alpha/2.); // faster than pow()
    theta2=atan2(y,x); 

    r=r1*r2;
    theta=theta1*(1-alpha)+theta2*alpha;
    x=r*cos(theta);
    y=r*sin(theta);

    if (y<-1.e-10) {
        printf("err in tilted_slit_cmap: result is below half plane\n");
        exit(1);
    }

    zz=x+y*I;
    return zz;
} 

complex<double> arc_ch_zip(double a, double b, complex<double> z)
// map_choice=3
{
    complex<double> sqroot,ztemp;
    sqroot=sqrt(b*b-2*a*z-z*z);
    sqroot*=I;
    ztemp=sqroot/(b*b-a*z+a*sqroot);
    ztemp*=(a*a+b*b);
    return ztemp;
}

////////////////////////////////////////////////////////////////////////
// Real versions of zipping maps
////////////////////////////////////////////////////////////////////////

double tilted_slit_ch_zip_real(double alpha, double xl, double xr,double x)
// map_choice=2
// Applies the conformal map that produces a slit starting at the origin
// with an angle of alpha*pi from the positive real axis:
//    f(z)=(z+xl)^(1-alpha)*(z-xr)^alpha
// to a point x on the real axis outside the interval (-xl,xr). 
// So image is real.
// If x >= xr, f(x)=
// If x <= -xl, f(x)=(-x-xl)^(1-alpha)*(-x+xr)^alpha
{
    if (x >0.) return(pow(x+xl,1.-alpha)*pow(x-xr,alpha));
    if (x <0.) return(-pow(-x-xl,1.-alpha)*pow(-x+xr,alpha));
    printf("bad x in tilted_slit_ch_zip_real() %le \n",x);
    exit(1);
} 

double arc_ch_unzip_real(double a, double b, double xm)
// NOT IMPLEMENTED YET
{
    return(0.);
}


////////////////////////////////////////////////////////////////////////
// Maps that take half plane minus a curve to  half plane 
// These are unzipping maps
////////////////////////////////////////////////////////////////////////

complex<double> vertical_slit_ch_unzip(double dt,double dx,complex<double> z)
// map_choice=1
// applies the conformal map z -> i*sqrt(-(z-dx)^2-4*dt)
{
    if (z.imag()<-1.e-10) {
        printf("err vertical_slit_ch_unzip:z below half plane %le\n",z.imag());
        exit(1);
    }

    complex<double> zz,w;
    w=-(z-dx)*(z-dx)-4*dt;
    zz=sqrt(w);
    zz*=I; 
    return zz; 
} 



complex<double> tilted_slit_ch_unzip(
    double alpha, double xl, double xr,complex<double> z,complex<double> iz)
// Inverse of   f(z)=(z+xl)^(1-alpha)*(z-xr)^alpha
// This takes half plane minus slit onto half plane
// xl,xr,alpha must be such that tip goes to origin: xr/xl=alpha/(1-alpha)
// This cannot be done explicitly, so we use a Newton's method.
// We follow Marshal and Rohde and divide into four regions.
// Let L be the length of the slit. 
// 1. big region: |z| >= (9/8)*L 
// 2. tip region: |z-tip|<(1/4)*Im(tip)
// 3. right region: arg(z)<alpha*pi
// 4. left region: arg(z)>alpha*pi
// In region 1 we just use Newton. In regions 2,3,4 we consider instead of 
// f(z0)=z, the equation k(f(z0))=k(z) where 
// 2. k(z)=sqrt(z/ztip-1) 
// 3. k(z)=z^{1/alpha}
// 4. k(z)=z^{1/(1-alpha)}
// We scale out len at the start and put it back it at the end.
{
    complex<double> ztip,theta,fz,fpz,z0,z1,temp,zt;
    long region;

    if (fabs(xr/xl-alpha/(1.-alpha))>1.e-8) {
        printf("bad xr,xl,alpha in   tilted_slit_ch_unzip \n");
        exit(0);
    }

    double len=pow(xl,1.-alpha)*pow(xr,alpha);
    xl/=len;
    xr/=len;
    z/=len;

    // kludges:
    // map origin to -xl
    if (abs(z) < 1.e-8) {
        z0=-xl;
        return(z0*len);
    }
    // map tip to origin 
    if (abs(z-ztip) < 1.e-8) {
        z0=0.;
        return(z0*len);
    }

    if (abs(z)>1.e4) {
        // Newton method will fail for large z due to numerical error.
        // Use expansion about infinity to compute inverse of f(z):
        // f^{-1}(z)=z-ddrive+dcap/z + ...
        // dcap is 2*t, dx=calpha*sqrt(t)
        // ddrive=calpha*sqrt(dcap/2.)
        double dcap=0.5*pow(alpha,1-2*alpha)*pow(1-alpha,2*alpha-1);
        double calpha=2*(1-2*alpha)/sqrt(alpha*(1-alpha));  
        double ddrive=calpha*sqrt(dcap/2.);
        z0=z-ddrive+dcap/z;
        return(z0*len);
    }

    if (abs(z)>1.125) { // region 1 : large z
        region=1;
        z0=z;
        zt=z;
        for (long i=0;i<40;i++) { // Newton iterations
            // compute value of function and its derivative
            fz=tilted_slit_ch_zip(alpha,xl,xr,z0);
            fpz=fz*((1.-alpha)/(z0+xl)+alpha/(z0-xr));
            if (abs(fz-z)<1.e-10) return(z0*len);
            z0=z0-(fz-zt)/fpz; 
            // if z0 is out of half plane, push it back in
            if (z0.imag()<0.) z0=z0.real();
        } 
        printf("Newton failed to converge in region 1 %le\n",abs(z));
    }

    theta=I;
    theta*=M_PI*alpha;
    ztip=exp(theta); // note that |ztip|=1 because of the scaling

    for (long itry=1;itry<=4;itry++) {

        if (itry>1) printf(">>>> itry==%ld, region was %ld \n",itry,region);

        // determine region (we are not in region 1 if we got here)
        if (abs(z-ztip)<0.25*ztip.imag()) region=2;
        else {
            if (abs(z) < 1.e-6 || arg(z)>alpha*M_PI) region=4;
            else region=3;
        }
  
        // compute initial guess and target
        switch (region) {
            case 2: 
                switch(itry) {
                    case 1: z0=0.+0.05*I; break;
                    case 2: z0=0.+0.5*I; break;
                    case 3: z0=0.+2.*I; break;
                    case 4: z0=0.+0.1*I; break;
                }
                zt=sqrt(z/ztip-1.);
            break;
            case 3: 
                switch(itry) {
                    case 1: z0=0.5*(ztip+xr); break; 
                    case 2: z0=0.+0.05*I; break;
                    case 3: z0=0.5*xr; break;
                    case 4: 
                        z0=exp(log(z)/alpha-(1-alpha)*log(xl+xr)/alpha)+xr;
                    break;
                }
                zt=exp(log(z)/alpha);
            break;
            case 4: 
                switch(itry) {
                    case 1: z0=0.5*(ztip-xl); break; 
                    case 2: z0=0.+0.05*I; break;
                    case 3: z0=-0.5*xl; break;
                    case 4: z0=-1.01*xl; break;
                }
                zt=exp(log(z)/(1.-alpha));
            break;
        }

        // Newton iterations
    
        for (long i=0;i<40;i++) {
            // compute value of function and its derivative
            switch (region) {
                case 2: 
                    temp=tilted_slit_ch_zip(alpha,xl,xr,z0);
                    fz=sqrt(temp/ztip-1.);
                    fpz=temp*((1.-alpha)/(z0+xl)+alpha/(z0-xr))
                       /(2.*fz*ztip);
                break;
                case 3: 
                    temp=tilted_slit_ch_zip(alpha,xl,xr,z0);
                    fz=exp(log(temp)/alpha);
                    fpz=fz*((1.-alpha)/(z0+xl)+alpha/(z0-xr))/alpha;
                break;
                case 4: 
                    temp=tilted_slit_ch_zip(alpha,xl,xr,z0);
                    fz=exp(log(temp)/(1.-alpha));
                    fpz=fz*((1.-alpha)/(z0+xl)+alpha/(z0-xr))/(1.-alpha);
                break;
            }
            if (abs(temp-z)<1.e-10) return(z0*len);

            z0=z0-(fz-zt)/fpz; 
            // if z0 is out of half plane, push it back in
            if (z0.imag()<0.) z0=z0.real();
        } 
    }


    if (alpha<0.5) { 
        z0=exp(log(z)/alpha-(1-alpha)*log(xl+xr)/alpha)+xr;
        temp=exp((1-alpha)*log(z0+xl)+alpha*log(z0-xr));
        printf("!! temp=%lf %lf %le\n",temp.real(),temp.imag(),abs(temp-z));
        if (abs(temp-z)<1.e-4) return(z0*len);
    }
    else { 
        complex<double> ii;
        ii=I;
	z0=exp(log(z)/(1-alpha)-alpha*(log(xl+xr)+M_PI*ii)/(1-alpha))-xl;
        temp=exp((1-alpha)*log(z0+xl)+alpha*log(z0-xr));
        printf("!! temp=%lf %lf %le\n",temp.real(),temp.imag(),abs(temp-z));
        if (abs(temp-z)<1.e-4) return(z0*len);

    }

    printf("region=%ld\n",region);
    printf("alpha=%lf xl=%lf xr=%lf\n",alpha,xl,xr);
    printf("z=%lf + %lf i, arg=%lf\n",z.real(),z.imag(),arg(z));
    printf("ztip=%lf + %lf i, arg=%lf\n",ztip.real(),ztip.imag(),arg(ztip));
    printf("??? Newton failed to converge, regions 2,3 or 4 \n");
    exit(1);
}


complex<double> arc_ch_unzip(double a, double b, complex<double> z)
{
// conformal map that takes half plane minus an arc onto the half plane. 
// arc is from origin to a+ib, orthogonal to horizontal axis
// see cap.tex for derivation of the map 
    complex<double> sqroot,ztemp,ztemp2,ztemp3;
    double caprsq;
    caprsq=a*a+b*b;
    sqroot=sqrt(-(z*z-2.*a*z+caprsq)/((caprsq-a*z)*(caprsq-a*z)));
    sqroot*= I;
    ztemp=b*b*sqroot/(a*sqroot+1.);

    // watch out for z on the real axis, use sqrt(caprsq) as a scale 
    // Actually, z on the arc is a problem too, but this is unresolvable
    if (z.imag()*z.imag()>1.e-12*caprsq) return ztemp;
    sqroot*=-1.;
    ztemp2=b*b*sqroot/(a*sqroot+1.);
    // to choose between ztemp and ztemp2 we perturb z into the upper half
    // We perturb 0 so it is left of arc. 
    ztemp3=I;
    ztemp3-=1.;
    ztemp3*=0.01*sqrt(caprsq);
    z+=ztemp3;
    sqroot=sqrt(-(z*z-2.*a*z+caprsq)/((caprsq-a*z)*(caprsq-a*z)));
    sqroot*= I;
    ztemp3=b*b*sqroot/(a*sqroot+1.);
    if (abs(ztemp3-ztemp)<abs(ztemp3-ztemp2)) return ztemp;
    else return ztemp2;
} 



////////////////////////////////////////////////////////////////////////
// Real versions of unzipping maps
////////////////////////////////////////////////////////////////////////

double tilted_slit_ch_unzip_real(double alpha, double xl, double xr,
    double x)
{
// use simple bisection instead of newton
    double xlow,xhigh,xmid,fx;

    if (x>0.) {
        xlow=xr;
        xhigh=x+xr;
        while (tilted_slit_ch_zip_real(alpha,xl,xr,xhigh)<x) xhigh*=2.;
    }
    else {
        xhigh=-xl;
        xlow=x-xl;
        while (tilted_slit_ch_zip_real(alpha,xl,xr,xlow)>x) xlow*=2.;
    }

    while (xhigh-xlow> 1.e-10*fabs(x)) {
        xmid=(xlow+xhigh)/2.;
        fx=tilted_slit_ch_zip_real(alpha,xl,xr,xmid);
        if (fx<x) xlow=xmid;
        else xhigh=xmid;
    }

    return((xlow+xhigh)/2.);
}

void tilted_slit_ch_unzip_origin(double alpha, double xl, double xr,
			       double *xm, double *xp)
// returns two images of origin under the map
{
    *xm= -xl;
    *xp= xr;
}

void arc_ch_unzip_origin(double a, double b, double *xm, double *xp)
// NOT IMPLEMENTED YET
{
    *xm=0.;
    *xp=0.;
}


////////////////////////////////////////////////////////////////////////
// Maps that take disc to disc minus a curve. These are zipping maps
////////////////////////////////////////////////////////////////////////

complex<double> vertical_slit_rad_zip(
    double dt,double dx,complex<double> z)
// map_choice=1
{
    complex<double> zz;

    if (abs(z)>1+1.e-10) {
        printf("err in complex vertical_slit_rad_zip: z out of disc\n");
        exit(1);
    }

    // map to right half plane
    z=(1.-z)/(1.+z);

    // insert horizontal slit from 1 to left
    z=sqrt(exp(-dt)*z*z+1.-exp(-dt));

    // map to disc
    z=(1.-z)/(1.+z);

    // zz=exp(-i*dx) z
    zz=cos(dx)+sin(dx)*I;
    z*=zz;

    return z;
} 

complex<double> tilted_slit_rad_zip(double dx,
    double alpha, double xl, double xr, complex<double> ww, complex<double> z)
  // this map fixes 0 by : 0 -> i -> ww -> 0
{
    complex<double> ii,wwbar;
    ii=I;

    // map unit disc to half plane, 0->i, 1->0: i(1-z)/(1+z) 
    z=ii*(1.-z)/(1.+z);

    // map half plane to half plane minus slit
    z=tilted_slit_ch_zip(alpha,xl,xr,z);

    // map half plane back to disc, ww->0, 0->1
    // ww is image of i under above map; was computed in compute_parms_rad_zip()
    wwbar=ww.real()-ww.imag()*I;
    z=wwbar*(z-ww)/(ww*(z-wwbar));
    return(z);
}

////////////////////////////////////////////////////////////////////////
// Maps for radial case in half plane
////////////////////////////////////////////////////////////////////////

complex<double> vertical_slit_half_rad_zip(
    double dt,double dx,complex<double> z)
// map_choice=1
{
    complex<double> zz,c,ii;

    ii=I;

    if (z.imag()<-1.e-10) {
        printf("err in vertical_slit_half_rad_zip: z out of half plane\n");
        printf("z=%lf + %lf i \n",z.real(),z.imag());
        exit(1);
    }

    // rotate upper half plane to right half plane 
    z=-ii*z;

    // map to disc
    zz=(1.-z)/(1.+z);

    // rotate: zz=exp(-i*dx) z
    c=cos(dx)-sin(dx)*ii;
    zz*=c;

    // map to right half plane
    zz=(1.-zz)/(1.+zz);

    // insert horizontal slit from 1 to left
    zz=sqrt(exp(-dt)*zz*zz+1.-exp(-dt));

    // rotate right half plane to upper half plane 
    zz=ii*zz;

    return zz;
} 

complex<double> tilted_slit_half_rad_zip(double dx,
    double alpha, double xl, double xr, complex<double> ww, complex<double> z)
{
    printf("tilted_slit_half_rad_zip not implemented \n");
    exit(10);
    return(0.);
}

////////////////////////////////////////////////////////////////////////
// Maps that take disc to other domains
////////////////////////////////////////////////////////////////////////

complex<double> disc_to_half(complex<double> z,complex<double> z1)
// takes unit disc to upper half plane, 1 to 0, sending 0 to z1 
{
    z=(z-1.)*abs(z1)*abs(z1)/(z1*z-conj(z1));
    return(z);
}

////////////////////////////////////////////////////////////////////////
// Maps that take half plane to other domains 
////////////////////////////////////////////////////////////////////////

complex<double> half_to_strip(complex<double> z,complex<double> z1)
// maps half plane to strip: 0<Im(z)<1, sending 0 -> 0, infinity -> z1
// NB only z1=i is implemented so far
{
    complex<double> w;
    if (fabs(z1.real())+fabs(z1.imag()-1.)>1.e-6) {
        printf("half_to_strip only implement for infinity to i\n");
        exit(10);
    }
    w=(1.+z)/(1.-z);
    w=log(w)/M_PI;
    return(w);
}

////////////////////////////////////////////////////////////////////////
// Maps that take strip to strip minus a curve. These are zipping maps
////////////////////////////////////////////////////////////////////////

complex<double> vertical_slit_di_zip(double dt,double dx,complex<double> z)
{
    if (z.imag()<-1.e-10 || z.imag()>M_PI+1.e-10) {
        printf("err in complex vertical_slit_id_zip: z is out of strip\n");
        exit(1);
    }

    complex<double> zz,cz,root;
    cz=(exp(z/2.)+exp(-z/2.))/2.;
    root=sqrt(exp(-dt)*cz*cz-1.);
    zz=2.*log(exp(-dt/2.)*cz+root)+dx;
    if (zz.imag()>0.) return(zz);
    zz=2.*log(exp(-dt/2.)*cz-root)+dx;
    return zz;
}


///////////////////////////////////////////////////////////////
// SC map and its inverse for equilateral triangle and square
///////////////////////////////////////////////////////////////
// 
// rotation maps for triangle and square seem to be in different directions
// 
// 
// 
// 
// 



complex<double> assign(double x, double y)
{
    complex<double> w;
    w=I;
    w=x+y*w;
    return(w);
}

complex<double> sc_triangle_integrand(complex<double> z)
{
  complex<double> w;
  w=pow(z-1.,2./3.)*pow(z+1.,2./3.);
  w=1./w;
  return(w);
}

complex<double> sc_square_integrand(complex<double> z)
{
  complex<double> w;
  w=sqrt(z-1.)*sqrt(z)*sqrt(z+1.);
  w=1./w;
  return(w);
}

complex<double> scd_square_integrand(complex<double> z)
{
  complex<double> w;
  w=sqrt(1.-z*z*z*z);
  w=1./w;
  return(w);
}

complex<double> sc_triangle_rotate_120_hplane(complex<double> z)
{
    complex<double> w;
    w=(z+3.)/(1.-z);
    return(w); 
}

complex<double> sc_triangle_rotate_60_hplane(complex<double> z)
{
    complex<double> w;
    w=(3.*z+3.)/(3.-z);
    return(w); 
}

complex<double> sc_square_rotate_90_hplane(complex<double> z)
{
    complex<double> w;
    w=(z+1.)/(1.-z);
    return(w); 
}

complex<double> sc_square_rotate_45_hplane(complex<double> z)
{
    complex<double> w;
    w=(z+sqrt(2.)-1.)/(sqrt(2)+1.-z);
    w*=sqrt(2.)+1;
    return(w); 
}

complex<double> sc_triangle_integrate_ps(complex<double> z)
{
  complex<double> w,total,wpow,wsq;
  double c;
  long k;
  complex<double> prefactor;
  double pre_mag;

  // 1.2143... is approx int from 0 to sqrt(3) of (x^2+1)^{-2/3} 
  pre_mag=(1./3.)/1.2143253239438;
  prefactor=assign(pre_mag*cos(2*M_PI/3.),pre_mag*sin(2*M_PI/3.));

  w=1./z;
  wpow=pow(w,1./3.);
  wsq=w*w;
  total=0.;
  // c is coef of power series of (1-x)^{-2/3}
  for (k=0;k<=100;k++) {
    if (k==0) c=1.;
    if (k==1) c=2./3.;
    //    if (k>1) c=c*(k-1./3.)/k;
    if (k>1) c=c*(1-1./(3.*k));
    //    total+=c*pow(w,2*k+1./3.)/(2*k+1./3.);
    //    total+=c*wthird*pow(w,2*k)/(2*k+1./3.);
    if (abs(wpow)<1.e-15) break;
    total+=c*wpow/(2*k+1./3.);
    wpow*=wsq;
  }
  total*=prefactor;
  total=-total;
  total+=assign(0.,2./3.);
  return(total);
}

complex<double> sc_triangle_integrate_ps_test(complex<double> z)
{
  complex<double> w,total,wpow;
  double c;
  long k;
  complex<double> prefactor;
  double pre_mag;

  // 1.2143... is approx int from 0 to sqrt(3) of (x^2+1)^{-2/3} 
  pre_mag=(1./3.)/1.2143253239438;
  prefactor=assign(pre_mag*cos(2*M_PI/3.),pre_mag*sin(2*M_PI/3.));

  w=1./z;
  wpow=pow(w,1./3.);
  total=0.;
  // c is coef of power series of (1-x)^{-2/3}
  for (k=0;k<=100;k++) {
    if (k==0) c=1.;
    if (k==1) c=2./3.;
    if (k>1) c=c*(k-1./3.)/k;
    //    total+=c*pow(w,2*k+1./3.)/(2*k+1./3.);
    total+=c*wpow/(2*k+1./3.);
    wpow*=(w*w);
  }
  total*=prefactor;
  total=-total;
  total+=assign(0.,2./3.);
  return(total);
}

complex<double> sc_square_integrate_ps(complex<double> z)
{
  complex<double> w,total,wpow,wsq;
  double c;
  long k;
  w=1./z;
  wpow=pow(w,1./2.);
  wsq=w*w;
  total=0.;
  complex<double> prefactor;

  // prefactor is related to integral of 1/sqrt(x^4+1) from 0 to 1
  prefactor=assign(0.,1./(0.927037338651*sqrt(2.)));

  // c is coef of power series of (1-x)^{-1/2}
  for (k=0;k<=100;k++) {
    if (k==0) c=1.;
    if (k==1) c=1./2.;
    if (k>1) c=c*(k-1./2.)/k;
    //    total+=c*pow(w,2*k+1./2.)/(2*k+1./2.);
    if (abs(wpow)<1.e-15) break;
    total+=c*wpow/(2*k+1./2.);
    wpow*=wsq;
  }


  total*=prefactor;
  total=-total;
  total+=assign(1.,1.);
  return(total);
}

complex<double> scd_square_integrate_ps(complex<double> z)
{
  complex<double> total,wpow,wfourth;
  double c;
  long k;
  complex<double> scd_prefactor;

  wpow=z;
  wfourth=z*z*z*z;
  total=0.;

  // 0.927... is integral of 1/sqrt(x^4+1) from 0 to 1
  scd_prefactor=1./0.927037338651;
  scd_prefactor*=assign(1./sqrt(2.),1./sqrt(2.));

  // c is coef of power series of (1-x)^{-1/2}
  for (k=0;k<=100;k++) {
    if (k==0) c=1.;
    if (k==1) c=1./2.;
    if (k>1) c=c*(k-1./2.)/k;
    if (abs(wpow)<1.e-15) break;
    total+=c*wpow/(4*k+1.);
    wpow*=wfourth;
  }

  total*=scd_prefactor;
  return(total);
}

complex<double> sc_square_integrate_ps_test(complex<double> z)
{
  complex<double> w,total;
  double c;
  long k;
  w=1./z;
  total=0.;
  complex<double> prefactor;

  // prefactor is related to integral of 1/sqrt(x^4+1) from 0 to 1
  prefactor=assign(0.,1./(0.927037338651*sqrt(2.)));

  // c is coef of power series of (1-x)^{-1/2}
  for (k=0;k<=100;k++) {
    if (k==0) c=1.;
    if (k==1) c=1./2.;
    if (k>1) c=c*(k-1./2.)/k;
    total+=c*pow(w,2*k+1./2.)/(2*k+1./2.);
  }
  total*=prefactor;
  total=-total;
  total+=assign(1.,1.);

  return(total);
}

complex<double> sc_triangle_integrate(complex<double> z)
// integrate from isqrt(3) to z
{
    complex<double> total,w,z0;
    double coef;
    long i;
    long n;
    complex<double> prefactor;
    double pre_mag;

    // 1.2143... is approx int from 0 to sqrt(3) of (x^2+1)^{-2/3} 
    pre_mag=(1./3.)/1.2143253239438;
    prefactor=assign(pre_mag*cos(2*M_PI/3.),pre_mag*sin(2*M_PI/3.));

    n=10000;

    z0=assign(0.,sqrt(3.));
    total=0.;

    for (i=0;i<=n;i++) {
        coef=1.;
        if (i==0 || i==n) coef=3./8.;        
        if (i==1 || i==n-1) coef=7./6.;      
        if (i==2 || i==n-2) coef=23./24.;    
        w=z0+(z-z0)*double(i)/double(n);
        total+=coef*sc_triangle_integrand(w);
    } 
    total*=(z-z0)/double(n);
    total*=prefactor;
    return(total);
}

complex<double> sc_square_integrate(complex<double> z)
// integrate from i to z
{
    complex<double> total,w,z0;
    double coef;
    long i,n;

    complex<double> prefactor;
    // prefactor is related to integral of 1/sqrt(x^4+1) from 0 to 1
    prefactor=assign(0.,1./(0.927037338651*sqrt(2.)));

    z0=assign(0.,1.);

    n=1000; // qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq

    total=0.;

    for (i=0;i<=n;i++) {
        coef=1.;
        if (i==0 || i==n) coef=3./8.;        
        if (i==1 || i==n-1) coef=7./6.;      
        if (i==2 || i==n-2) coef=23./24.;    
        w=z0+(z-z0)*double(i)/double(n);
        total+=coef*sc_square_integrand(w);
    } 
    total*=(z-z0)/double(n);
    total*=prefactor;
    return(total);
}

complex<double> scd_square_integrate(complex<double> z)
// integrate from 0 to z
{
    complex<double> total,w,z0;
    double coef;
    long i,n;

    complex<double> scd_prefactor;
    // 0.927... is integral of 1/sqrt(x^4+1) from 0 to 1
    scd_prefactor=1./0.927037338651;
    scd_prefactor*=assign(1./sqrt(2.),1./sqrt(2.));

    z0=assign(0.,0.);

    n=1000; // qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq

    total=0.;

    for (i=0;i<=n;i++) {
        coef=1.;
        if (i==0 || i==n) coef=3./8.;        
        if (i==1 || i==n-1) coef=7./6.;      
        if (i==2 || i==n-2) coef=23./24.;    
        w=z0+(z-z0)*double(i)/double(n);
        total+=coef*scd_square_integrand(w);
    } 
    total*=(z-z0)/double(n);
    total*=scd_prefactor;
    return(total);
}

complex<double> sc_square_integrate_test(complex<double> z)
// integrate from i to z
{
    complex<double> total,w,z0;
    double coef;
    long i,n;

    complex<double> prefactor;
    // prefactor is related to integral of 1/sqrt(x^4+1) from 0 to 1
    prefactor=assign(0.,1./(0.927037338651*sqrt(2.)));

    z0=assign(0.,1.);

    n=10000;

    total=0.;

    for (i=0;i<=n;i++) {
        coef=1.;
        if (i==0 || i==n) coef=3./8.;        
        if (i==1 || i==n-1) coef=7./6.;      
        if (i==2 || i==n-2) coef=23./24.;    
        w=z0+(z-z0)*double(i)/double(n);
        total+=coef*sc_square_integrand(w);
    } 
    total*=(z-z0)/double(n);
    total*=prefactor;
    return(total);
}

complex<double> sc_triangle_map(complex<double> z)
{
    complex<double> w0,w1,w2,w,fw,rot;
    long nrotate;     
    complex<double> prefactor;
    double pre_mag;
    // 1.2143... is approx int from 0 to sqrt(3) of (x^2+1)^{-2/3} 
    //    pre_mag=(1./3.)/1.2143253239438;
    //    prefactor=assign(pre_mag*cos(2*M_PI/3.),pre_mag*sin(2*M_PI/3.));

    w0=z;
    w1=sc_triangle_rotate_120_hplane(w0);
    w2=sc_triangle_rotate_120_hplane(w1);

    if (abs(w0)>=abs(w1) && abs(w0)>=abs(w2)) {
      nrotate=0;
      w=w0;
    }
    if (abs(w1)>=abs(w0) && abs(w1)>=abs(w2)) {
      nrotate=1;
      w=w1;
    }
    if (abs(w2)>=abs(w0) && abs(w2)>=abs(w1)) {
      nrotate=2;
      w=w2;
    }

    if (abs(w)<(sqrt(3.)-0.001)) {
      printf("integrate called: abs(w)=%lf < %lf, %lf %lf %lf \n",
	     abs(w),sqrt(3.)-0.1,abs(w0),abs(w1),abs(w2));
      fw=sc_triangle_integrate(w);
    }
    else {
      fw=sc_triangle_integrate_ps(w);
    }

    rot=assign(cos(-nrotate*2*M_PI/3),sin(-nrotate*2*M_PI/3));
    fw=fw*rot;
    return(fw);
}

complex<double> sc_triangle_map_test(complex<double> z)
{
    complex<double> w0,w1,w2,w,fw,rot;
    long nrotate;     
    complex<double> prefactor;
    double pre_mag;
    // 1.2143... is approx int from 0 to sqrt(3) of (x^2+1)^{-2/3} 
    pre_mag=(1./3.)/1.2143253239438;
    prefactor=assign(pre_mag*cos(2*M_PI/3.),pre_mag*sin(2*M_PI/3.));

    w0=z;
    w1=sc_triangle_rotate_120_hplane(w0);
    w2=sc_triangle_rotate_120_hplane(w1);

    if (abs(w0)>=abs(w1) && abs(w0)>=abs(w2)) {
      nrotate=0;
      w=w0;
    }
    if (abs(w1)>=abs(w0) && abs(w1)>=abs(w2)) {
      nrotate=1;
      w=w1;
    }
    if (abs(w2)>=abs(w0) && abs(w2)>=abs(w1)) {
      nrotate=2;
      w=w2;
    }

    if (abs(w)<(sqrt(3.)-0.001)) {
      printf("integrate called: abs(w)=%lf < %lf, %lf %lf %lf \n",
	     abs(w),sqrt(3.)-0.1,abs(w0),abs(w1),abs(w2));
      fw=sc_triangle_integrate(w);
    }
    else {
      fw=sc_triangle_integrate_ps_test(w);
    }

    rot=assign(cos(-nrotate*2*M_PI/3),sin(-nrotate*2*M_PI/3));
    fw=fw*rot;
    return(fw);
}

complex<double> sc_triangle_map_inv(complex<double> z)
// computes inverse of sc_tri_map() by Newton's method
{
    complex<double> w0,w,fw0,fprime,dif,cent,zz;
    long i,iguess,nguess,nrotate;
    double theta;    
    complex<double> prefactor;
    double pre_mag;
    // 1.2143... is approx int from 0 to sqrt(3) of (x^2+1)^{-2/3} 
    pre_mag=(1./3.)/1.2143253239438;
    prefactor=assign(pre_mag*cos(2*M_PI/3.),pre_mag*sin(2*M_PI/3.));

    // rotate point 
    theta=180.*arg(z)/M_PI;    
    if (theta<0.) theta+=360.;
    nrotate=0;
    if (theta>=150. && theta<270.) {
      z*=assign(cos(-2*M_PI/3),sin(-2*M_PI/3));
      nrotate=1;
    }
    else if (theta>=270. || theta<30.) {
      z*=assign(cos(-4*M_PI/3),sin(-4*M_PI/3));
      nrotate=2;
    }    

    theta=180.*arg(z)/M_PI;    
    if (theta<29.999999999 || theta>150.0000001) {
      printf("ERR theta=%lf \n",theta);
      exit(10);
    }

    nguess=5;
    for (iguess=1;iguess<=nguess;iguess++) {
      // initial guess 
      switch(iguess) {
        case 1: w0=assign(0.,sqrt(3.));  break;  
        case 2: w0=assign(-2.,1.); break;  
        case 3: w0=assign(2.,1.); break;  
        case 4: w0=assign(-3.,0.5); break;  
        case 5: w0=assign(3.,0.5); break;  
      }

      //if (iguess>1) printf("iguess=%ld \n",iguess);

      w=w0;

      for (i=0;i<20;i++) {
        fw0=sc_triangle_map(w0);
        if (isnan(abs(fw0)) || abs(fw0)>1.e10) break;
	//   printf("i=%ld, error=%le \n",i,abs(fw0-z));
        if (abs(fw0-z)<1.e-11) {
          for (long irotate=0;irotate<nrotate;irotate++) 
            w=sc_triangle_rotate_120_hplane(w);
          return(w); 
        }
        fprime=prefactor*sc_triangle_integrand(w0);
        dif=(z-fw0)/fprime;
        w=w0+dif;
        // if newton leaves half plane, push it back in 
        if (imag(w)<0.) w=real(w);
        if (isnan(abs(w)) || abs(w)>1.e10) break;
        w0=w;
      }
    }

    printf("FAILURE sc_triangle_map_inv(), z=%lf+%lfi\n",real(z),imag(z));
    return(acos(2.));
}

complex<double> sc_triangle_map_inv_test(complex<double> z)
// computes inverse of sc_tri_map() by Newton's method
{
    complex<double> w0,w,fw0,fprime,dif,cent,zz;
    long i,iguess,nguess,nrotate;
    double theta;    
    complex<double> prefactor;
    double pre_mag;
    // 1.2143... is approx int from 0 to sqrt(3) of (x^2+1)^{-2/3} 
    pre_mag=(1./3.)/1.2143253239438;
    prefactor=assign(pre_mag*cos(2*M_PI/3.),pre_mag*sin(2*M_PI/3.));

    // rotate point 
    theta=180.*arg(z)/M_PI;    
    if (theta<0.) theta+=360.;
    nrotate=0;
    if (theta>=150. && theta<270.) {
      z*=assign(cos(-2*M_PI/3),sin(-2*M_PI/3));
      nrotate=1;
    }
    else if (theta>=270. || theta<30.) {
      z*=assign(cos(-4*M_PI/3),sin(-4*M_PI/3));
      nrotate=2;
    }    

    theta=180.*arg(z)/M_PI;    
    if (theta<29.999999999 || theta>150.0000001) {
      printf("ERR theta=%lf \n",theta);
      exit(10);
    }

    nguess=5;
    for (iguess=1;iguess<=nguess;iguess++) {
      // initial guess 
      switch(iguess) {
        case 1: w0=assign(0.,sqrt(3.));  break;  
        case 2: w0=assign(-2.,1.); break;  
        case 3: w0=assign(2.,1.); break;  
        case 4: w0=assign(-3.,0.5); break;  
        case 5: w0=assign(3.,0.5); break;  
      }

      //      if (iguess>1) printf("iguess=%ld \n",iguess);

      w=w0;

      for (i=0;i<20;i++) {
        fw0=sc_triangle_map_test(w0);
        if (isnan(abs(fw0)) || abs(fw0)>1.e10) break;
	//   printf("i=%ld, error=%le \n",i,abs(fw0-z));
        if (abs(fw0-z)<1.e-11) {
          for (long irotate=0;irotate<nrotate;irotate++) 
            w=sc_triangle_rotate_120_hplane(w);
          return(w); 
        }
        fprime=prefactor*sc_triangle_integrand(w0);
        dif=(z-fw0)/fprime;
        w=w0+dif;
        // if newton leaves half plane, push it back in 
        if (imag(w)<0.) w=real(w);
        if (isnan(abs(w)) || abs(w)>1.e10) break;
        w0=w;
      }
    }

    printf("FAILURE sc_triangle_map_inv(), z=%lf+%lfi\n",real(z),imag(z));
    return(acos(2.));
}

complex<double> sc_square_map(complex<double> z)
{
    complex<double> w0,w1,w2,w3,w,fw,rot,ii,zd,fz;
    long nrotate;     
    complex<double> prefactor;
 
    // prefactor is related to integral of 1/sqrt(x^4+1) from 0 to 1
    prefactor=assign(0.,1./(0.927037338651*sqrt(2.)));

    w0=z;
    w1=sc_square_rotate_90_hplane(w0);
    w2=sc_square_rotate_90_hplane(w1);
    w3=sc_square_rotate_90_hplane(w2);

    if (abs(w0)>=abs(w1) && abs(w0)>=abs(w2) && abs(w0)>=abs(w3)) {
      nrotate=0;
      w=w0;
    }
    if (abs(w1)>=abs(w0) && abs(w1)>=abs(w2) && abs(w1)>=abs(w3)) {
      nrotate=1;
      w=w1;
    }
    if (abs(w2)>=abs(w0) && abs(w2)>=abs(w1) && abs(w2)>=abs(w3)) {
      nrotate=2;
      w=w2;
    }
    if (abs(w3)>=abs(w0) && abs(w3)>=abs(w1) && abs(w3)>=abs(w2)) {
      nrotate=3;
      w=w3;
    }

    if (abs(w)<2.) {
      ii=assign(0.,1.);
      zd=(w-ii)/(w+ii);
      fw=scd_square_integrate_ps(zd);
    }
    else fw=sc_square_integrate_ps(w);

    rot=assign(cos(-nrotate*2*M_PI/4),sin(-nrotate*2*M_PI/4));
    fw=fw*rot;
    return(fw);
}

complex<double> scd_square_map(complex<double> z)
{
    complex<double> w,fw;

    w=z;
    //    if (abs(w)<2.) fw=sc_square_integrate(w);
    //    else fw=sc_square_integrate_ps(w);
    fw=scd_square_integrate(w);
    return(fw);
}

complex<double> sc_square_map_test(complex<double> z)
{
    complex<double> w0,w1,w2,w3,w,fw,rot;
    long nrotate;     
    complex<double> prefactor;
 
    // prefactor is related to integral of 1/sqrt(x^4+1) from 0 to 1
    prefactor=assign(0.,1./(0.927037338651*sqrt(2.)));

    w0=z;
    w1=sc_square_rotate_90_hplane(w0);
    w2=sc_square_rotate_90_hplane(w1);
    w3=sc_square_rotate_90_hplane(w2);

    if (abs(w0)>=abs(w1) && abs(w0)>=abs(w2) && abs(w0)>=abs(w3)) {
      nrotate=0;
      w=w0;
    }
    if (abs(w1)>=abs(w0) && abs(w1)>=abs(w2) && abs(w1)>=abs(w3)) {
      nrotate=1;
      w=w1;
    }
    if (abs(w2)>=abs(w0) && abs(w2)>=abs(w1) && abs(w2)>=abs(w3)) {
      nrotate=2;
      w=w2;
    }
    if (abs(w3)>=abs(w0) && abs(w3)>=abs(w1) && abs(w3)>=abs(w2)) {
      nrotate=3;
      w=w3;
    }

    if (abs(w)<2.) fw=sc_square_integrate_test(w);
    else fw=sc_square_integrate_ps_test(w);

    rot=assign(cos(-nrotate*2*M_PI/4),sin(-nrotate*2*M_PI/4));
    fw=fw*rot;
    return(fw);
}

complex<double> sc_square_map_inv(complex<double> z)
// computes inverse of SC map for square by Newton's method
{
    complex<double> w0,w,fw0,fprime,dif,cent,zz;
    long i,iguess,nguess,nrotate;
    double theta;
    complex<double> prefactor;
 
    // prefactor is related to integral of 1/sqrt(x^4+1) from 0 to 1
    prefactor=assign(0.,1./(0.927037338651*sqrt(2.)));

    // rotate point 
    theta=180.*arg(z)/M_PI;    
    if (theta<0.) theta+=360.;
    nrotate=0;
    if (theta>=90. && theta<180.) {
      z*=assign(cos(-M_PI/2.),sin(-M_PI/2.));
      nrotate=1;
    }
    else if (theta>=180. && theta<270.) {
      z*=assign(cos(-2*M_PI/2.),sin(-2*M_PI/2.));
      nrotate=2;
    }    
    else if (theta>=270.) {
      z*=assign(cos(-3*M_PI/2.),sin(-3*M_PI/2.));
      nrotate=3;
    }    

    theta=180.*arg(z)/M_PI;    
    if (theta>90.0000001) {
      printf("ERR theta=%lf, nrotate=%ld\n",theta,nrotate);
      exit(10);
    }

    nguess=5;
    for (iguess=1;iguess<=nguess;iguess++) {
      // initial guess 
      switch(iguess) {
        case 1: w0=assign(0.,1.); break;  
        case 2: w0=assign(-0.5,0.5); break;  
        case 3: w0=assign(0.5,0.5); break;  
        case 4: w0=assign(-1.5,0.5); break;  
        case 5: w0=assign(-1.5,0.5); break;  
      }

      //      if (iguess>1) printf("iguess=%ld \n",iguess);

      w=w0;

      for (i=0;i<20;i++) {
        fw0=sc_square_map(w0);
        if (isnan(abs(fw0)) || abs(fw0)>1.e10) break;
        //printf("i=%ld, error=%le \n",i,abs(fw0-z));
        if (abs(fw0-z)<1.e-11) {
          for (long irotate=0;irotate<nrotate;irotate++) 
            w=sc_square_rotate_90_hplane(w);
          return(w); 
        }
        fprime=prefactor*sc_square_integrand(w0);
        dif=(z-fw0)/fprime;
        w=w0+dif;
        // if newton leaves half plane, push it back in 
        if (imag(w)<0.) w=real(w);
        if (isnan(abs(w)) || abs(w)>1.e10) break;
        w0=w;
      }

    }

    printf("FAILURE sc_square_map_inv(), z=%lf+%lfi\n",real(z),imag(z));
    exit(10);
    return(acos(2.));
}

complex<double> scd_square_map_inv(complex<double> z)
// computes inverse of SC map for square by Newton's method
{
    complex<double> w0,w,fw0,fprime,dif,cent,zz;
    long i,iguess,nguess,nrotate;
    double theta;
    complex<double> scd_prefactor;
 
    // 0.927... is integral of 1/sqrt(x^4+1) from 0 to 1
    scd_prefactor=1./0.927037338651;
    scd_prefactor*=assign(1./sqrt(2.),1./sqrt(2.));

    // rotate point 
    theta=180.*arg(z)/M_PI;    
    if (theta<0.) theta+=360.;
    nrotate=0;
    if (theta>=90. && theta<180.) {
      z*=assign(cos(-M_PI/2.),sin(-M_PI/2.));
      nrotate=3;
    }
    else if (theta>=180. && theta<270.) {
      z*=assign(cos(-2*M_PI/2.),sin(-2*M_PI/2.));
      nrotate=2;
    }    
    else if (theta>=270.) {
      z*=assign(cos(-3*M_PI/2.),sin(-3*M_PI/2.));
      nrotate=1;
    }    

    theta=180.*arg(z)/M_PI;    
    if (theta>90.0000001) {
      printf("ERR theta=%lf, nrotate=%ld\n",theta,nrotate);
      exit(10);
    }

    nguess=5;
    for (iguess=1;iguess<=nguess;iguess++) {
      // initial guess 
      switch(iguess) {
        case 1: w0=assign(0.,1.); break;  
        case 2: w0=assign(-0.5,0.5); break;  
        case 3: w0=assign(0.5,0.5); break;  
        case 4: w0=assign(-1.5,0.5); break;  
        case 5: w0=assign(-1.5,0.5); break;  
      }

      //      if (iguess>1) printf("iguess=%ld \n",iguess);

      w=w0;

      for (i=0;i<20;i++) {
        fw0=scd_square_map(w0);
        if (isnan(abs(fw0)) || abs(fw0)>1.e10) break;
        //printf("i=%ld, error=%le \n",i,abs(fw0-z));
        if (abs(fw0-z)<1.e-11) {
          for (long irotate=0;irotate<nrotate;irotate++) 
            w*=assign(cos(-M_PI/2.),sin(-M_PI/2.));
          return(w); 
        }
        fprime=scd_prefactor*scd_square_integrand(w0);
        dif=(z-fw0)/fprime;
        w=w0+dif;
        // if newton leaves half plane, push it back in 
        if (imag(w)<0.) w=real(w);
        if (isnan(abs(w)) || abs(w)>1.e10) break;
        w0=w;
      }

    }

    printf("FAILURE scd_square_map_inv(), z=%lf+%lfi\n",real(z),imag(z));
    return(acos(2.));
}

complex<double> sc_square_map_inv_test(complex<double> z)
// computes inverse of SC map for square by Newton's method
{
    complex<double> w0,w,fw0,fprime,dif,cent,zz;
    long i,iguess,nguess,nrotate;
    double theta;
    complex<double> prefactor;
 
    // prefactor is related to integral of 1/sqrt(x^4+1) from 0 to 1
    prefactor=assign(0.,1./(0.927037338651*sqrt(2.)));

    // rotate point 
    theta=180.*arg(z)/M_PI;    
    if (theta<0.) theta+=360.;
    nrotate=0;
    if (theta>=90. && theta<180.) {
      z*=assign(cos(-M_PI/2.),sin(-M_PI/2.));
      nrotate=3;
    }
    else if (theta>=180. && theta<270.) {
      z*=assign(cos(-2*M_PI/2.),sin(-2*M_PI/2.));
      nrotate=2;
    }    
    else if (theta>=270.) {
      z*=assign(cos(-3*M_PI/2.),sin(-3*M_PI/2.));
      nrotate=1;
    }    

    theta=180.*arg(z)/M_PI;    
    if (theta>90.0000001) {
      printf("ERR theta=%lf, nrotate=%ld\n",theta,nrotate);
      exit(10);
    }

    nguess=5;
    for (iguess=1;iguess<=nguess;iguess++) {
      // initial guess 
      switch(iguess) {
        case 1: w0=assign(0.,1.); break;  
        case 2: w0=assign(-0.5,0.5); break;  
        case 3: w0=assign(0.5,0.5); break;  
        case 4: w0=assign(-1.5,0.5); break;  
        case 5: w0=assign(-1.5,0.5); break;  
      }

      //      if (iguess>1) printf("iguess=%ld \n",iguess);

      w=w0;

      for (i=0;i<20;i++) {
        fw0=sc_square_map_test(w0);
        if (isnan(abs(fw0)) || abs(fw0)>1.e10) break;
        //printf("i=%ld, error=%le \n",i,abs(fw0-z));
        if (abs(fw0-z)<1.e-11) {
          for (long irotate=0;irotate<nrotate;irotate++) 
            w=sc_square_rotate_90_hplane(w);
          return(w); 
        }
        fprime=prefactor*sc_square_integrand(w0);
        dif=(z-fw0)/fprime;
        w=w0+dif;
        // if newton leaves half plane, push it back in 
        if (imag(w)<0.) w=real(w);
        if (isnan(abs(w)) || abs(w)>1.e10) break;
        w0=w;
      }

    }

    printf("FAILURE sc_square_map_inv_test(), z=%lf+%lfi\n",real(z),imag(z));
    return(acos(2.));
}

double sc_triangle_theta(double frac, double t)
{
  complex<double> z,z1,z2,vtop,vleft,vright,zinv;
  double theta;
  //    z1=(1-frac)*vleft+frac*vtop;  
  //    z2=(1-frac)*vright+frac*vtop;  
  //    z=(1-t)*z1+t*z2;
  //  vtop=assign(0.,2./3.);
  //  vright=assign(1./sqrt(3.),-1./3.);
  //  vleft=assign(-1./sqrt(3.),-1/3.);n

  // inverse SC can fail on the boundary, so the following KLUDGE
  if (t<1.e-5) return(0.);
  if (t>1.-1.e-5) return(1.);

  // OLD:  z=assign((2*t-1)*(1-frac)/sqrt(3.),(3*frac-1)/3.);
  z=assign((1-2*t)*(1-frac)/sqrt(3.),(3*frac-1)/3.);
  zinv=sc_triangle_map_inv(z);
  if (isnan(abs(zinv))) {
    if (t<0.5) return(0.);
    else return(1.);
  }
  theta=arg(zinv)/M_PI;
  if (theta<0) theta+=2.;
  if (theta>1.9999) theta-=2.;
  return(theta);
 }

double sc_triangle_theta_test(double frac, double t)
{
  complex<double> z,z1,z2,vtop,vleft,vright,zinv;
  double theta;
  //    z1=(1-frac)*vleft+frac*vtop;  
  //    z2=(1-frac)*vright+frac*vtop;  
  //    z=(1-t)*z1+t*z2;
  //  vtop=assign(0.,2./3.);
  //  vright=assign(1./sqrt(3.),-1./3.);
  //  vleft=assign(-1./sqrt(3.),-1/3.);n

  // inverse SC can fail on the boundary, so the following KLUDGE
  if (t<1.e-5) return(0.);
  if (t>1.-1.e-5) return(1.);

  // OLD:  z=assign((2*t-1)*(1-frac)/sqrt(3.),(3*frac-1)/3.);
  z=assign((1-2*t)*(1-frac)/sqrt(3.),(3*frac-1)/3.);
  zinv=sc_triangle_map_inv_test(z);
  if (isnan(abs(zinv))) {
    if (t<0.5) return(0.);
    else return(1.);
  }
  theta=arg(zinv)/M_PI;
  if (theta<0) theta+=2.;
  if (theta>1.9999) theta-=2.;
  return(theta);
 }

double sc_square_theta(double frac, double t)
{
  complex<double> z,z1,z2,zinv;
  double theta;
  /*
  if (frac<0.5) {
    z1=assign(-1.,4*frac-1);
    z2=assign(4*frac-1,-1.);
  }
  else {
    z1=assign(4*frac-3.,1.);
    z2=assign(1.,4*frac-3.);
  }
  z=(1-t)*z1+t*z2;
  */

  // inverse SC can fail on the boundary, so the following KLUDGE
  if (t<1.e-5) return(0.);
  if (t>1.-1.e-5) return(1.);

  // OLD:  if (frac<0.5) z=assign(-1 +t*4*frac,(1-t)*4*frac-1);
  // OLD:  else z=assign((1-t)*(4*frac-3.)+t,1-t+t*(4*frac-3.)); 
  if (frac<0.5) z=assign(-1+(1-t)*4*frac,t*4*frac-1);
  else z=assign(t*(4*frac-3.)+(1-t),t+(1-t)*(4*frac-3.)); 

  zinv=sc_square_map_inv(z); 
  if (isnan(abs(zinv))) {
    if (t<0.5) return(0.);
    else return(1.);
  }
  theta=arg(zinv)/M_PI;
  if (theta<0) theta+=2.;
  if (theta>1.9999) theta-=2.;
  return(theta);
 }

double sc_square_theta_test(double frac, double t)
{
  complex<double> z,z1,z2,zinv;
  double theta;
  /*
  if (frac<0.5) {
    z1=assign(-1.,4*frac-1);
    z2=assign(4*frac-1,-1.);
  }
  else {
    z1=assign(4*frac-3.,1.);
    z2=assign(1.,4*frac-3.);
  }
  z=(1-t)*z1+t*z2;
  */

  // inverse SC can fail on the boundary, so the following KLUDGE
  if (t<1.e-5) return(0.);
  if (t>1.-1.e-5) return(1.);

  // OLD:  if (frac<0.5) z=assign(-1 +t*4*frac,(1-t)*4*frac-1);
  // OLD:  else z=assign((1-t)*(4*frac-3.)+t,1-t+t*(4*frac-3.)); 
  if (frac<0.5) z=assign(-1+(1-t)*4*frac,t*4*frac-1);
  else z=assign(t*(4*frac-3.)+(1-t),t+(1-t)*(4*frac-3.)); 
  zinv=sc_square_map_inv_test(z); 
  if (isnan(abs(zinv))) {
    if (t<0.5) return(0.);
    else return(1.);
  }
  theta=arg(zinv)/M_PI;
  if (theta<0) theta+=2.;
  if (theta>1.9999) theta-=2.;
  return(theta);
 }

///////////////////////////////////////////////////////////////
//  checkk checkk checkk checkk checkk checkk checkk checkk  //
///////////////////////////////////////////////////////////////
