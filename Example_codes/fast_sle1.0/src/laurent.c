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


// file: src/laurent.c

// This class is slightly mis-named. The coefs of the Laurent series
// f(z) about infinity are the coefs of the power series of f(1/z) 
// about 0. This is NOT what the class contains. It contains the 
// power series of f^(z) about 0, where f^(z)=1/f(1/z). 
// A big advantage of working with f^ is that it intertwines with 
// composition in a natural way:
//   (f o g)^ = f^ o g^
// Throughout this file f denotes a function of upper half into itself.

// ?????????????????????????????
// laurent series always have real coefs ???????
// ?????????????????????????????

#include <cstdio> 
using namespace std; 
#include "math.h"
#include <iostream>
using namespace std; 
#include <cstdlib>
#include "sys/time.h"
#include "time.h"
#include "limits.h"
#include <complex>
using namespace std;

#define I complex<double>(0.,1.);


#include "laurent.h"

long laurent::get_nterms()
{
    return nterms;
}

void laurent::fprint(FILE *fptr)
{
    fprintf(fptr,"nterms=%ld\n",nterms);
    for (long i=0;i<=nterms;i++) 
        fprintf(fptr,"coef[%2ld]=%lf+%lf i\n",i,coef[i].real(),coef[i].imag());
}

void laurent::assign(long n,complex<double> *c)
{
    nterms=n;
    for (long i=0;i<=n;i++) coef[i]=c[i];
}

/////////////////////////////////////////////////////////////////////////
// overloaded arithmetic operators:
// NB: these represent addition, subtraction ... of f^, not f
// where f is the map in the subset of the upper half plane 
/////////////////////////////////////////////////////////////////////////

// The overloaded arithemtic operators are all inlined and so are in laurent.h

/////////////////////////////////////////////////////////////////////////
// following are generic in that they do not involve particular series
/////////////////////////////////////////////////////////////////////////

complex<double> laurent::evaluate(complex<double> z)
// evaluates the power series at z; NB: to evaluate f use 
// f(z)=1./fl.evaluate(1./z)), where fl is laurent for f(z)
{
    complex<double> zpow,total;
    total=coef[0];
    zpow=1.;
    for (long k=1;k<=nterms;k++) {
        zpow*=z; // zpow is z^k
        total+=coef[k]*zpow;
    }
    return total;
}


void laurent::compose(laurent f1,laurent f2)
// computes f1(f2(z))
{
    long k;
    laurent f2pow,temp;

    (*this).constant(f1.nterms,f1.coef[0]);
    f2pow.constant(nterms,complex<double>(1.));
    for (k=1;k<=nterms;k++) {
        f2pow*=f2; // f2pow is f2^k
        (*this)+=(f1.coef[k]*f2pow);
    } 
}

void laurent::inverse(laurent f1)
// computes inverse of f1. f1 must start with z. 
// We use f^{-1}(f(z))=z and computer the "remainder" z-f^{-1}(f(z)) 
// order by order
// NB: This routine is not used by the fast_sle package. It is included 
// for possible future utility.
{
    laurent remain,fpow;
    if (abs(f1.coef[0])>1.e-8 || abs(f1.coef[1]-1.)>1.e-8) {
        printf("error in laurent::inverse\n");
        exit(1);
    }
    nterms=f1.nterms;
    remain.identity(nterms);
    coef[0]=0.;
    coef[1]=1.;

    fpow=f1;
    for (long k=2;k<=nterms;k++) {  
        remain-=(coef[k-1]*fpow); 
	// remain now equals z-f^{-1}(f(z)) to order k-1
        coef[k]=remain.coef[k];
        fpow*=f1;
    } 
}

laurent laurent::shift(complex<double> c)
// If (*this) is the power series of f^(z), then shift() returns
// power series of g^(z) where g(z)=f(z)+c. 
// We use g^(z)=h(f^(z)) with h(z)=z/(1+cz). 
{
    laurent f,h;
    h.geometric(-c,nterms);
    f.identity(nterms);
    h*=f;
    f.compose(h,*this);
    return f;
} // end laurent::shift()

/////////////////////////////////////////////////////////////////////////
// following create specific power series 
/////////////////////////////////////////////////////////////////////////

void laurent::constant(long n,complex<double> c)
// creates power series equal to the constant c
{
    nterms=n;
    coef[0]=c;
    for (long i=1;i<=n;i++) coef[i]=0.;
}

void laurent::identity(long n)
// creates power series z
{
    nterms=n;
    coef[0]=0.;
    coef[1]=1.;
    for (long i=2;i<=n;i++) coef[i]=0.;
}

void laurent::geometric(complex<double> c,long n)
// creates geometric series, sum of c^k z^k up to order n
{
    complex<double> cprod;
    nterms=n;
    coef[0]=1.;
    cprod=1.;
    for (long k=1;k<=n;k++) {
        cprod*=c; // cprod is c^k
        coef[k]=cprod;
    }
} // end laurent::geometric()

void laurent::approximate_frac_pow(double alpha, double c,long n)
// power series to order n of the map z -> (1+cz)^alpha. NB no ^
// Warning: this used to compute  z -> (1-cz)^{-alpha}. NB no ^
{
    nterms=n;
    double a=1.;
    coef[0]=1.;
    for (long k=1;k<=nterms;k++) {
        a*=(alpha-k+1)*c/k;
        coef[k]=a;
    }
} // end laurent::approximate_frac_pow()

void laurent::approximate_vertical_slit_ch_zip(double dt,double dx,long n)
// power series to order n of f^(z) where f(z) = sqrt(z^2-4*dt)+dx 
// For dx=0 we have f^(z)=z/sqrt(1-4*dt*z^2)
{
    double a;
    nterms=n;
    coef[0]=0.;
    coef[1]=1.;
    a=1.;
    for (long k=1;2*k+1<=nterms;k++) {
        a*=4*dt*(2.*k-1)/(2.*k);
        coef[2*k+1]=a;
    }
    for (long k=1;2*k<=nterms;k++) coef[2*k]=0.;
    (*this)=(*this).shift(complex<double>(dx));
}

void laurent::approximate_tilted_slit_ch_zip(
    double alpha,double xl,double xr,long n)
// power series to order n of f^(z), f(z) = (z+xl)^(1-alpha)*(z-xr)*alpha 
// f^(z)= z (1+xl*z)^{-(1-alpha)} (1-xr*z)^{-alpha}
{
    laurent f;
    nterms=n;
    (*this).identity(n);
    f.approximate_frac_pow(alpha-1,xl,n);
    (*this)*=f;
    f.approximate_frac_pow(-alpha,-xr,n);
    (*this)*=f;
} 

void laurent::approximate_arc_ch_zip(double a,double b, long n)
{
    // see cap_drive.tex for derivation of following:
    // h^(z)= [ a + (b^2 z-a) [1+2az-b^2 z^2]^{-1/2}] /R^2 
    double caprsq=a*a+b*b;
    laurent f1,f2,f3;

    // f1= -2az+b^2 z^2
    f1.nterms=nterms;
    f1.coef[0]=0.;
    f1.coef[1]= -2.*a;
    f1.coef[2]= b*b;
    for (long i=3;i<=n;i++) f1.coef[i]=0.;

    // f2=[1-z]^{-1/2} 
    f2.approximate_frac_pow(-0.5,-1.,n);

    // (*this)=f2 o f1 = [1+2az - b^2 z^2]^{-1/2}
    (*this).compose(f2,f1);

    // f1=(b^2 z-a) 
    f1.nterms=nterms;
    f1.coef[0]= -a;
    f1.coef[1]= b*b;
    for (long i=2;i<=n;i++) f1.coef[i]=0.;

    (*this) *=f1;
    (*this).coef[0]+=a;
    (*this) *= 1./caprsq;

    // delete following test eventually
    if (abs((*this).coef[0])>1.e-6) {
        printf("error in laurent::approximate_arc_ch_zip()\n");
	printf("%le a=%le b=%le\n",abs((*this).coef[0]),a,b);
        exit(1);
    }
}

void laurent::approximate_vertical_slit_ch_unzip(double dt,double dx,long n)
// power series to order n of f^(z), f(z) = sqrt((z-dx)^2+4 dt)
// f^(z)= ????
// KLUDGE : just uses inverse function - can probably speed this up
// KLUDGE : just uses inverse function - can probably speed this up
// KLUDGE : just uses inverse function - can probably speed this up
// KLUDGE : just uses inverse function - can probably speed this up
{
    (*this).approximate_vertical_slit_ch_zip(dt, dx, n);
    (*this).inverse(*this);
} 

void laurent::approximate_tilted_slit_ch_unzip(
    double alpha,double xl,double xr,long n)
// power series to order n of f^(z), f(z) = (z+xl)^(1-alpha)*(z-xr)*alpha 
// f^(z)= z (1+xl*z)^{-(1-alpha)} (1-xr*z)^{-alpha}
{
    (*this).approximate_tilted_slit_ch_zip(alpha, xl, xr, n);
    (*this).inverse(*this);
} 

void laurent::approximate_arc_ch_unzip(double a,double b,long n)
{
    // see cap_drive.tex for derivation of following:
    // g^(z)=a/b^2 + (R^2 z - a) [1-2 a z + R^2 z^2]^{-1/2} / b^2 
    double caprsq=a*a+b*b;
    laurent f1,f2,f3;
    complex<double> c;
    long i;

    // If b is very small we will take a different approach
    if (b>1.e-4*fabs(a)) { // not sure this cutoff is optimal ?????
        // f1= 2 a z - R^2 z^2
        f1.nterms=nterms;
        f1.coef[0]=0.;
        f1.coef[1]= 2.*a;
        f1.coef[2]= -caprsq;
        for (i=3;i<=n;i++) f1.coef[i]=0.;

        // f2=[1-z]^{-1/2} 
        f2.approximate_frac_pow(-0.5,-1.,n);

        // (*this)=f2 o f1 = [1-2 a z + R^2 z^2]^{-1/2}
        (*this).compose(f2,f1);

        // f1=(R^2 z - a) 
        f1.nterms=nterms;
        f1.coef[0]= -a;
        f1.coef[1]= caprsq;
        for (i=2;i<=n;i++) f1.coef[i]=0.;

        (*this) *=f1;

        (*this) *= 1./(b*b);
        (*this).coef[0]+=a/(b*b);
        // delete following test eventually
        if (abs((*this).coef[0])>1.e-6) {
            printf("error in laurent::approximate_arc_ch_unzip()\n");
            printf("coef=%le a=%le b=%le\n",abs((*this).coef[0]),a,b);
            // exit(1);
        }
    }
    else {  
        // for small b we do an expansion in b:
        // g^(z)=z/(1-az) + (a/2) z^2/(1-az)^2-(b^2/2) z^3/(1-az)^3
        //       -(3ab^2/8) z^4/(1-az)^4
        c=0.;
	(*this).constant(n,c); 

        // z/(1-az) 
        f1.nterms=nterms;
        for (i=0;i<=n;i++) f1.coef[i]=0.;
        f1.coef[1]=1.;
        f2.approximate_frac_pow(-1.,-a,n);
        (*this)+=f1*f2;

        //  (a/2) z^2/(1-az)^2
        f1.nterms=nterms;
        for (i=0;i<=n;i++) f1.coef[i]=0.;
        f1.coef[2]=a/2.;
        f2.approximate_frac_pow(-2.,-a,n);
        (*this)+=f1*f2;

        //  -(b^2/2) z^3/(1-az)^3
        f1.nterms=nterms;
        for (i=0;i<=n;i++) f1.coef[i]=0.;
        f1.coef[3]=-b*b/2.;
        f2.approximate_frac_pow(-3.,-a,n);
        (*this)+=f1*f2;

        // -(3ab^2/8) z^4/(1-az)^4
        f1.nterms=nterms;
        for (i=0;i<=n;i++) f1.coef[i]=0.;
        f1.coef[4]=-3.*a*b*b/8.;
        f2.approximate_frac_pow(-4.,-a,n);
        (*this)+=f1*f2;

    }
}


void laurent::approximate_vertical_slit_rad_zip(double dt,double dx,long n)
// map_choice=1
// power series of 
// Is this used?????????????????????????????????
// Is this used?????????????????????????????????
// Is this used?????????????????????????????????
// Is this used?????????????????????????????????
// Is this used?????????????????????????????????
// Is this used?????????????????????????????????
// I don't think so 
{
    laurent f1,f2,f3;
    complex<double> a,c,d;

    c= -dx*I;
    c=exp(c); 
    d=exp(-dt);
    // z -> exp(i*dx) z, then map to right half plane: z -> (1-z)/(1+z) 
    // f1 = (1-c*z)/(1+c*z) = 1-2*c*z/(1+c*z)
    f1.nterms=n;
    f1.coef[0]=1.;
    a=2.;
    for (long k=1;k<=f1.nterms;k++) {
        a*= -c;
        f1.coef[k]=a;
    }
    // insert hor slit from 1 to left: z -> sqrt(exp(-dt)*zz*zz+1.-exp(-dt));
    // f2=exp(-dt)*f1*f1-exp(-dt), so above is sqrt(1+f2)
    f2=d*f1*f1;
    f2.coef[0]-=d;
    // f3=sqrt(1+f2)
    f1.approximate_frac_pow(0.5,1.,n);
    f3.compose(f1,f2);
    // map to disc:  z -> (1-z)/(1+z). So we have (1-f3)/(1+f3) 
    // f3=1 at zero order.Let f3=1-f3; we now have f3/(2-f3)=0.5*f3/(1-0.5*f3) 
    f3= -1.*f3;
    f3.coef[0]+=1.;    

    // f2=1/(1-0.5*f3)
    f1.geometric(0.5,n);
    f2.compose(f1,f3);
    (*this)=f3*f2;
    (*this)*=0.5;

    //    f2.fprint(stdout); exit(0);
}

// see saw/versionN/tex/radial.tex for documentation
void laurent::approximate_vertical_slit_half_rad_zip(double dt,double dx,long n)
// map_choice=1
// hat power series of S(f(z)) where 
// f(z) = (z+alpha)/(1+alpha*z) 
//  alpha=(1-exp(-i dx))/(1+exp(-i dx))
//  S(z)=sqrt(exp(-dt)*z*z+1-exp(-dt))
//      =exp(-dt/2)*sqrt(z*z+exp(dt)-1)
{
    laurent f1,f2,f3;
    complex<double> c,alpha,d;

    c= I;
    c=exp(-c*dx);    
    alpha=(1.-c)/(1.+c);

    // f1 is f above 
    f1.nterms=n;
    f1.coef[0]=alpha;
    d=1.;
    for (long k=1;k<=f1.nterms;k++) {
        f1.coef[k]=d*alpha*alpha-d;
        if (k%2==1) f1.coef[k]*=-1.;
        d*=alpha;
    }

    f2.approximate_vertical_slit_ch_zip(0.25*(1.-exp(dt)),0.,n);
    f2*=exp(dt/2.);

    (*this).compose(f2,f1);

    // we now have map in right half plane. To get it for upper half plane 
    // we need to multiply nth coef by -i^{n+1}
    c=-I;
    for (long k=0;k<=nterms;k++) {
        coef[k]*=c;
        c*=I;
    }
}


///////////////////////////////////////////////////////////////
//  checkk checkk checkk checkk checkk checkk checkk checkk  //
///////////////////////////////////////////////////////////////



