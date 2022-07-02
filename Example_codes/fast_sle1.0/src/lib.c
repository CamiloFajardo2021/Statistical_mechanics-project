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

// file name : src/lib.c
//
// simple routines that don't involve any classes

// system includes
#include <cstdlib>
#include "math.h"
#include "sys/time.h"

// local includes
#include "random_number_generator.h"
#include "lib.h"

extern unsigned long long xkiss, ykiss, zkiss, ckiss;

////////////////////////////////////////////////////////////////
//                      simple routines                       //
////////////////////////////////////////////////////////////////

double cauchy_rv()
// Generates a cauchy RV
{
    double x,y;
    x=RNG();
    y=tan(M_PI*(x-0.5));
    return(y);
}

double normal_rv()
// Generates a standard normal RV
// Uses the polar form of the Box-Muller transformation.
{
    double x1, x2, w, y1;

    do {
        x1=2.0*RNG() - 1.0;
        x2=2.0*RNG() - 1.0;
        w=x1*x1 + x2*x2;
    } while (w>=1.0);
    w=sqrt((-2.0*log(w ))/w);
    y1=x1*w;
    return(y1);
}

double mytime() 
// returns time in secs using gettimeofday
// only used for timing, does not affect the computation
{
    double temp;
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp,&tzp);
    temp=tp.tv_sec+tp.tv_usec/1.e6;
    return(temp);
} // end mytime()


void least_squares(long nx,double* y,double* x1,double *sd,
    double* beta,double* beta_err,double *rss) 
// fit to y=beta[1]*x1+beta[0]
// inputs are nx, y, x1, sd 
// output are beta, beta_err, rss
// nx is number of x,y pairs. They should be in x[i],y[i] with i=0,1,...,nx-1
// sd[i] are the standard deviations of the y[i]. 
// NB: the number returned in beta_err is not yet the error in beta
//  It needs to be multiplied by the t-value for alpha/2, nx-2
//
// Confidence intervals: 
//   We use Bonferroni confidence intervals 
//    beta-hat +/- t(alpha/(2k),n-k-1)*sqrt(g_jj)
//    where g_jj is jj entry of    (X^t V^-1 X)^-1 
//    k is number of parameters (=1 or 2 here)
//    n is number of x values 

{
double var_sum=0.,mx=0.,my=0.,mxx=0.,mxy=0.;
double a,b,c,aa,bb,cc;
long i;

for (i=0;i<nx;i++) var_sum+=1./(sd[i]*sd[i]);
for (i=0;i<nx;i++) mx+=x1[i]/(sd[i]*sd[i]);
for (i=0;i<nx;i++) my+=y[i]/(sd[i]*sd[i]);
mxx=0.; for (i=0;i<nx;i++) mxx+=x1[i]*x1[i]/(sd[i]*sd[i]);
mxy=0.; for (i=0;i<nx;i++) mxy+=x1[i]*y[i]/(sd[i]*sd[i]);

a=var_sum;
b=mxx;
c=mx;
// invert 2 by 2 
double det=a*b-c*c;
aa=b/det;
bb=a/det;
cc=-c/det;

beta[0]=aa*my+cc*mxy;
beta[1]=cc*my+bb*mxy;

// still need to multiply these by t value
beta_err[0]=sqrt(aa);
beta_err[1]=sqrt(bb);
 
*rss=0.;
for (i=0;i<nx;i++)  (*rss)+=
  (y[i]-beta[1]*x1[i]-beta[0])*(y[i]-beta[1]*x1[i]-beta[0])/(sd[i]*sd[i]);

} 

///////////////////////////////////////////////////////////////
//  checkk checkk checkk checkk checkk checkk checkk checkk  //
///////////////////////////////////////////////////////////////


