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


// driver_sle.c
//
// Minimal program to run sle simulation. 
// User must input all the parameters, but program offers sample values.
//
// We use nsteps*dsteps iterations of the conformal map, but 
// we only compute the point on the SLE every dsteps maps. So the 
// resulting walk has nsteps. It is written to the file "sle.plt"

// system includes
#include <cstdio> 
using namespace std; 
#include "math.h"
#include <iostream>
using namespace std; 
#include <cstdlib>
#include "limits.h"
#include <complex>
using namespace std;
#define I complex<double>(0.,1.);

// local includes
#include "random_number_generator.h"
#include "complex_walk.h"
#include "conformal_maps.h"

unsigned long long xkiss,ykiss,zkiss,ckiss;

int main()
{
long ns;
complex_walk cw;           // the SLE path 
double kappa;              // the SLE parameter you know and love

// parameters that control how SLE is discretized
long nsteps;               // number of points computed on the SLE path
long max_nsteps;           // used for memory allocation
long dsteps;               // number of time intervals between points
long map_choice;           // conformal map - vertical vs tilted slits
long dt_choice;            // time discretization - uniform vs. non-uniform
long dx_choice;            // random walk steps - bernouli vs normal

// parameters that control how laurent series are implemented
long laurent_flag=0;       // exact vs laurent
long blength;              // number of functions in a block
long nterms;               // order of laurent series 
double laurent_factor=1.6; // controls when laurent series are used

long iseed; //random number generator seed, -1 makes system generate it

// irrelevant parameters - these play no role in this program
double nu;

complex<double> z1,z2;

long chordal=1;

printf("In the following a sample value is given in (...) \n");

printf("kappa (2.6667) \n");
ns=scanf("%lf",&kappa);
nu=1./(1.+kappa/8.);

printf("number of points on SLE trace you wish to compute (2000) \n");
ns=scanf("%ld",&nsteps);
max_nsteps=nsteps;

printf("number of time intervals in between points (5)\n");
ns=scanf("%ld",&dsteps);

printf("choice of conformal map: 1=vertical slit, 2=tilted slit  (2) \n");
ns=scanf("%ld",&map_choice);

printf("type of discretization of time: 1=uniform, 2=non-uniform  (2)\n");
ns=scanf("%ld",&dt_choice);

printf("type of steps for the random walk: 1=bernoulli, 2=normal (2) \n");
ns=scanf("%ld",&dx_choice);

printf("0 for an exact computation, 1 to use laurent series (1) ");
ns=scanf("%ld",&laurent_flag);
if (laurent_flag==1) {
    printf("Following control just how laurent series are used \n");
    printf("enter number of functions in a block (20) \n");
    ns=scanf("%ld",&blength);

    printf("enter number of terms in laurent series (10) \n");
    ns=scanf("%ld",&nterms);

    printf("enter laurent_factor  (3.) \n");
    ns=scanf("%lf",&laurent_factor);
}

printf("enter seed for random number generator, (-1)\n");
ns=scanf("%ld",&iseed);

initialize_random_number_generator(iseed);

// print the parameters 
printf(
    "kappa=%lf, will compute %ld points, %ld time intervals in between\n",
     kappa,nsteps,dsteps);
switch (map_choice) {
    case 1: printf("Conformal map produces vertical slit, "); break;
    case 2: printf("Conformal map produces tilted slit, "); break;
}
switch (dt_choice) {
    case 1: printf("uniform time discretization, "); break;
    case 2: printf("non-uniform time discretization, "); break;
}
switch (dx_choice) {
    case 1: printf("Bernoulli steps\n"); break;
    case 2: printf("normally distributed steps\n"); break;
}
if (laurent_flag) printf("Using laurent: blength=%ld, nterms=%ld, L=%lf\n",
    blength,nterms,laurent_factor);
else  printf("Not using laurent series\n");

printf("Beginning computation \n");

cw.allocate(nsteps,dsteps,map_choice,blength,chordal);   // allocate memory  

long bm_approx=0;
// equal dt, bernoulli steps
if (dt_choice==1 && dx_choice==1) bm_approx=1;

// dt as d(t^1.5), bernoulli steps
if (dt_choice==2 && dx_choice==1) bm_approx=2;

// equal dt, normal distribution for steps 
if (dt_choice==1 && dx_choice==2) bm_approx=3;

// dt as d(t^1.5), normal distribution for steps 
if (dt_choice==2 && dx_choice==2) bm_approx=4;

double min_dt=0.;
double max_dt=0.;
double max_dx=0.;
long sle_domain=1;

if (chordal) {
    printf("map_choice=%ld\n",map_choice);
    printf("bm_approx=%ld\n",bm_approx);
    cw.adapt_chordal_sle_laurent(sle_domain,z1,z2,nsteps,dsteps,
        map_choice, bm_approx,kappa,0.,0L,0.,nu,min_dt,max_dt,max_dx,
        blength,nterms,laurent_factor);
} 
else { // radial sle 
    printf("radial sle not implemented in adapt_sle \n");    
    exit(1);
}

// write the SLE trace to the file "sle.plt" 
FILE *fptr;
fptr=fopen("sle.plt","w");
cw.plot(fptr);
fclose(fptr);

printf("Finished computation\n");
printf("File sle.plt contains the sle trace \n");

ns++;
} // end main

///////////////////////////////////////////////////////////////
//  checkk checkk checkk checkk checkk checkk checkk checkk  //
///////////////////////////////////////////////////////////////
