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


// file : src/complex_walk.c 

// the class complex_walk represents a walk in the complex plane

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
#include <complex>
using namespace std;
 
#define I complex<double>(0.,1.);

// local includes, only include what is needed, 
#include "local_constants.h"
#include "random_number_generator.h"
#include "lib.h"
#include "laurent.h"
#include "conformal_maps.h"
#include "complex_walk.h"

//////////////////////////////////////////////////////////////////
//        non-member functions for class complex_walk           //
//////////////////////////////////////////////////////////////////

double max_distance(complex_walk *cw1,complex_walk *cw2) 
// computes maximum separation of points on two walks
{
    double temp,sup=0.;
    for (long i=0;i<=(*cw1).nsteps;i++) {
        temp=abs((*cw1).steps[i]-(*cw2).steps[i]);
        if (temp>sup) sup=temp;
    }
    return(sup);
}

double average_distance(complex_walk *cw1,complex_walk *cw2) 
// computes average separation of points on two walks; not l1
{
    double sum=0.;
    for (long i=1;i<=(*cw1).nsteps;i++) {
        sum+=abs((*cw1).steps[i]-(*cw2).steps[i]);
    }
    sum/=double((*cw1).nsteps);
    return(sum);
}

////////////////////////////////////////////////////////////////
//             member functions for class complex_walk        //
////////////////////////////////////////////////////////////////

complex_walk::complex_walk() 
// NB this constructor does no memory allocation. Must call allocate(...)
{
} 

void complex_walk::allocate(long n, long ds, long mc, long bl, long st) 
{
    num_unzip=0;
    sle_type=st;

    steps = new complex<double>[n+1];
    if (steps==0) {
        cout << "allocate error in complex_walk::allocate \n";
        exit(1);
    }

    dx = new double[(n+1)*ds];
    dt = new double[(n+1)*ds];
    t  = new double[n+1];

    sewing_minus = new double[(n+1)*ds];
    sewing_plus = new double[(n+1)*ds];

    map_choice=mc;
    max_nsteps=n; 
    dsteps=ds;
    nsteps=0;
    bm_nsteps=0;
    nforce_pt=-1; // to flag that nforce_pt is not set

    blength=bl;
    if (blength==0) nblock=0; // blength==0 means don't use laurent
    else nblock=(max_nsteps*dsteps)/blength;

    switch(sle_type) {
    case 1: // chordal
        switch(map_choice) {
        case 1:
            fblock= new laurent[nblock+1];
            radius= new double[nblock+2];
        break;
        case 2:
            alpha = new double[max_nsteps*dsteps+1];
            xl = new double[max_nsteps*dsteps+1];
            xr = new double[max_nsteps*dsteps+1];
            fblock= new laurent[nblock+1];
            radius= new double[nblock+2];
        break;
        case 3:
            aa = new double[max_nsteps*dsteps+1];
            bb = new double[max_nsteps*dsteps+1];
            fblock= new laurent[nblock+1];
            radius= new double[nblock+2];
        break;
        }
    break;
    case 2: 
        switch(map_choice) {
        case 1:
            fblock= new laurent[nblock+1]; // delete this ??????
        break;
        case 2:
            alpha = new double[max_nsteps*dsteps+1];
            xl = new double[max_nsteps*dsteps+1];
            xr = new double[max_nsteps*dsteps+1];
            ww = new complex<double>[max_nsteps*dsteps+1];
            fblock= new laurent[nblock+1]; // delete this ??????
        break;
        case 3:
        break;
        }
    break;
    case -2:
        switch(map_choice) {
        case 1:
            fblock= new laurent[nblock+1];
            radius= new double[nblock+2];
        break;
        case 2:
            alpha = new double[max_nsteps*dsteps+1];
            xl = new double[max_nsteps*dsteps+1];
            xr = new double[max_nsteps*dsteps+1];
            ww = new complex<double>[max_nsteps*dsteps+1];
            fblock= new laurent[nblock+1];
            radius= new double[nblock+2];
        break;
        case 3:
            printf("map choice not implemented asefdi\n");
            exit(10);
        break;
        }
    break;
    case 3: 
        switch(map_choice) {
        case 1:
            fblock= new laurent[nblock+1];
        break;
        }
    break;
    }
} 

void complex_walk::deallocate()
{
    delete [] steps;
    delete [] dx;
    delete [] dt;
    delete [] t;
    delete [] sewing_minus;
    delete [] sewing_plus;

    switch(sle_type) {
    case 1:
        switch(map_choice) {
        case 1:
            delete [] radius;
            delete [] fblock;
        break;
        case 2:
            delete [] radius;
            delete [] fblock;
            delete [] alpha;
            delete [] xl;
            delete [] xr;
        break;
        case 3:
            delete [] radius;
            delete [] fblock;
            delete [] aa;
            delete [] bb;
        break;
        }
    break;
    case 2:
    case -2:
        switch(map_choice) {
        case 1:
            delete [] fblock;
        break;
        case 2:
            delete [] fblock;
            delete [] alpha;
            delete [] xl;
            delete [] xr;
            delete [] ww;
        break;
        case 3:
        break;
	}
    break;
    case 3:
        switch(map_choice) {
        case 1:
            delete [] fblock;
        break;
	}
    break;
    }
}; 

void complex_walk::print(FILE *fptr) 
// prints the number of steps and the complex values of the complex_walk
{
    fprintf(fptr,"nsteps=%ld\n",nsteps);
    for (long i=0;i<=nsteps;i++) {
        fprintf(fptr,"%lf + %lf i\n",steps[i].real(),steps[i].imag());
    }
} 

void complex_walk::plot(FILE *fptr)
// prints only the coordinates of the complex_walk for plotting
{ 
    for (long i=0;i<=nsteps;i++) {  
        fprintf(fptr,"%14.10le %14.10le\n",steps[i].real(),steps[i].imag());
    }
} 

void complex_walk::plot_skip(FILE *fptr,long dn)
// prints only the coordinates of the complex_walk for plotting 
// only plots every dn
{ 
    for (long i=0;i<=nsteps;i+=dn) {  
        fprintf(fptr,"%14.10le %14.10le\n",steps[i].real(),steps[i].imag());
    }
} 

long complex_walk::get_nsteps() 
{
    return(nsteps);
}

void complex_walk::set_nsteps(long n) 
{
    nsteps=n;
}

void complex_walk::set_dsteps(long n) 
{
    dsteps=n;
}

void complex_walk::set_max_nsteps(long n) 
{
    max_nsteps=n;
}

complex<double> complex_walk::get_step(long i) 
{
    return(steps[i]);
}

void complex_walk::assign_step(long i,complex<double> z)
{
    steps[i]=z;
}

void complex_walk::invert() 
// applies the map z -> -1/z to the walk and reverses it order
// assumes the walk starts at 0, and we just leave this 0 there
{
    // reverse the order of the walk
    complex<double> temp;
    for (long i=1;i<=nsteps/2;i++) {
        temp=steps[i];
	steps[i]=steps[nsteps-i+1];
        steps[nsteps-i+1]=temp;
    }
    // invert it
    for (long i=1;i<=nsteps;i++) steps[i]=-1./steps[i];
} 

/////////////////////////////////////////////////////////////////////////////
// Functions to generate approximation of root(kappa)*brownian motion.
// num_dt is the number of time intervals 
// dt[k], k=1,2,...,num_dt are the intervals.
// dx[k] is change in driving function over that time interval.
// NB : num_dt is not the number of points computed for the SLE path
// Instead, num_dt=nsteps*dsteps where nsteps is the number of SLE points
// and dsteps is the number of time intervals in between points.
/////////////////////////////////////////////////////////////////////////////

void complex_walk::read_bm(FILE *fptr)
{
    long ns;
    ns=fscanf(fptr,"%ld\n",&dsteps);
    ns=fscanf(fptr,"%ld\n",&bm_nsteps);

    for (long i=1;i<=bm_nsteps*dsteps;i++) {
        if (fscanf(fptr,"%le %le",dt+i,dx+i)!=2) {
                printf("error in read_bm\n");
                exit(10);
            }
        else {
            dx[i]*=sqrt(kappa);
	}
    }
    ns++;
}


void complex_walk::write_bm(FILE *fptr)
{
    fprintf(fptr,"%ld\n",dsteps);
    fprintf(fptr,"%ld\n",bm_nsteps);
    for (long i=1;i<=bm_nsteps*dsteps;i++) 
        fprintf(fptr,"%le %le\n",dt[i],dx[i]/sqrt(kappa));
}

void complex_walk::cauchy_deterministic_dt()
{
    long istep,num_dt;
    double t,incr;
    num_dt=dsteps*nsteps;

    for (istep=1;istep<=num_dt;istep++) {
        dt[istep]=1./double(num_dt);
    }
    t=0.;
    for (istep=1;istep<=num_dt;istep++) {
        t+=dt[istep];
        incr=cauchy_rv();
        dx[istep]=0.5*incr*kappa*dt[istep]/sqrt(t);
    } // end loop on istep
} 

void complex_walk::bm_drift_deterministic_dt()
{
    long istep,num_dt;
    double t;
    num_dt=dsteps*nsteps;
    bm_nsteps=num_dt;

    for (istep=1;istep<=num_dt;istep++) {
        dt[istep]=1./double(num_dt);
    }
    t=0.;
    for (istep=1;istep<=num_dt;istep++) {
        t+=dt[istep];
        dx[istep]=normal_rv()*sqrt(kappa*dt[istep])+drift*dt[istep];
    } // end loop on istep
} 

void complex_walk::bm_deterministic_dt(long dt_choice,long dx_choice)
// set up array of dt for sle generation
// dt_choice==1: time intervals are all the same
// dt_choice==2: times are at (k/num_dt)^{3/2}, k=1,2,...,num_dt
//     This makes points on SLE roughly equally spaced.
// set up array of dx for sle generation
// dx_choice==1: dx is just +/- sqrt(kappa*dt) (Bernoulli)
// dx_choice==2: dx is normal with variance kappa*dt
{
    double temp;
    long istep,num_dt;
    double incr;
    num_dt=dsteps*nsteps;
    bm_nsteps=num_dt;

    switch (dt_choice) {
        case 1: 
            for (istep=1;istep<=num_dt;istep++) {
                dt[istep]=1./double(num_dt);
	    }
        break;
        case 2: 
            temp=0.;
            for (istep=1;istep<=num_dt;istep++) {
                dt[istep]=-temp;
                temp=sqrt(double(istep)/double(num_dt))
                    *double(istep)/double(num_dt);
                dt[istep]+=temp;
            }
        break;
        default: 
            printf("bad dt_choice - bm_deterministic_dt() %ld \n",dt_choice);
            exit(1);
        break;
    }

    for (istep=1;istep<=num_dt;istep++) {
        switch(dx_choice) {
            case 1: // Bernoulli
                if (RNG()>0.5) incr=1.;
                else incr=-1.;
            break;
            case 2: // normal 
                incr=normal_rv();
            break;
            default: 
               printf("bad dx_choice \n"); 
               exit(0);
            break;
        }
        dx[istep]=incr*sqrt(kappa*dt[istep]);
    } // end loop on istep
} 

/////////////////////////////////////////////////////////////////////////////
// member functions to generate chordal SLE by iterating conformal maps
//
// NB: we use nsteps*dsteps iterations of the conformal map, but 
//  we only compute the point on the SLE every dsteps maps. So the 
//  resulting walk has nsteps points. 
//  For chordal SLE we usually use the time interval [0,1] since SLE for 
//  the time interval [0,T] is trivially related by scaling.  
//  For some rv's we can stop compute before t reaches 1.
//
// *var is used to return the fractal (1/nu) variation.
//
// Parameters: following routines use following parms:
//  long nsteps             // number of points computed on SLE trace
//  long dsteps             // number of time intervals between above point
//  double kappa            // the usual SLE parameter
//  double *var             // returns 1/nu variation of walk
//  double nu               // only used in computation of 1/nu variation
//  long bm_approx          // how bm is discretized 
// Following only used for laurent series routines
//  long blength            // length of blocks of functions in composition
//  long nterms             // order of Laurent series
//  double laurent_factor   // controls when Laurent series are used

void complex_walk::sample_bm(double max_dt, double *next_dt)
{
    // arrays dx[] and dt[] contain the increments of the driving function,
    // i.e., root(kappa)*Browian motion. 
    // Their indices start at 1; the final index is always multiple of dsteps.
    // When routine is called, the BM contains dsteps*bm_nsteps increments.
    // When routine is done, the number of increments is either the same or 
    // bm_nsteps is increased by 1. 
    // nsteps is the index of the point we are trying to add to the SLE.
    // The smallest bm_nsteps should be is nsteps-1, indicating the 
    // BM has not been sampled at all in the future and we just add on 
    // gaussian steps. We take the dt to be max_dt.
    // If bm_nsteps is at least nsteps, then the BM has already been 
    // sampled in the future. What we do depends on next_dt. If it is 
    // no smaller than the first dt in the future, we just use the 
    // next dsteps increments in the BM. We also reduce next_dt to 
    // the first dt in the future.
    // If it is smaller than the first dt in the future, we need
    // to do a more refined sampling of the BM in a time interval where it 
    // has already been sampled and so must use a Brownian bridge. 

    long istep;

    if (bm_nsteps==nsteps-1) {
        // bm has not been sampled beyond t; sample it with max_dt
        for (istep=(nsteps-1)*dsteps+1;istep<=nsteps*dsteps;istep++) {
            dt[istep]=max_dt;
            dx[istep]=normal_rv()*sqrt(kappa*dt[istep]);
        }
        bm_nsteps++; 
        return;
    }

    if (bm_nsteps>=nsteps && (*next_dt)/dt[(nsteps-1)*dsteps+1]>1.-1.e-10) {
        // bm has been sampled beyond t, next_dt >= next dt in bm 
        // no sampling of BM - use what is there
        *next_dt=dt[(nsteps-1)*dsteps+1];
        return;
    }

    if (bm_nsteps>=nsteps 
	&& fabs((*next_dt)/dt[(nsteps-1)*dsteps+1]-0.5)<1.e-10) {
        // bm has been sampled beyond t, next_dt < next dt in bm 
        // sample BM using Brownian bridge
        for (istep=bm_nsteps*dsteps;istep>=nsteps*dsteps+1;istep--) {
            dt[istep+dsteps]=dt[istep];
            dx[istep+dsteps]=dx[istep];
        }

        long n=(nsteps-1)*dsteps;
        for (long j=dsteps;j>=1;j--) {
            dt[n+2*j-1]=*next_dt;
            dt[n+2*j]=*next_dt;
            dx[n+2*j]=dx[n+j];
            //         <----  dt[n+2*j-1] -----> <------- dt[n+2*j] ------>
            //         |                        |                         |
            // before  <------------------ dx[n+2*j] --------------------->
            // after   <------  dx[n+2*j-1] ------> <----- dx[n+2*j] ----->
            double mean=dx[n+2*j]/2.;
            double variance=kappa*dt[n+2*j-1]/2.;
            dx[n+2*j-1]=mean+normal_rv()*sqrt(variance);
	    dx[n+2*j]-=dx[n+2*j-1];
        }

        bm_nsteps++; 
        return;
    }

    printf("error in sample_bm,nsteps=%ld bm_nsteps=%ld \n",nsteps,bm_nsteps);
    printf("next_dt=%le, dt=%le\n",(*next_dt),dt[(nsteps-1)*dsteps+1]);
    exit(10);
}

/////////////////////////////////////////////////////////////////////////////
// General types of computation: 
//    CHORDAL ZIP = ch_zip
//    CHORDAL UNZIP = ch_unzip
//    RADIAL ZIP = rad_zip
//    DIPOLAR ZIP = di_zip
//
// CHORDAL ZIP: conformal map takes half plane to half plane minus a curve. 
//  The important example of this is generating chordal SLE:
//    adapt_chordal_sle_laurent()
//  The above creates the driving function as well as zipping.
//  If we already have some driving function, the zipping is done by
//    ch_zip()
//  It calls 
//    compute_parms_ch_zip
//    conformal_map_ch_zip
//    compute_block_ch_zip
//  which in turn call 
//    (approximate_)vertical_slit_ch_zip 
//    (approximate_)tilted_slit_ch_zip 
//    (approximate_)arc_ch_zip
//
//
// CHORDAL UNZIP: conformal map takes half plane minus a curve to half plane.
//  The important example of this is computing the capacity or driving 
//  function of a given curve:
//    ch_unzip() 
//  It calls
//    compute_parms_ch_unzip 
//    conformal_map_ch_unzip 
//    compute_block_ch_unzip 
//  which in turn call 
//    (approximate_)vertical_slit_ch_unzip 
//    (approximate_)tilted_slit_ch_unzip
//    (approximate_)arc_ch_unzip 
//
// RADIAL ZIP: the conformal map takes disc to disc minus a curve. 
//  The important example of this is generating radial SLE: 
//    adapt_radial_sle_laurent() 
//  It calls 
//    compute_parms_rad_zip   
//    conformal_map_rad_zip  
//    compute_block_rad_zip  <<<<<<<<<<<<<<<<<<<<<<<< not done
//  which in turn call 
//    (approximate_)vertical_slit_rad_zip 
//    (approximate_)tilted_slit_rad_zip 
//    (approximate_)arc_rad_zip 
//
// HALF_RADIAL ZIP: radial case in the half plane
//  So the conformal map takes half plane to half plane minus a curve,
//  fixing i rather than infinity.  
//  Note that for vertical slit we use right half plane rather than upper 
//  half plane. 
//  The important example of this is generating radial SLE: 
//    adapt_half_radial_sle_laurent() 
//  It calls 
//    compute_parms_half_rad_zip   
//    conformal_map_half_rad_zip  
//    compute_block_half_rad_zip  
//  which in turn call 
//    (approximate_)vertical_slit_half_rad_zip 
//    (approximate_)tilted_slit_half_rad_zip 
//
// DIPOLAR ZIP : the conformal map takes strip to strip minus a curve
//    adapt_dipolar_sle_laurent() 
//  It calls 
//    compute_parms_di_zip  DONE
//    conformal_map_di_zip  DONE
//    compute_block_di_zip
//  which in turn call 
//    (approximate_)vertical_slit_di_zip
//    (approximate_)tilted_slit_di_zip
//    (approximate_)arc_di_zip
//
/////////////////////////////////////////////////////////////////////////////

void complex_walk::compute_block_ch_unzip(long ib)
{
    long i;
    laurent ftemp;

    fblock[ib].identity(nterms);
    for (i=(ib-1)*blength+1;i<=ib*blength;i++) {
        switch(map_choice) {
            case 1: // vertical slit
                ftemp.approximate_vertical_slit_ch_unzip(dt[i],dx[i],nterms);
            break;
            case 2: // tilted slit
                 ftemp.approximate_tilted_slit_ch_unzip(
                     alpha[i],xl[i],xr[i],nterms);
            break;
            case 3: // arc
               if (fabs(aa[i])+fabs(bb[i])>1.e-8) 
                   ftemp.approximate_arc_ch_unzip(aa[i],bb[i],nterms);
               else // KLUDGE to avoid calling with aa=bb=0
                   ftemp.identity(nterms);
            break;
        }
        fblock[ib].compose(ftemp,fblock[ib]);
    }
}

void complex_walk::compute_block_ch_zip(long ib)
{
    long i,real_flag;
    laurent ftemp;
    complex<double> exact;
    double xreal=0.,xcomplex,xmid=0.;

    fblock[ib].identity(nterms);
    for (i=ib*blength;i>=(ib-1)*blength+1;i--) {
        switch(map_choice) {
            case 1:
               ftemp.approximate_vertical_slit_ch_zip(dt[i],dx[i],nterms);
            break;
            case 2:
               ftemp.approximate_tilted_slit_ch_zip(alpha[i],xl[i],xr[i],
                   nterms);
            break;
            case 3:
               ftemp.approximate_arc_ch_zip(aa[i],bb[i],nterms);
            break;
        }
        fblock[ib].compose(ftemp,fblock[ib]);
    }
   
    if (map_choice>1) {// find radius: smallest x s.t. x, -x have real images
        xcomplex=0.;
        xreal=1.;
        real_flag=0; 
        while (xreal-xcomplex>1.e-5*xreal) {
            if (!real_flag) xreal*=2.;
            xmid=(xreal+xcomplex)/2.;
            long complex_flag=0;
            for (long sgn=-1;sgn<=1 && complex_flag==0;sgn+=2) {
                exact=sgn*xmid;
                for (i=ib*blength;i>=(ib-1)*blength+1;i--) {
                    switch(map_choice) {
                        case 2:
                            exact=tilted_slit_ch_zip(alpha[i],xl[i],xr[i],
                                exact);
                        break;
                        case 3:
                            exact=arc_ch_zip(aa[i],bb[i],exact);
                        break;
                    }
                    if (exact.imag()>1.e-8*xreal) break;
                }
                if (exact.imag()>1.e-8*xreal) complex_flag=1;
            } // end loop on sgn
            if (!complex_flag) real_flag=1;

            if (complex_flag) xcomplex=xmid;
            else xreal=xmid;
        }
        radius[ib]=fabs(xmid); 
    }
    else {    //  brute force comp of radius 
        double maxx=1.e-5;
        double xfactor=1.05;
        for (long sgn=-1;sgn<=1;sgn+=2) 
        for (double xx=maxx;xx<=10.;xx*=xfactor) {
            exact=sgn*xx;
            for (i=ib*blength;i>=(ib-1)*blength+1;i--) {
                exact=vertical_slit_ch_zip(dt[i],dx[i],exact);
// ???????????????????????????????????????????????????
// xreal is 0. how did following ever work? replace xreal by xx ???
                if (exact.imag()>1.e-8*xreal) break;
            }
// ???????????????????????????????????????????????????
// xreal is 0. how did following ever work? replace xreal by xx ???
            if (exact.imag()>1.e-8*xreal && fabs(xx)>maxx) maxx=fabs(xx); 
        } 
        maxx*=xfactor;
        radius[ib]=maxx;
        printf("radius=%lf \n",radius[ib]);
    }
}

void complex_walk::compute_block_half_rad_zip(long ib)
{
    long i,real_flag;
    laurent ftemp;
    complex<double> exact;
    double xreal=0.,xcomplex,xmid=0.;

    fblock[ib].identity(nterms);
    for (i=(ib-1)*blength+1;i<=ib*blength;i++) {
        switch(map_choice) {
            case 1: // vertical slit
                ftemp.approximate_vertical_slit_half_rad_zip(
                    dt[i],dx[i],nterms);
            break;
/*
NOT IMPLEMENTED YET
            case 2: // tilted slit
                ftemp.approximate_tilted_slit_half_rad_zip(
                    dx[i],alpha[i],xl[i],xr[i],ww[i],nterms);
            break;
*/
            default: 
                printf("bad map_choice in compute_block_half_rad_zip()\n");
                exit(10);
            break;
        }
        fblock[ib].compose(ftemp,fblock[ib]);
    }

    if (map_choice>1) {// find radius: smallest x s.t. x, -x have real images
        xcomplex=0.;
        xreal=1.;
        real_flag=0; 
        while (xreal-xcomplex>1.e-5*xreal) {
            if (!real_flag) xreal*=2.;
            xmid=(xreal+xcomplex)/2.;
            long complex_flag=0;
            for (long sgn=-1;sgn<=1 && complex_flag==0;sgn+=2) {
                exact=sgn*xmid;
                for (i=ib*blength;i>=(ib-1)*blength+1;i--) {
                    switch(map_choice) {
                        case 2:
                            exact=tilted_slit_half_rad_zip(dx[i],alpha[i],
                                xl[i],xr[i],ww[i],exact);
                        break;
                    }
                    if (exact.imag()>1.e-8*xreal) break;
                }
                if (exact.imag()>1.e-8*xreal) complex_flag=1;
            } // end loop on sgn
            if (!complex_flag) real_flag=1;

            if (complex_flag) xcomplex=xmid;
            else xreal=xmid;
        }
        radius[ib]=fabs(xmid); 
    }
    else {    //  brute force comp of radius 
        double maxx=1.e-5;
        double xfactor=1.05;
        for (long sgn=-1;sgn<=1;sgn+=2) 
        for (double xx=maxx;xx<=10.;xx*=xfactor) {
            exact=sgn*xx;
            for (i=ib*blength;i>=(ib-1)*blength+1;i--) {
                exact=vertical_slit_half_rad_zip(dt[i],dx[i],exact);
                if (exact.imag()>1.e-8*xx) break;
            }
            if (0) if (exact.imag()>1.e-8*xx && fabs(xx)>maxx)  // delete me 
              printf("maxx changed %lf -> %lf \n",maxx,fabs(xx)); // delete me
            if (exact.imag()>1.e-8*xx && fabs(xx)>maxx) maxx=fabs(xx); 
	    if (0) printf("%lf -> %le + %le i maxx=%lf\n",        // delete me 
		   sgn*xx,exact.real(),exact.imag(),maxx); // delete me 
        } 
        maxx*=xfactor;
        radius [ib]=maxx;

if (0){
	printf(">>>>> ib=%ld radius=%lf\n",ib,radius[ib]);       // delete me 
	printf("dt=%le \n",dt[ib*blength]);
	printf(">>>>> ib=%ld radius=%lf\n",ib,radius[ib]);       // delete me 
	printf(">>>>> ib=%ld radius=%lf\n",ib,radius[ib]);       // delete me 
}
    }
}

void complex_walk::compute_block_rad_zip(long ib)
{
    printf("compute_block_rad_zip not implemented\n");
    exit(0);
}

void complex_walk::compute_block_di_zip(long ib)
{
    printf("compute_block_di_zip not implemented\n");
    exit(0);
}

void complex_walk::compute_parms_ch_zip(long istep,double ddt,double ddx)
{
    double alp,kappa_eff;
    switch(map_choice) {
        case 1: // vertical slit - no parms to be computed
        break;
        case 2:
            kappa_eff=ddx*ddx/ddt; 
            if (ddx>0.) alp=0.5-0.5*sqrt(kappa_eff/(16+kappa_eff));
            else        alp=0.5+0.5*sqrt(kappa_eff/(16+kappa_eff));
            alpha[istep]=alp;
            xl[istep]=2*sqrt(ddt*(1-alp)/alp);
            xr[istep]=2*sqrt(ddt*alp/(1-alp));
        break;
        case 3: // dt=(a^2+2b^2)/8  dx=3a/2
            aa[istep]=2.*ddx/3.;
            bb[istep]=sqrt(4*ddt-aa[istep]*aa[istep]/2.);
        break;
    }
}

void complex_walk::compute_parms_rad_zip(long is,double ddt,double ddx)
{
    double alp,kappa_eff;
    complex<double> z;
    z=I;

    switch(map_choice) {
        case 1:
        break;
        case 2:
            kappa_eff=ddx*ddx/ddt; 
            if (ddx>0.) alp=0.5-0.5*sqrt(kappa_eff/(16+kappa_eff));
            else alp = 0.5 + 0.5*sqrt(kappa_eff/(16+kappa_eff));
            alpha[is]=alp;
            xl[is]=sqrt(ddt*(1-alp)/alp); // factor of 2 wrt chordal ?
            xr[is]=sqrt(ddt*alp/(1-alp)); // factor of 2 wrt chordal ?
	    z=tilted_slit_ch_zip(alpha[is],xl[is],xr[is],z);
            ww[is]=z;
        break;
        case 3: 
            printf("case 3 not implemented in compute_parms_rad_zip\n");
            exit(0);
        break;
    }
}

void complex_walk::compute_parms_ch_unzip(long ipt,complex<double> z)
{
    double r,alp,dcap=0.,ddrive=0.,calpha,x,y;
    switch(map_choice) {
        case 1: // vertical slit 
            ddrive=z.real();
            y=z.imag();
	    dcap=y*y/2.;
        break;
        case 2: // tilted slit
            r=abs(z);
            alp=atan2(z.imag(),z.real())/M_PI;
            alpha[ipt]=alp;
            xl[ipt]=r*pow((1-alp)/alp,alp);
            xr[ipt]=r*pow(alp/(1-alp),1-alp);
            dcap=0.5*r*r*pow(alp,1-2*alp)*pow(1-alp,2*alp-1);
            // dcap is 2*t, dx=calpha*sqrt(t)
            // dx=calpha*sqrt(dcap/2.)
            calpha=2*(1-2*alp)/sqrt(alp*(1-alp));  
            ddrive=calpha*sqrt(dcap/2.);
            if (alp<1.e-10 || alp>1.-1.e-10) {
                printf("%ld alp=%le z=%le+%lei\n",ipt,alp,z.real(),z.imag());
                exit(0);
            }
        break;
        case 3: // arc
            x=z.real();
            y=z.imag();
            aa[ipt]=x;
            bb[ipt]=y;
            dcap=(x*x+2.*y*y)/4.;
            ddrive=1.5*x; 
        break;
    }
    dt[ipt]=dcap/2.;
    dx[ipt]=ddrive;
}

void complex_walk::compute_parms_di_zip(long istep,double ddt,double ddx)
{
    switch(map_choice) {
        case 1: // vertical slit - no parms to be computed
        break;
        default:
            printf("bad map_choice in compute_parms_di_zip \n");
            exit(0);
        break;
    }
}


void complex_walk::ch_zip(long ns,long ds,long mc,long bl,long nt,double lf)
//
// This function uses laurent series. Set blength=0 to prevent their use, i.e
// do "exact" computaion.
//
// Class must already contain dt and dx. We compute the curve that 
// Loewner eq. produces.
{
    long istep,nb,i;
    complex<double> z;

    sle_domain=1;     // needed ??
    nterms=nt;
    laurent_factor=lf;

    // check that mem allocation was consistent with ds,mc,bl 
    if (ds!=dsteps || mc!=map_choice ||bl!=blength) { 
        printf("dsteps/map_choice/blength err in ch_zip\n");
        exit(1);     
    }
    if (laurent_factor<1. && blength>0) {
        printf("laurent_factor <1,chordal_sle_laurent()\n");
        exit(1);
    }

    // Iterate the conformal maps, using the blocks when possible
    t[0]=0.;
    z=0.;
    steps[0]=z;
    for (nsteps=1;nsteps<=ns;nsteps++) {
        if (nsteps%1000==0) printf("nsteps=%ld/%ld, t=%14.10lf, z=%lf+%lf i\n",
	      nsteps,ns,t[nsteps-1],steps[nsteps-1].real(),
              steps[nsteps-1].imag());
        if (blength==0) nb=0; // blength==0 means don't use laurent
        else nb=(nsteps-1)*dsteps/blength;
        if (nb<0) nb=0;
        for (istep=(nsteps-1)*dsteps+1;istep<=nsteps*dsteps;istep++) 
            compute_parms_ch_zip(istep,dt[istep],dx[istep]);
        z=0.;
        z=conformal_map_ch_zip(nsteps,0L,z);

        steps[nsteps]=z;
        t[nsteps]=t[nsteps-1];
        for (i=(nsteps-1)*dsteps+1;i<=nsteps*dsteps;i++) t[nsteps]+=dt[i];

        for (istep=(nsteps-1)*dsteps+1;istep<=nsteps*dsteps;istep++) 
           if (blength>0 && istep%blength==0) 
               compute_block_ch_zip(istep/blength);
    } // end loop on nsteps

} 

void complex_walk::adapt_chordal_sle_laurent(long sd, 
    complex<double> z1, complex<double> z2, long ns, long ds,
    long mc, long bm_approx, double kap, double drft, 
    long stop_choice, double stop_max, 
    double nu, double min_dt, double max_dt, double max_dx, 
    long bl,long nt,double lf)
//
// This function uses laurent series. Set blength=0 to prevent their use, i.e
// do "exact" computaion.
//
// value of bm_approx determines whether we use a deterministic set 
// of dt's or choose them adaptively. 
//
// For bm_approx=1,2,3,4 we use deterministic dt[]. 
// The SLE will have nsteps.
// min_dt, max_dt, max_dx are not used.  
// bm_approx==1 equal dt, bernoulli steps
// bm_approx==2 dt as d(t^1.5), bernoulli steps
// bm_approx==3 equal dt, normal distribution for steps 
// bm_approx==4 dt as d(t^1.5), normal distribution for steps 

// For bm_approx=10 we choose dt[] in an adaptive way to try to make 
// points on SLE roughly equally spaced. The number of steps in the SLE
// is not predetermined, so the values passed in nsteps is not used. 
// Each dt is less than max_dt and (hopefully) small enough that the 
// resulting step in the SLE has distance less than max_dx. But we 
// don't divide dt below min_dt. 
//
// NB: bm_approx=20,21,22,23 have been deleted. See version57 to recover them.
{
    long istep,nb,i;
    double next_dt,delta,var;
    complex<double> z,zdif;
    // SLE(kappa,rho) uses drive for running value of driving function W_t
    double drive=0.;

    sle_domain=sd;    
    nterms=nt;
    laurent_factor=lf;
    kappa=kap;
    drift=drft;

    // check that mem allocation was consistent with ds,mc,bl 
    if (ds!=dsteps || mc!=map_choice ||bl!=blength) { 
        printf("dsteps/map_choice/blength err in adapt_chordal_sle_laurent\n");
        exit(1);     
    }
    if (laurent_factor<1. && blength>0) {
        printf("laurent_factor <1,chordal_sle_laurent()\n");
        exit(1);
    }

    if (bm_approx==300 && dsteps>1) {
        printf("dsteps>1 not implemented for SLE(kappa,rho)\n");
        exit(1);
    }

    nsteps=ns; // needed by deterministic_dt
    switch(bm_approx) {
        // ==1 equal dt, bernoulli steps
        case 1: bm_deterministic_dt(1L,1L); break;
        // ==2 dt as d(t^1.5), bernoulli steps
        case 2: bm_deterministic_dt(2L,1L); break;
        // ==3 equal dt, normal distribution for steps 
        case 3: bm_deterministic_dt(1L,2L); break;
        // ==4 dt as d(t^1.5), normal distribution for steps 
        case 4: bm_deterministic_dt(2L,2L); break;
        case 10: break;
        case 100: cauchy_deterministic_dt(); break;
        case 200: bm_drift_deterministic_dt(); break;
        case 300: break;
    }

    /////////////////////////////////////////////////////////////
    // Iterate the conformal maps, using the blocks when possible
    /////////////////////////////////////////////////////////////
    var=0.;
    t[0]=0.;
    z=0.;
    assign_step(0L,z);
    long count=0;
    next_dt=max_dt;
    for (nsteps=1;;nsteps++) {
        if (nsteps%1000==0)  
            printf(
              "nsteps=%ld/%ld, t=%14.10lf, var=%lf kappa=%6.4lf z=%lf+%lf i\n",
	      nsteps,ns,t[nsteps-1],var,kappa,steps[nsteps-1].real(),
              steps[nsteps-1].imag());
        if (blength==0) nb=0; // blength==0 means don't use laurent
        else nb=(nsteps-1)*dsteps/blength;
        if (nb<0) nb=0;
        switch(bm_approx) {
            case 1:       
            case 2:       
            case 3:       
            case 4:       
            case 100:       
            case 200:       
                for (istep=(nsteps-1)*dsteps+1;istep<=nsteps*dsteps;istep++) 
                    compute_parms_ch_zip(istep,dt[istep],dx[istep]);
                count++;
                z=0.;
                z=conformal_map_ch_zip(nsteps,0L,z);
            break;            
            case 10:       
                next_dt=max_dt;
                do {
                    double old_next_dt=next_dt;
                    sample_bm(max_dt,&next_dt);
                    // NB sample_bm can change next_dt
                   for (istep=(nsteps-1)*dsteps+1;istep<=nsteps*dsteps;istep++)
                        compute_parms_ch_zip(istep,dt[istep],dx[istep]);
                    count++;
                    z=0.;
                    z=conformal_map_ch_zip(nsteps,0L,z);
                    if (sle_domain==3) z=half_to_strip(z,z1);
                    delta=abs(z-steps[nsteps-1]);
                    // only half next_dt if it was not changed by sample_bm
                    if (fabs(1.-next_dt/old_next_dt)<1.e-8) next_dt/=2.;
                } while (delta>max_dx && next_dt>0.99*min_dt);
            break;            
            case 300:       
                next_dt=max_dt;
                do {
                    double old_next_dt=next_dt;
                    sample_bm(max_dt,&next_dt);
                    // NB sample_bm can change next_dt
                    // Add contribution from force pts, and evolve force pts
                    for (i=0;i<nforce_pt;i++) {
                        dx[nsteps]+= rho[i]*dt[nsteps]/(drive-force_pt[i]);
                        force_pt[i]+= 2*dt[nsteps]/(force_pt[i]-drive);
                        if (nsteps%100==0) 
                            printf("force_pt[%1ld]=%lf\n",i,force_pt[i]);
                    }
                    drive+=dx[nsteps];
                    compute_parms_ch_zip(nsteps,dt[nsteps],dx[nsteps]);
                    count++;
                    z=0.;
                    z=conformal_map_ch_zip(nsteps,0L,z);
                    if (sle_domain==3) z=half_to_strip(z,z1);
                    delta=abs(z-steps[nsteps-1]);
                    // only half next_dt if it was not changed by sample_bm
                    if (fabs(1.-next_dt/old_next_dt)<1.e-8) next_dt/=2.;
                } while (delta>max_dx && next_dt>0.99*min_dt);
            break;            
        }

        // NB for sle_domain==3, z has already been mapped to the strip
        assign_step(nsteps,z);
        zdif=z-I;

        t[nsteps]=t[nsteps-1];
        for (i=(nsteps-1)*dsteps+1;i<=nsteps*dsteps;i++) t[nsteps]+=dt[i];

        // compute increment for fractional variation
        var+=pow(abs(steps[nsteps]-steps[nsteps-1]),1./nu);

        if (bm_nsteps>dsteps*max_nsteps) break;
        double stop_var=0.; 
        switch (stop_choice) {
            case 1: stop_var=t[nsteps]; break;
            case 2: stop_var=var; break;
            case 3: stop_var=abs(steps[nsteps]); break;
	}
        if (stop_choice>0 && stop_var>stop_max) break;

        // deterministic dt should stop at given number of steps
        if (bm_approx<10 && nsteps>=ns) break;
        if ((bm_approx==100 || bm_approx==200 || bm_approx>300) 
            && nsteps>=ns) break;

        // stop adaptive SLE if it gets close to finite terminal point
	// NB : this is not always correct
        // For example, sle_domain=1 and terminal pt not infinity
	// NB : this is not always correct
	// NB : this is not always correct
	// NB : this is not always correct
	// NB : this is not always correct
	// NB : this is not always correct
	// NB : this is not always correct
	// NB : this is not always correct
	// NB : this is not always correct
	// NB : this is not always correct
	// NB : this is not always correct
	// NB : this is not always correct
	// NB : this is not always correct
	// NB : this is not always correct
	// NB : this is not always correct
	// NB : this is not always correct
	// NB : this is not always correct
        if (bm_approx==10 && sle_domain!=1 && abs(zdif)<max_dx) break;

        for (istep=(nsteps-1)*dsteps+1;istep<=nsteps*dsteps;istep++) 
           if (blength>0 && istep%blength==0) 
               compute_block_ch_zip(istep/blength);
    } // end loop on nsteps

	// NB : this is not always correct
	// NB : this is not always correct
	// NB : this is not always correct
	// NB : this is not always correct
    if (sle_domain!=1 && stop_choice==0) { // add terminal point 
        nsteps++;
        z=I;
        assign_step(nsteps,z);
        t[nsteps]=1.e10;
    }

    printf("count=%ld nsteps=%ld  percent kept=%lf var=%lf t=%lf\n",
     count,nsteps,100*double(nsteps)/double(count),var,t[nsteps]);
} // end adapt_chordal_sle_laurent()


void complex_walk::adapt_radial_sle_laurent(long sd, 
    complex<double> z1, complex<double> z2,long ns, long ds,
    long mc, long bm_approx, double kap, long stop_choice, double stop_max, 
    double nu, double min_dt, double max_dt, double max_dx, 
    long bl,long nt,double lf)
//
// This function uses laurent series. Set blength=0 to prevent their use, i.e
// do "exact" computaion.
//
// value of bm_approx determines whether we use a deterministic set 
// of dt's or choose them adaptively. 
//
// For bm_approx=1,2,3,4 we use deterministic dt[]. 
// The SLE will have nsteps.
// min_dt, max_dt, max_dx are not used.  

// For bm_approx=10 we choose dt[] in an adaptive way to try to make 
// points on SLE roughly equally spaced. The number of steps in the SLE
// is not predetermined, so the values passed in nsteps is not used. 
// Each dt is less than max_dt and (hopefully) small enough that the 
// resulting step in the SLE has distance less than max_dx. But we 
// don't divide dt below min_dt. 
//
{
    long istep,nb,i;
    double next_dt,var,delta;
    complex<double> z,zdif,zend;

    sle_domain=sd;    
    nterms=nt;
    laurent_factor=lf;
    kappa=kap;

    zend=0.; // terminal point for radial SLE 
    if (sle_domain==1) zend=disc_to_half(zend,z1);

    // check that mem allocation was consistent with ds,mc,bl 
    if (ds!=dsteps || mc!=map_choice ||bl!=blength) { 
        printf("dsteps/map_choice/blength err in adapt_radial_sle_laurent\n");
        exit(1);     
    }
    if (laurent_factor<1. && blength>0) {
        printf("laurent_factor <1,adapt_radial_sle_laurent()\n");
        exit(1);
    }

    nsteps=ns; // needed by deterministic_dt
    switch(bm_approx) {
        // ==1 equal dt, bernoulli steps
        case 1: bm_deterministic_dt(1L,1L); break;
        // ==2 dt as d(t^1.5), bernoulli steps
        case 2: bm_deterministic_dt(2L,1L); break;
        // ==3 equal dt, normal distribution for steps 
        case 3: bm_deterministic_dt(1L,2L); break;
        // ==4 dt as d(t^1.5), normal distribution for steps 
        case 4: bm_deterministic_dt(2L,2L); break;
        case 10: break;
    }

    /////////////////////////////////////////////////////////////
    // Iterate the conformal maps, using the blocks when possible
    /////////////////////////////////////////////////////////////
    var=0.;
    t[0]=0.;
    z=1.;
    if (sle_domain==1) z=disc_to_half(z,z1);
    assign_step(0L,z);
    long count=0;
    next_dt=max_dt;
    for (nsteps=1;;nsteps++) {
        if (nsteps%100==0) 
          printf(
           "ns=%ld/%ld, %6.2le t=%13.10lf, var=%lf kap=%6.2lf z-zend=%lf RAD\n",
            nsteps,ns,next_dt,t[nsteps-1],var,kappa,abs(steps[nsteps-1]-zend));
        if (blength==0) nb=0; // blength==0 means don't use laurent
        else nb=(nsteps-1)*dsteps/blength;
        if (nb<0) nb=0;
        switch(bm_approx) {
            case 1:       
            case 2:       
            case 3:       
            case 4:       
                for (istep=(nsteps-1)*dsteps+1;istep<=nsteps*dsteps;istep++) 
                    compute_parms_rad_zip(istep,dt[istep],dx[istep]);
                count++;
		// kludge to avoid numerical err pushing you over branch cut
                z=1.-1.e-10;
                z=conformal_map_rad_zip(nsteps,z);
                if (sle_domain==1) z=disc_to_half(z,z1);
            break;            
            case 10:       
                next_dt=max_dt;
                do {
                    double old_next_dt=next_dt;
                    sample_bm(max_dt,&next_dt);
                    // NB sample_bm can change next_dt
                   for (istep=(nsteps-1)*dsteps+1;istep<=nsteps*dsteps;istep++)
                        compute_parms_rad_zip(istep,dt[istep],dx[istep]);
                    count++;
                    z=1.-1.e-10;
                    z=conformal_map_rad_zip(nsteps,z);
                    if (sle_domain==1) z=disc_to_half(z,z1);
                    delta=abs(z-steps[nsteps-1]);
                    // only half next_dt if it was not changed by sample_bm
                    if (fabs(1.-next_dt/old_next_dt)<1.e-8) next_dt/=2.;
                } while (delta>max_dx && next_dt>0.99*min_dt);
            break;            
        }

        assign_step(nsteps,z);
        t[nsteps]=t[nsteps-1];
        for (i=(nsteps-1)*dsteps+1;i<=nsteps*dsteps;i++) t[nsteps]+=dt[i];

        // compute increment for fractional variation
        var+=pow(abs(steps[nsteps]-steps[nsteps-1]),1./nu);

        // was >= ??????????????????????????????????
        if (bm_nsteps>max_nsteps) break;

        double stop_var=0.; 
        switch (stop_choice) {
            case 1: stop_var=t[nsteps]; break;
            case 2: stop_var=var; break;
            case 3: stop_var=abs(steps[nsteps]); break;
	}
        if (stop_choice>0 && stop_var>stop_max) break;

        // deterministic dt should stop at given number of steps
        if (bm_approx<10 && nsteps>=ns) break;

        // stop adaptive SLE if it gets close to terminal point
	if ((bm_approx>=10 && bm_approx<30) && abs(z-zend)<max_dx) break;

        for (istep=(nsteps-1)*dsteps+1;istep<=nsteps*dsteps;istep++) 
           if (blength>0 && istep%blength==0) 
               compute_block_rad_zip(istep/blength);
    } // end loop on nsteps

    if (stop_choice==0) { // add terminal point 
        nsteps++;
        assign_step(nsteps,zend);
        t[nsteps]=1.e10;
    }

    printf("count=%ld nsteps=%ld  percent kept=%lf var=%lf t=%lf\n",
     count,nsteps,100*double(nsteps)/double(count),var,t[nsteps]);
} // end adapt_radial_sle_laurent()

void complex_walk::adapt_half_radial_sle_laurent(long sd, 
    complex<double> z1, complex<double> z2,long ns, long ds,
    long mc, long bm_approx, double kap, long stop_choice, double stop_max, 
    double nu, double min_dt, double max_dt, double max_dx, 
    long bl,long nt,double lf)
//
// This function uses laurent series. Set blength=0 to prevent their use, i.e
// do "exact" computaion.
//
// value of bm_approx determines whether we use a deterministic set 
// of dt's or choose them adaptively. 
//
// For bm_approx=1,2,3,4 we use deterministic dt[]. 
// The SLE will have nsteps.
// min_dt, max_dt, max_dx are not used.  

// For bm_approx=10 we choose dt[] in an adaptive way to try to make 
// points on SLE roughly equally spaced. The number of steps in the SLE
// is not predetermined, so the values passed in nsteps is not used. 
// Each dt is less than max_dt and (hopefully) small enough that the 
// resulting step in the SLE has distance less than max_dx. But we 
// don't divide dt below min_dt. 
//
{
    long istep,nb,i;
    double next_dt,var,delta;
    complex<double> z,zdif,zend;

    sle_domain=sd;    
    nterms=nt;
    laurent_factor=lf;
    kappa=kap;

    zend=I; // terminal point for radial SLE in half plane
    if (sle_domain!=1) {
        printf("sle_domain!=1 in adapt_half_radial_sle_laurent()\n"); 
        exit(0);
    }

    // check that mem allocation was consistent with ds,mc,bl 
    if (ds!=dsteps || mc!=map_choice ||bl!=blength) { 
        printf("dsteps/map_choice/blength er:adapt_half_radial_sle_laurent\n");
        exit(1);     
    }
    if (laurent_factor<1. && blength>0) {
        printf("laurent_factor <1,adapt_half_radial_sle_laurent()\n");
        exit(1);
    }

    nsteps=ns; // needed by deterministic_dt
    switch(bm_approx) {
        // ==1 equal dt, bernoulli steps
        case 1: bm_deterministic_dt(1L,1L); break;
        // ==2 dt as d(t^1.5), bernoulli steps
        case 2: bm_deterministic_dt(2L,1L); break;
        // ==3 equal dt, normal distribution for steps 
        case 3: bm_deterministic_dt(1L,2L); break;
        // ==4 dt as d(t^1.5), normal distribution for steps 
        case 4: bm_deterministic_dt(2L,2L); break;
        case 10: break;
    }

    /////////////////////////////////////////////////////////////
    // Iterate the conformal maps, using the blocks when possible
    /////////////////////////////////////////////////////////////
    var=0.;
    t[0]=0.;
    z=0.;
    // for other domains need to map z here  <<<<<<<<<<<<<<<<<<<<<<<
    assign_step(0L,z);
    long count=0;
    next_dt=max_dt;
    for (nsteps=1;;nsteps++) {
        if (nsteps%100==0) 
          printf(
          "ns=%ld/%ld, %6.2le t=%13.10lf, var=%lf kap=%6.2lf z-zend=%lf RAD\n",
            nsteps,ns,next_dt,t[nsteps-1],var,kappa,abs(steps[nsteps-1]-zend));
        if (blength==0) nb=0; // blength==0 means don't use laurent
        else nb=(nsteps-1)*dsteps/blength;
        if (nb<0) nb=0;
        switch(bm_approx) {
            case 1:       
            case 2:       
            case 3:       
            case 4:       
                for (istep=(nsteps-1)*dsteps+1;istep<=nsteps*dsteps;istep++) 
                    compute_parms_rad_zip(istep,dt[istep],dx[istep]);
                count++;
		// kludge to avoid numerical err pushing you over branch cut
                z=1.e-10*I;
                z=conformal_map_half_rad_zip(nsteps,z);
                // for other domains need to map z here  <<<<<<<<<<<<<<<<<<<<<
            break;            
            case 10:       
                next_dt=max_dt;
                do {
                    double old_next_dt=next_dt;
                    sample_bm(max_dt,&next_dt);
                    // NB sample_bm can change next_dt
                   for (istep=(nsteps-1)*dsteps+1;istep<=nsteps*dsteps;istep++)
                        compute_parms_rad_zip(istep,dt[istep],dx[istep]);
                    count++;
                    z=1.e-10*I;
                    z=conformal_map_half_rad_zip(nsteps,z);
                    // for other domains need to map z here  <<<<<<<<<<<<<<<<<<
                    delta=abs(z-steps[nsteps-1]);
                    // only half next_dt if it was not changed by sample_bm
                    if (fabs(1.-next_dt/old_next_dt)<1.e-8) next_dt/=2.;
                } while (delta>max_dx && next_dt>0.99*min_dt);
            break;            
        }

        assign_step(nsteps,z);
        t[nsteps]=t[nsteps-1];
        for (i=(nsteps-1)*dsteps+1;i<=nsteps*dsteps;i++) t[nsteps]+=dt[i];

        // compute increment for fractional variation
        var+=pow(abs(steps[nsteps]-steps[nsteps-1]),1./nu);

        // was >= ??????????????????????????????????
        if (bm_nsteps>max_nsteps) break;

        double stop_var=0.; 
        switch (stop_choice) {
            case 1: stop_var=t[nsteps]; break;
            case 2: stop_var=var; break;
            case 3: stop_var=abs(steps[nsteps]); break;
	}
        if (stop_choice>0 && stop_var>stop_max) break;

        // deterministic dt should stop at given number of steps
        if (bm_approx<10 && nsteps>=ns) break;

        // stop adaptive SLE if it gets close to terminal point
	if ((bm_approx>=10 && bm_approx<30) && abs(z-zend)<max_dx) break;

        // always stop if number of steps exceeds max_nsteps
        // this is useful for giving up when algorithm gets stuck
        if (nsteps>=max_nsteps) break;

        for (istep=(nsteps-1)*dsteps+1;istep<=nsteps*dsteps;istep++) 
           if (blength>0 && istep%blength==0) 
               compute_block_half_rad_zip(istep/blength);
    } // end loop on nsteps

    if (stop_choice==0) { // add terminal point 
        nsteps++;
        assign_step(nsteps,zend);
        t[nsteps]=1.e10;
    }

    printf("count=%ld nsteps=%ld  percent kept=%lf var=%lf t=%lf\n",
     count,nsteps,100*double(nsteps)/double(count),var,t[nsteps]);
} // end adapt_half_radial_sle_laurent()

void complex_walk::adapt_dipolar_sle_laurent(long sd, 
    complex<double> z1, complex<double> z2,long ns, long ds,
    long mc, long bm_approx, double kap, 
    long stop_choice, double stop_max, 
    double nu, double min_dt, double max_dt, double max_dx, 
    long bl,long nt,double lf)
//
// This function uses laurent series. Set blength=0 to prevent their use, i.e
// do "exact" computaion.
//
// value of bm_approx determines whether we use a deterministic set 
// of dt's or choose them adaptively. 
// Values are same as for adapt_chordal_sle_laurent, except we only 
// implement bm_approx=1,2,3,4,10
//
{
    long istep,nb,i;
    double next_dt,delta,var,dif;
    complex<double> z;
    sle_domain=sd;    
    nterms=nt;
    laurent_factor=lf;
    kappa=kap;

    // check that mem allocation was consistent with ds,mc,bl 
    if (ds!=dsteps || mc!=map_choice ||bl!=blength) { 
        printf("dsteps/map_choice/blength err in adapt_chordal_sle_laurent\n");
        exit(1);     
    }
    if (laurent_factor<1. && blength>0) {
        printf("laurent_factor <1,chordal_sle_laurent()\n");
        exit(1);
    }

    printf("nsteps=%ld\n",nsteps);
    nsteps=ns; // needed by deterministic_dt
    switch(bm_approx) {
        // ==1 equal dt, bernoulli steps
        case 1: bm_deterministic_dt(1L,1L); break;
        // ==2 dt as d(t^1.5), bernoulli steps
        case 2: bm_deterministic_dt(2L,1L); break;
        // ==3 equal dt, normal distribution for steps 
        case 3: bm_deterministic_dt(1L,2L); break;
        // ==4 dt as d(t^1.5), normal distribution for steps 
        case 4: bm_deterministic_dt(2L,2L); break;
        case 10: break;
    }

    /////////////////////////////////////////////////////////////
    // Iterate the conformal maps, using the blocks when possible
    /////////////////////////////////////////////////////////////
    var=0.;
    t[0]=0.;
    z=0.;
    assign_step(0L,z);
    long count=0;
    next_dt=max_dt;
    for (nsteps=1;;nsteps++) {
        if (nsteps%100==0) 
            printf(
              "nsteps=%ld/%ld, t=%14.10lf, var=%lf kappa=%6.4lf z=%lf+%lf i\n",
	      nsteps,ns,t[nsteps-1],var,kappa,steps[nsteps-1].real(),
              steps[nsteps-1].imag());
        if (blength==0) nb=0; // blength==0 means don't use laurent
        else nb=(nsteps-1)*dsteps/blength;
        if (nb<0) nb=0;
        switch(bm_approx) {
            case 1:       
            case 2:       
            case 3:       
            case 4:       
                for (istep=(nsteps-1)*dsteps+1;istep<=nsteps*dsteps;istep++) 
                    compute_parms_di_zip(istep,dt[istep],dx[istep]);
                count++;
                z=0.;
                z=conformal_map_di_zip(nsteps,0L,z);
printf("conformal_map_di_zip (nsteps=%ld) => z=%lf+%lf i\n",
nsteps,z.real(),z.imag());
            break;            
            case 10:       
                next_dt=max_dt;
                do {
                    double old_next_dt=next_dt;
                    sample_bm(max_dt,&next_dt);
                    // NB sample_bm can change next_dt
                   for (istep=(nsteps-1)*dsteps+1;istep<=nsteps*dsteps;istep++)
                        compute_parms_di_zip(istep,dt[istep],dx[istep]);
                    count++;
                    z=0.;
                    z=conformal_map_di_zip(nsteps,0L,z);
                    delta=abs(z-steps[nsteps-1]);
                    // only half next_dt if it was not changed by sample_bm
                    if (fabs(1.-next_dt/old_next_dt)<1.e-8) next_dt/=2.;
                } while (delta>max_dx && next_dt>0.99*min_dt);
            break;            
        }

        assign_step(nsteps,z);
        dif=M_PI-z.imag();

        t[nsteps]=t[nsteps-1];
        for (i=(nsteps-1)*dsteps+1;i<=nsteps*dsteps;i++) t[nsteps]+=dt[i];

        // compute increment for fractional variation
        var+=pow(abs(steps[nsteps]-steps[nsteps-1]),1./nu);

        if (bm_nsteps>=max_nsteps) break;
        double stop_var=0.; 
        switch (stop_choice) {
            case 1: stop_var=t[nsteps]; break;
            case 2: stop_var=var; break;
            case 3: stop_var=steps[nsteps].imag()/M_PI; break;
	}
        if (stop_choice>0 && stop_var>stop_max) break;

        // deterministic dt should stop at given number of steps
        if (bm_approx<10 && nsteps>=ns) break;

        // NB : may not alway be correct
        // NB : may not alway be correct
        // NB : may not alway be correct
        // NB : may not alway be correct
        // NB : may not alway be correct
        // stop adaptive SLE if it gets close to finite terminal point
	if (bm_approx==10 && dif<max_dx) break;

        for (istep=(nsteps-1)*dsteps+1;istep<=nsteps*dsteps;istep++) 
           if (blength>0 && istep%blength==0) 
               compute_block_di_zip(istep/blength);
    } // end loop on nsteps

    printf("count=%ld nsteps=%ld  percent kept=%lf var=%lf t=%lf\n",
     count,nsteps,100*double(nsteps)/double(count),var,t[nsteps]);
} // end adapt_dipolar_sle_laurent()

void complex_walk::ch_unzip(long npt,long bl,long nt,double lf,long mc,
    double max_t) 
// Given a walk in the complex plane, this routine finds a sequence of 
// simple conformal maps that unzips its. 
// map_choice determines what simple conformal map is used:
// The parameters that represent these maps are stored in the class 
//  map_choice=1: vertical slit; parameters none
//  map_choice=2: tilted slit; parameters xl,xr,alpha
//  map_choice=3: circular arc; aa,bb
// The increments of the driving function, dt[],dx[] are also put in the class 
// 
// Stopping: we only unzip up to at most point npt which can be less than 
// nsteps. Furthermore, we stop if the time (cap/2) exceed max_t.
{
    complex<double> z,ii;
    long ib,nb,ipt;
    laurent ftemp;
    double rad,t=0.;
    ii=I;

    // KLUDGE: mc,nt,bl,lf should not be parameters in the function
    map_choice=mc;    
    nterms=nt;
    blength=bl;
    laurent_factor=lf;

    for (ib=1;ib<=nblock;ib++) radius[ib]=0.;
    // the timing lines should be uncommented only to run timing tests, e.g.,
    // run_loewner_timing
    double etime=-mytime();
    FILE *fptr; 
    fptr=fopen("time.plt","w");
    fclose(fptr); 

    for (ipt=1;ipt<=npt;ipt++) {
        // the timing lines should be uncommented only to run timing tests
        if (ipt%1000==0) {
             printf("ipt=%ld time=%lf \n",ipt,etime+mytime());
             fptr=fopen("time.plt","a");
             fprintf(fptr,"%ld %lf \n",ipt,etime+mytime());
             fclose(fptr);
        }    
        z=steps[ipt];
        // apply conformal maps 1 to (ipt-1); z is now next point to map to 0
        z=conformal_map_ch_unzip(ipt-1,&rad,z);

        // determine radius: it is maximal distance from "curve" to origin
        // curve for radius[nb+1] is image of points nb*blength+1 to
        // nb*(blength+1) after applying block maps 1 to nb
        if (blength>0) {
            nb=(ipt-1)/blength;
            if (rad>radius[nb+1]) radius[nb+1]=rad;
        }

        compute_parms_ch_unzip(ipt,z);
 
        // compute inverted laurent series for block ib
        if (blength>0) {
            ib=ipt/blength;
            if (ipt%blength==0) compute_block_ch_unzip(ib);
        }
       
        t+=dt[ipt];
        if (t>max_t) {
            num_unzip=ipt;
            return;
        }
    } // end loop on ipt
    num_unzip=npt;
} 

// ???????????????????????????? location ?????????????????????????
void complex_walk::ch_sewing(long npts)
{
    for (long ipt=0;ipt<=npts;ipt++) {
       conformal_map_ch_sewing(ipt,npts,sewing_minus+ipt,sewing_plus+ipt);
    }
}

complex<double> complex_walk::conformal_map_ch_zip(long ns, long imag_flag, 
    complex<double> z)
//  A complex walk defines a curve in the half plane. So it determines 
//  a conformal map - the map that takes the half plane onto the half plane 
//  minus the curve. The map is given by a composition of many simple 
//  maps. The parameters that determine these maps are stored in the class. 
//  The routine computes the image of z under the conformal map.
{
    long istep,ib,nb;

    if (blength==0) nb=0; // blength==0 means don't use laurent
    else nb=(ns-1)*dsteps/blength;
    if (nb<0) nb=0;

    // partial block: must use exact computation
    for (istep=ns*dsteps;istep>nb*blength;istep--) {
        switch(map_choice) {
            case 1:
                z=vertical_slit_ch_zip(dt[istep],dx[istep],z);
            break;
            case 2:
                z=tilted_slit_ch_zip(alpha[istep],xl[istep],xr[istep],z);
            break;
            case 3:
                z=arc_ch_zip(aa[istep],bb[istep],z);
            break;
        }
    }

    // full blocks: must decide exact computation vs laurent series
    for (ib=nb;ib>=1;ib--) {
        // test if we can use laurent series
        if (abs(z)>laurent_factor*radius[ib]) {// use laurent
            z=1./fblock[ib].evaluate(1./z);   
            // error from laurent series can result in z in lower half plane
            if (z.imag()<0.) z=z.real();
        }
        else {  // use exact
          for (istep=ib*blength;istep>=(ib-1)*blength+1;istep--) {
            switch(map_choice) {
                case 1:
                    z=vertical_slit_ch_zip(dt[istep],dx[istep],z);
                break;
                case 2:
                    z=tilted_slit_ch_zip(alpha[istep],xl[istep],xr[istep],z);
                break;
                case 3:
                    z=arc_ch_zip(aa[istep],bb[istep],z);
                break;
            }
          }
        } // end else
    } // end loop on ib
    return(z);
}

complex<double> complex_walk::conformal_map_rad_zip(long ns, complex<double> z)
{
    long is,ib,nb;

    if (blength==0) nb=0; // blength==0 means don't use laurent
    else nb=(ns-1)*dsteps/blength;
    if (nb<0) nb=0;

    // partial block: must use exact computation
    for (is=ns*dsteps;is>nb*blength;is--) {
        switch(map_choice) {
            case 1:
                z=vertical_slit_rad_zip(dt[is],dx[is],z);
            break;
            case 2:
                z=tilted_slit_rad_zip(dx[is],alpha[is],xl[is],xr[is],
                    ww[is],z);
            break;
            case 3:
	      // z=arc_rad_zip(aa[is],bb[is],z);
                printf("case 3 not implemented in conformal_map_rad_zip\n"); 
                exit(0); 
            break;
        }
    }

    // full blocks: must decide exact computation vs laurent series
    for (ib=nb;ib>=1;ib--) {
        // test if we can use laurent series
        if (abs(z)>laurent_factor*radius[ib]) {  // use laurent
            z=1./fblock[ib].evaluate(1./z);   
            // error from laurent series can result in z in lower half plane
            if (z.imag()<0.) z=z.real();
        }
        else {  // use exact
          for (is=ib*blength;is>=(ib-1)*blength+1;is--) {
            switch(map_choice) {
                case 1:
                    z=vertical_slit_rad_zip(dt[is],dx[is],z);
                break;
                case 2:
		    z=tilted_slit_rad_zip(dx[is],alpha[is],xl[is],xr[is],
                        ww[is],z);
                break;
                case 3:
		  // z=arc_rad_zip(aa[is],bb[is],z);
                   printf("case 3 not implemnted in conformal_map_rad_zip\n"); 
                   exit(0); 
                break;
            }
          }
        } // end else
    } // end loop on ib
    return(z);
}

complex<double> complex_walk::conformal_map_half_rad_zip(long ns, 
    complex<double> z)
{
    long is,ib,nb;

    if (blength==0) nb=0; // blength==0 means don't use laurent
    else nb=(ns-1)*dsteps/blength;
    if (nb<0) nb=0;

    // partial block: must use exact computation
    for (is=ns*dsteps;is>nb*blength;is--) {
        switch(map_choice) {
            case 1:
                z=vertical_slit_half_rad_zip(dt[is],dx[is],z);
            break;
            case 2:
                z=tilted_slit_half_rad_zip(dx[is],alpha[is],xl[is],xr[is],
                    ww[is],z);
            break;
            case 3:
	      // z=arc_rad_zip(aa[is],bb[is],z);
                printf("case 3 not implemented:conformal_map_half_rad_zip\n"); 
                exit(0); 
            break;
        }
    }

    // full blocks: must decide exact computation vs laurent series
    for (ib=nb;ib>=1;ib--) {


      /*
complex<double> zlaurent;
double ratio; 
ratio=abs(z)/radius[ib];
if (ratio>4.) printf("\n>>>>>>> orig z=%lf + %lf i\n",z.real(),z.imag());
zlaurent=1./fblock[ib].evaluate(1./z);   
// error from laurent series can result in z in lower half plane
if (zlaurent.imag()<0.) zlaurent=z.real();
if (ratio>4.) 
  printf("   zlaurent=%lf + %lf i\n",zlaurent.real(),zlaurent.imag());

for (is=ib*blength;is>=(ib-1)*blength+1;is--) 
                z=vertical_slit_half_rad_zip(dt[is],dx[is],z);
if (ratio>4.) printf("   zexact  =%lf + %lf i\n",z.real(),z.imag());
if (ratio>4.) printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>ratio=%lf err  =%le\n",
		     ratio,abs(z-zlaurent));
      */



        // test if we can use laurent series
        if (abs(z)>laurent_factor*radius[ib]) {  // use laurent
            z=1./fblock[ib].evaluate(1./z);   
            // error from laurent series can result in z in lower half plane
            if (z.imag()<0.) z=z.real();
        }
        else {  // use exact
          for (is=ib*blength;is>=(ib-1)*blength+1;is--) {
            switch(map_choice) {
                case 1:
                    z=vertical_slit_half_rad_zip(dt[is],dx[is],z);
                break;
                case 2:
		    z=tilted_slit_half_rad_zip(dx[is],alpha[is],xl[is],xr[is],
                        ww[is],z);
                break;
                case 3:
		  // z=arc_rad_zip(aa[is],bb[is],z);
                   printf("cse 3 not implemnted:conformal_map_half_rad_zip\n"); 
                   exit(0); 
                break;
            }
          }
        } // end else
    } // end loop on ib
    return(z);
}

complex<double> complex_walk::conformal_map_ch_unzip(long ipt, double *rad,
    complex<double> z)
//  A complex walk defines a curve in the half plane. So it determines 
//  a conformal map - the map that takes the half plane minus the curve onto 
//  the half plane. The map is given by a composition of many simple 
//  maps. The parameters that determine these maps are stored in the class. 
//  The routine computes the image of z under the conformal map.
{
    long ib,nb,i;

    if (blength==0) nb=0;
    else nb=ipt/blength;
    // full blocks: must decide exact computation vs laurent series
    for (ib=1;ib<=nb;ib++) {
        // test if we can use laurent series
        if (abs(z)>laurent_factor*radius[ib]) {  // use laurent
            z=1./fblock[ib].evaluate(1./z);   
            // error from laurent series can result in z in lower half plane
            if (z.imag()<0.) z=z.real();
        }
        else {  // use exact
            switch(map_choice) {
                case 1: // vertical slit
                    for (i=(ib-1)*blength+1;i<=ib*blength;i++) {
                        z=vertical_slit_ch_unzip(dt[i],dx[i],z);
                    }
                break;
                case 2: // tilted slit
                    for (i=(ib-1)*blength+1;i<=ib*blength;i++) {
                        z=tilted_slit_ch_unzip(alpha[i],xl[i],xr[i],z,z);
                    }
                break;
                case 3: // arc
                    for (i=(ib-1)*blength+1;i<=ib*blength;i++) {
                        if (fabs(aa[i])+fabs(bb[i])>1.e-8) 
                            z=arc_ch_unzip(aa[i],bb[i],z);
                        else // KLUDGE to avoid calling with aa=bb=0
                            z=z;
                    }
                break;
            }
        } // end else
    } // end loop on ib

    // radius is maximal distance from "curve" to origin
    // curve for radius[nb+1] is image of points nb*blength+1 to
    // nb*(blength+1) after applying block maps 1 to nb
    // So we need to return this distance
    *rad=abs(z);

    // partial block: must use exact computation
    for (i=nb*blength+1;i<=ipt;i++) {
        switch(map_choice) {
            case 1: // vertical slit
                z=vertical_slit_ch_unzip(dt[i],dx[i],z);
            break;
            case 2: // tilted slit
                z=tilted_slit_ch_unzip(alpha[i],xl[i],xr[i],z,z);
            break;
            case 3: // arc
                if (fabs(aa[i])+fabs(bb[i])>1.e-8) 
                    z=arc_ch_unzip(aa[i],bb[i],z);
                else // KLUDGE to avoid calling with aa=bb=0
                    z=z;
            break;
        }
    }
    return(z);
}

void complex_walk::conformal_map_ch_sewing(long ipt, long npts,
    double *xm,double *xp)
// Suppose curve from points 0 to npts is fully unzipped. 
// Then the point ipt on the curve is mapped to points on the real axis:
//   *xm and *xp
// Conformal maps 1 to ipt take point ipt to the origin.
// Then conformal map ipt+1 splits 0+ and 0- 
// Then conformal maps ipt+2 to npts map them out onto pos/neg real axis
// Present version does not use laurent series speed up.
{
    // always use exact computation 
    switch(map_choice) {
        case 1: // vertical slit
            printf("vertical slit not allowed for sewing\n");
            exit(1);
        break;
        case 2: // tilted slit
            tilted_slit_ch_unzip_origin(alpha[ipt+1],xl[ipt+1],xr[ipt+1],
                xm,xp);
        break;
        case 3: // arc
            if (fabs(aa[ipt+1])+fabs(bb[ipt+1])>1.e-8) 
                arc_ch_unzip_origin(aa[ipt+1],bb[ipt+1],xm,xp);
            else {// KLUDGE to avoid calling with aa=bb=0
                *xm=-1.e-8;  // KLUDGE
                *xp=1.e-8;  // KLUDGE
            }
        break;
    }

    for (long i=ipt+2;i<=npts;i++) {
        switch(map_choice) {
            case 2: // tilted slit
                *xm=tilted_slit_ch_unzip_real(alpha[i],xl[i],xr[i],*xm);
                *xp=tilted_slit_ch_unzip_real(alpha[i],xl[i],xr[i],*xp);
            break;
            case 3: // arc
                // KLUDGE to avoid calling with aa=bb=0
                if (fabs(aa[i])+fabs(bb[i])>1.e-8) {
                    *xm=arc_ch_unzip_real(aa[i],bb[i],*xm);
                    *xp=arc_ch_unzip_real(aa[i],bb[i],*xp);
                }
            break;
        }
    }

 printf("conformal_map_ch_sewing, ipt=%ld/%ld xm=%lf xp=%lf\n",
	ipt,npts,*xm,*xp);
}

complex<double> complex_walk::conformal_map_di_zip(long ns, long imag_flag, 
    complex<double> z)
{
    long istep,ib,nb;

    if (blength==0) nb=0; // blength==0 means don't use laurent
    else nb=(ns-1)*dsteps/blength;
    if (nb<0) nb=0;

    // partial block: must use exact computation
    for (istep=ns*dsteps;istep>nb*blength;istep--) {
        switch(map_choice) {
            case 1:
                z=vertical_slit_di_zip(dt[istep],dx[istep],z);
            break;
        }
    }

    // full blocks: must decide exact computation vs laurent series
    for (ib=nb;ib>=1;ib--) {
        // test if we can use laurent series
        if (abs(z)>laurent_factor*radius[ib]) {// use laurent
            z=1./fblock[ib].evaluate(1./z);   
            // error from laurent series can result in z in lower half plane
            if (z.imag()<0.) z=z.real();
        }
        else {  // use exact
          for (istep=ib*blength;istep>=(ib-1)*blength+1;istep--) {
            switch(map_choice) {
                case 1:
                    z=vertical_slit_di_zip(dt[istep],dx[istep],z);
                break;
            }
          }
        } // end else
    } // end loop on ib
    return(z);
}

double complex_walk::max_radius()
// computes max of |z| along the walk
{
    double max=0.,temp;
    for (long i=0;i<=nsteps;i++) {
        temp=abs(steps[i]);
        if (temp>max) max=temp; 
    }
    return(max);
}


complex<double> segment_intersect(complex<double> z1,complex<double> z2,
    complex<double> w1,complex<double> w2)
// Used by brownian_outer_loop
// see if the segments [z1,z2] and [w1,w2] intersect. 
// If they intersect return the intersection, else return 1.e10
// we return +1 or -1 depending on a right/left handed crossing. 
// If they do not intersect we return 0
// If they intersect in an endpoint we say they don't intersect
// NB: this means there is a cutoff in the routine.
{
    double x1,x2,y1,y2,t,x;

    // considering the quadrilateral with vertices z1,z2,w1,w2 we
    // see that if they intersect then |z1-w1|+|z2-w2|<=|z1-z2|+|w1-w2|
    if (fabs(z1.real()-w1.real())+fabs(z2.real()-w2.real())
        >fabs(z1.real()-z2.real())+fabs(w1.real()-w2.real())) return(1.e10);

    // shift z1 to 0
    z2-=z1;
    w1-=z1;
    w2-=z1;
    // rotate and scale so z2=1
    w1/=z2;
    w2/=z2;

    x1=w1.real();
    y1=w1.imag();
    x2=w2.real();
    y2=w2.imag();
    if (y1*y2>-1.e-10) return(1.e10);
    
    // intersection when (1-t)*w1+t*w2 in [0,1] 
    // In particular, (1-t)*y1+t*y2=0, i.e., t=y1/(y1-y2)
    // then we need (1-t)*x1+t*x2 in [0,1]
    t=y1/(y1-y2);
    x=(1-t)*x1+t*x2;
    if (x<1.e-10 || x>1.-1.e-10) return(1.e10);

    w1*=z2;
    w2*=z2;
    w1+=z1;
    w2+=z1;

    complex<double> zz1;
    zz1=(1-t)*w1+t*w2;
    return(zz1);         

}

void complex_walk::brownian_loop(long ns)
{
    // Brownian loop is obtained by taking independent Brownian bridges
    // for x and y components.
    // W_t=B_t - t B_1, 0 <=t<=1, B_t is standard BM

    long istep;
    double *x,*y,stdev;
    complex<double> ztotal,z,dz,z2,zz,znear,zstart,zprev;

    ///////////////////  Generate the loop ////////////////////////

    stdev=sqrt(1./double(ns));
    x = new double[ns];
    y = new double[ns];

    for (istep=0;istep<ns;istep++) {
        x[istep]=normal_rv()*stdev;
        y[istep]=normal_rv()*stdev;
    }

    steps[0]=0.;
    ztotal=0.;
    for (istep=1;istep<=ns;istep++) {
        ztotal+=x[istep-1]+y[istep-1]*I;
        steps[istep]=ztotal;
    }

    for (istep=1;istep<=ns;istep++) {
        steps[istep]-=ztotal*double(istep)/double(ns);
    }

    nsteps=ns;
    delete [] x;
    delete [] y;
}

void complex_walk::brownian_outer_loop(long ns)
{
    // Brownian loop is obtained by taking independent Brownian bridges
    // for x and y components.
    // W_t=B_t - t B_1, 0 <=t<=1, B_t is standard BM

double etime= -mytime(); // delete me

    long istep,jstep,inext,inear=0,jjump=0;
    double *x,*y,stdev;
    complex<double> *loop,ztotal,z,dz,z2,zz,znear,zstart,zprev;
    double intersect_cut=1.e-8;

    ///////////////////  Generate the loop ////////////////////////

    loop = new complex<double>[ns+1]; // the Brownian loop

    stdev=sqrt(1./double(ns));
    x = new double[ns];
    y = new double[ns];

    for (istep=0;istep<ns;istep++) {
        x[istep]=normal_rv()*stdev;
        y[istep]=normal_rv()*stdev;
    }

    loop[0]=0.;
    ztotal=0.;
    for (istep=1;istep<=ns;istep++) {
        ztotal+=x[istep-1]+y[istep-1]*I;
        loop[istep]=ztotal;
    }

    double max_dist=0.,dist;
    long istep_max=0;
    for (istep=1;istep<=ns;istep++) {
        loop[istep]-=ztotal*double(istep)/double(ns);
        dist=abs(loop[istep]-loop[istep-1]);
        if (dist>max_dist) {max_dist=dist; istep_max=istep;}
    }

    double max_dist2=0.;
    for (istep=1;istep<=ns;istep++) if (istep!=istep_max) {
        dist=abs(loop[istep]-loop[istep-1]);
        if (dist>max_dist2) max_dist2=dist; 
    }

printf("BM loop done %lf, max=%lf max2=%lf \n",mytime()+etime,
max_dist,max_dist2);

    FILE *fptr;
    fptr=fopen("temp/bloop.plt","w");   
    for (istep=1;istep<=ns;istep++)
       fprintf(fptr,"%lf %lf \n",loop[istep].real(),loop[istep].imag());
    fclose(fptr);

printf("BM loop done %lf \n",mytime()+etime);


    /////////////////////  Find outer boundary ////////////////////////

    // find rightmost site on the loop 
    double x1,y1,x2,y2;
    double xmax=-1.e10;
    long imax=0,incr_flag;
    for (istep=0;istep<ns;istep++) if (loop[istep].real()>xmax) {
        xmax=loop[istep].real();
        imax=istep;
    }
    zstart=loop[imax];
    zprev=zstart;

    // we will traverse outer boundary in counterclockwise direction
    // incr_flag=1 (-1) means the index increases (decreases)
    // (x1,y1)=loop[imax]-loop[imax-1]
    // (x2,y2)=loop[imax]-loop[imax+1]
    // (x1,y1,0) x (x2,y2,0)=x1*y2-x2*y1
    x1=loop[imax].real()-loop[(imax+ns-1)%ns].real();
    y1=loop[imax].imag()-loop[(imax+ns-1)%ns].imag();
    x2=loop[imax].real()-loop[(imax+1)%ns].real();
    y2=loop[imax].imag()-loop[(imax+1)%ns].imag();
    incr_flag = (x1*y2-x2*y1<0);
    if (incr_flag) inext=(imax+1)%ns;
    else inext=(imax+ns-1)%ns;

    zprev=zstart;
    nsteps=0;

    steps[0]=zstart;
    while (abs(zprev-zstart)>1.e-10 || nsteps==0) {
        // Integers inext,incr_flag should be defined here.
        // If incr_flag, we are moving along loop[inext-1] -> loop[inext].
        // Else we are moving along loop[inext+1] -> loop[inext].
        z2=loop[inext];

	// find intersection (if any) on this edge closest to zprev
	// and in the direction we are going
        znear=1.e10;
long count=0; // delete me
        for (istep=1;istep<=ns;istep+=jjump) {
count++; // delete me
            zz=segment_intersect(zprev,z2,loop[istep],loop[istep-1]);

            if (zz.real()<1.e9) {
                // KLUDGE  KLUDGE  KLUDGE  KLUDGE  KLUDGE  KLUDGE  KLUDGE 
                if (abs(zz-zprev)<intersect_cut || abs(zz-z2)<intersect_cut ||
                    abs(zz-loop[istep])<intersect_cut ||
                    abs(zz-loop[istep-1])<intersect_cut) zz=1.e10; 
	        if (abs(zz-zprev)<abs(znear-zprev)) {
                    znear=zz;
                    inear=istep;
                }
            }
            // speed up search for intersections by taking advantage of 
            // fact that steps in the brownian loop are in order and 
            // go no further than max_dist.
            // So if dist=|loop[istep]-zprev|, the next possible intersection
            // is from [loop[istep+jjump],loop[istep+jjump-1]]
            // where jjump=dist/max_dist-1
            dist=abs(loop[istep]-zprev);
            jjump=long(floor(dist/max_dist))-1;
            if (jjump<1) jjump=1;
        }

        if (znear.real()<1.e9) { 
            // Next point is the closest intersection
            // inear-1, inear will be the next edge
            x1=z2.real()-zprev.real();
            y1=z2.imag()-zprev.imag();
            x2=loop[inear].real()-loop[(inear+ns-1)%ns].real();
            y2=loop[inear].imag()-loop[(inear+ns-1)%ns].imag();
            incr_flag = (x1*y2-x2*y1<0);
            if (incr_flag) inext=inear;
            else inext=(inear+ns-1)%ns;
            zprev=znear;
	}
        else { 
            // Next point is endpoint of the edge
            // Next edge is next edge in the loop
            if (incr_flag) inext=(inext+1)%ns;
            else inext=(inext+ns-1)%ns;
            zprev=z2;
	}

        nsteps++;
        if (nsteps>max_nsteps) {
            printf("mem exceed in brownian loop %ld\n",max_nsteps);
            nsteps=0; 
            delete [] x;
            delete [] y;
            delete [] loop;
            return;
	}

        steps[nsteps]=zprev;

	if (nsteps%10000==0) printf("boundary has %ld points so far\n",nsteps);
    }

printf("outer loop done %lf \n",mytime()+etime);

    /////////////////////  Eliminate cut points  ////////////////////////

    // Eliminate cut points by trimming small dangling loops 
    // NB : first site in the loop may belong to such a loop. 
    // If so, must find a starting point that does not.
    long ishift=0,jshift=0,icut=0,jcut=0,loop_flag;
    do {
        loop_flag=0;
        for (istep=ishift+1;istep<=ishift+nsteps/10;istep++) 
        for (jstep=jshift-1;jstep>=jshift-nsteps/10;jstep--) 
        if (abs(steps[istep]-steps[(jstep+nsteps)%nsteps])<1.e-10) {
            icut=istep;
            jcut=jstep;
            loop_flag=1;
        }
        if (loop_flag) {
            ishift=icut;
            jshift=jcut;
	}
    } while (loop_flag);
    if (jshift<0) nsteps=(jshift+nsteps)%nsteps;



    // remove the loops by shifting steps[] 

    steps[0]=steps[ishift];
    for (istep=ishift+1;istep<=nsteps;) {
        loop_flag=0;
        for (jstep=1;jstep<=nsteps/10 && (istep+jstep)<=nsteps;jstep++) {
            if (abs(steps[istep]-steps[istep+jstep])<1.e-10) {
                loop_flag=1;
                jjump=jstep;
            }
        }
        if (loop_flag) {
            istep+=jjump;
            ishift+=jjump;
        }
        else {
            steps[istep-ishift]=steps[istep];
            istep++;
	}
    }
    nsteps-=ishift;

    delete [] x;
    delete [] y;
    delete [] loop;

printf("cut points done %lf \n",mytime()+etime);

}

void complex_walk::rw_loop(long ns)
{
    // Run an ordinary rw for ns steps. Find last time it was at origin.
    // Discard walk after this point. So we have a loop.
    // Is this legit?
    // Find outer boundary of loop.

    long istep,jstep,inext,inear=0,jjump=0;
    double *x,*y,stdev;
    complex<double> *loop,ztotal,z,dz,z2,zz,znear,zstart,zprev;
    double intersect_cut=1.e-8;

    ///////////////////  Generate the loop ////////////////////////

    loop = new complex<double>[ns+1]; // the Brownian loop

    stdev=sqrt(1./double(ns));
    x = new double[ns];
    y = new double[ns];

    complex<double> ii;
    ii=I;

    loop[0]=0.;
    for (istep=1;istep<=ns;istep++) {
        do {
            loop[istep]=loop[istep-1];
            double xrand=RNG();
            if (xrand<0.25) loop[istep]+=stdev;
            else if (xrand<0.5) loop[istep]-=stdev;
            else if (xrand<0.75) loop[istep]+=stdev*ii;
            else loop[istep]-=stdev*ii;
        } while (istep>1 && 
            abs(loop[istep]-loop[istep-2])<1.e-8);
    }

    // find last visit to origin
    for (istep=ns;istep>=0;istep--) if (abs(loop[istep])<1.e-10) 
      {ns=istep;
       break;
      }

    printf("last visit to origin is at %ld \n",ns);

    double max_dist=0.,dist;
    for (istep=1;istep<=ns;istep++) {
        loop[istep]-=ztotal*double(istep)/double(ns);
        dist=abs(loop[istep]-loop[istep-1]);
        if (dist>max_dist) max_dist=dist;        
    }

    FILE *fptr;
    fptr=fopen("temp/bloop.plt","w");   
    for (istep=1;istep<=ns;istep++)
       fprintf(fptr,"%lf %lf \n",loop[istep].real(),loop[istep].imag());
    fclose(fptr);


    /////////////////////  Find outer boundary ////////////////////////

    // find rightmost site on the loop 
    double x1,y1,x2,y2;
    double xmax=-1.e10;
    long imax=0,incr_flag;
    for (istep=0;istep<ns;istep++) if (loop[istep].real()>xmax) {
        xmax=loop[istep].real();
        imax=istep;
    }
    zstart=loop[imax];
    zprev=zstart;

    // we will traverse outer boundary in counterclockwise direction
    // incr_flag=1 (-1) means the index increases (decreases)
    // (x1,y1)=loop[imax]-loop[imax-1]
    // (x2,y2)=loop[imax]-loop[imax+1]
    // (x1,y1,0) x (x2,y2,0)=x1*y2-x2*y1
    x1=loop[imax].real()-loop[(imax+ns-1)%ns].real();
    y1=loop[imax].imag()-loop[(imax+ns-1)%ns].imag();
    x2=loop[imax].real()-loop[(imax+1)%ns].real();
    y2=loop[imax].imag()-loop[(imax+1)%ns].imag();
    incr_flag = (x1*y2-x2*y1<0);
    if (incr_flag) inext=(imax+1)%ns;
    else inext=(imax+ns-1)%ns;

    zprev=zstart;
    nsteps=0;

    steps[0]=zstart;
    while (abs(zprev-zstart)>1.e-10 || nsteps==0) {
        // Integers inext,incr_flag should be defined here.
        // If incr_flag, we are moving along loop[inext-1] -> loop[inext].
        // Else we are moving along loop[inext+1] -> loop[inext].
        z2=loop[inext];

	// find intersection (if any) on this edge closest to zprev
	// and in the direction we are going
        znear=1.e10;
        for (istep=1;istep<=ns;istep+=jjump) {
            zz=segment_intersect(zprev,z2,loop[istep],loop[istep-1]);

            if (zz.real()<1.e9) {
                // KLUDGE  KLUDGE  KLUDGE  KLUDGE  KLUDGE  KLUDGE  KLUDGE 
                if (abs(zz-zprev)<intersect_cut || abs(zz-z2)<intersect_cut ||
                    abs(zz-loop[istep])<intersect_cut ||
                    abs(zz-loop[istep-1])<intersect_cut) zz=1.e10; 
	        if (abs(zz-zprev)<abs(znear-zprev)) {
                    znear=zz;
                    inear=istep;
                }
            }
            // speed up search for intersections by taking advantage of 
            // fact that steps in the brownian loop are in order and 
            // go no further than max_dist.
            // So if dist=|loop[istep]-zprev|, the next possible intersection
            // is from [loop[istep+jjump],loop[istep+jjump-1]]
            // where jjump=dist/max_dist-1
            dist=abs(loop[istep]-zprev);
            jjump=long(floor(dist/max_dist))-1;
            if (jjump<1) jjump=1;
        }

        if (znear.real()<1.e9) { 
            // Next point is the closest intersection
            // inear-1, inear will be the next edge
            x1=z2.real()-zprev.real();
            y1=z2.imag()-zprev.imag();
            x2=loop[inear].real()-loop[(inear+ns-1)%ns].real();
            y2=loop[inear].imag()-loop[(inear+ns-1)%ns].imag();
            incr_flag = (x1*y2-x2*y1<0);
            if (incr_flag) inext=inear;
            else inext=(inear+ns-1)%ns;
            zprev=znear;
	}
        else { 
            // Next point is endpoint of the edge
            // Next edge is next edge in the loop
            if (incr_flag) inext=(inext+1)%ns;
            else inext=(inext+ns-1)%ns;
            zprev=z2;
	}

        nsteps++;
        if (nsteps>max_nsteps) {
            printf("mem exceed in brownian loop \n");
            exit(0);
	}

        steps[nsteps]=zprev;

	if (nsteps%10000==0) printf("boundary has %ld points so far\n",nsteps);
    }

    /////////////////////  Eliminate cut points  ////////////////////////

    // Eliminate cut points by trimming small dangling loops 
    // NB : first site in the loop may belong to such a loop. 
    // If so, must find a starting point that does not.
    long ishift=0,jshift=0,icut=0,jcut=0,loop_flag;
    do {
        loop_flag=0;
        for (istep=ishift+1;istep<=ishift+nsteps/10;istep++) 
        for (jstep=jshift-1;jstep>=jshift-nsteps/10;jstep--) 
        if (abs(steps[istep]-steps[(jstep+nsteps)%nsteps])<1.e-10) {
            icut=istep;
            jcut=jstep;
            loop_flag=1;
        }
        if (loop_flag) {
            ishift=icut;
            jshift=jcut;
	}
    } while (loop_flag);
    if (jshift<0) nsteps=(jshift+nsteps)%nsteps;



    // remove the loops by shifting steps[] 

    steps[0]=steps[ishift];
    for (istep=ishift+1;istep<=nsteps;) {
        loop_flag=0;
        for (jstep=1;jstep<=nsteps/10 && (istep+jstep)<=nsteps;jstep++) {
            if (abs(steps[istep]-steps[istep+jstep])<1.e-10) {
                loop_flag=1;
                jjump=jstep;
            }
        }
        if (loop_flag) {
            istep+=jjump;
            ishift+=jjump;
        }
        else {
            steps[istep-ishift]=steps[istep];
            istep++;
	}
    }
    nsteps-=ishift;

    delete [] x;
    delete [] y;
    delete [] loop;
}

void complex_walk::circle(long ns)
// justs creates a circle centered at origin with radius 1
{
    double theta;
    complex<double> ii;
    ii= -1.;
    ii=sqrt(ii);
    nsteps=ns;
    for (long i=0;i<ns;i++) {
        theta=2*M_PI*double(i)/double(ns);
        steps[i]=cos(theta)+sin(theta)*ii;
    }
    steps[ns]=steps[0];    
}


long complex_walk::interior(complex<double> z)
// look at intersections of horizontal line through z with loop
{
    long nleft=0,i;
    double x,y,x1,y1,x2,y2,t;
    y=z.imag();
    x1=steps[nsteps].real();
    y1=steps[nsteps].imag();

    for (i=0;i<=nsteps;i++) {
        x2=steps[i].real();
        y2=steps[i].imag();
        if ((y1<y && y<y2) || (y1>y && y>y2)) {
            t=(y-y2)/(y1-y2);
            x=t*x1+(1.-t)*x2;
            if (x<z.real()) nleft++;
        }
        x1=x2;
        y1=y2;
    }
    if (nleft%2==1) return(1L);
    else return(0L);


	////////////////////////////////////////////
	////////////////////////////////////////////
	////////////////////////////////////////////
// tests if z is in the interior of the walk. Treats walk as a loop.
// we compute winding number of loop around the point

    double angle,dangle;
    long di=1;
    x1=(steps[nsteps]-z).real();
    y1=(steps[nsteps]-z).imag();
    x2=(steps[0]-z).real();
    y2=(steps[0]-z).imag();
    angle=asin((x2*y1-x1*y2)/sqrt((x1*x1+y1*y1)*(x2*x2+y2*y2)));

    for (long i=0;i+di<=nsteps;i+=di) {
        x1=(steps[i]-z).real();
        y1=(steps[i]-z).imag();
        x2=(steps[i+di]-z).real();
        y2=(steps[i+di]-z).imag();
        dangle=asin((x2*y1-x1*y2)/sqrt((x1*x1+y1*y1)*(x2*x2+y2*y2)));
        if (!isnan(dangle)) angle+=dangle;
    }
    angle=fabs(angle)/(2*M_PI);
    if (fabs(angle-1)<1.e-2) return(1L);    
    else return(0L); 
}

int area_cmp(const void *x1,const void *x2)
{
    if ((*(double *)x1) < (*(double *)x2)) return(-1);
    if ((*(double *)x1) > (*(double *)x2)) return(1);
    return(0);
}

double complex_walk::area()
// computes area enclosed by walk that is a loop 
// Green's thm says this area is (up to sign) integral of y dx and is also 
// integral of x dy. We compute both to check they agree.
{
    long is;
    double x1,x2,y1,y2,area1,area2;
 
    // area1 is integral of x dy;  area2 is integral of y dx
    area1=0.;
    area2=0.;
    x1=steps[0].real();
    y1=steps[0].imag();
    for (is=1;is<=nsteps;is++) {
        x2=steps[is].real();
        y2=steps[is].imag();
        area1+=(x1+x2)*(y2-y1)/2.;                            
        area2+=(y1+y2)*(x2-x1)/2.;                            
        x1=x2;
        y1=y2;
    } 

    // last step; if steps[nsteps]==steps[0] this simply adds 0 to areas
    x2=steps[0].real();
    y2=steps[0].imag();
    area1+=(x1+x2)*(y2-y1)/2.;                            
    area2+=(y1+y2)*(x2-x1)/2.;                            

    area1=fabs(area1);
    area2=fabs(area2);

    if (fabs(area1-area2)/area1>1.e-10) {
        printf(">>>>> ERROR: error in area %lf %lf\n",area1,area2);
    }
    return((area1+area2)/2.);
}

double complex_walk::diameter()
// computes diameter of walk that is a loop
{
    double delta,dist,max_dist=0.;
    long np,*p,max_np=100000,ip,jp,i,j,k;
    p = new long[max_np+1];

    // use fact=0.04;

for (double fact=0.04;fact<0.05;fact*=2) {

    delta=fact*max_radius();
    p[0]=0;
    np=1;
    for (k=1;k<=nsteps;k++) {
        if (abs(steps[k]-steps[p[np-1]])>=delta) {
            if (np>=max_np-1) {
                printf("max_np too small in diameter \n"); 
                exit(1L);
            }
            p[np++]=k;
        }
    }
    p[np]=nsteps+1;
    double cdiam=0.;
    for (ip=0;ip<np;ip++) for (jp=ip+1;jp<np;jp++) {
        dist=abs(steps[p[ip]]-steps[p[jp]]);
        if (dist>cdiam) cdiam=dist;
    }


    for (ip=0;ip<np;ip++) for (jp=ip+1;jp<np;jp++) {
        dist=abs(steps[p[ip]]-steps[p[jp]]);
        if (dist>cdiam-2*delta) 
        for (i=p[ip];i<p[ip+1];i++) for (j=p[jp];j<p[jp+1];j++) {
            dist=abs(steps[i]-steps[j]);
            if (dist>max_dist) max_dist=dist;
        }
    }
}
    ////////////////////////////////////////////////////
  
    delete [] p;
    return(max_dist);

    max_dist=0.;
    long di=10,i1,i2;
    for (i1=0;i1+di<=nsteps;i1+=di) for (i2=i1+di;i2+di<=nsteps;i2+=di) {
        dist=abs(steps[i1]-steps[i2]);
        if (dist>max_dist) max_dist=dist;
    }
    printf("di=%ld max_dist=%lf \n",di,max_dist);

    di=5;
    for (i1=0;i1+di<=nsteps;i1+=di) for (i2=i1+di;i2+di<=nsteps;i2+=di) {
        dist=abs(steps[i1]-steps[i2]);
        if (dist>max_dist) max_dist=dist;
    }
    printf("di=%ld max_dist=%lf \n",di,max_dist);

    di=2;
    for (i1=0;i1+di<=nsteps;i1+=di) for (i2=i1+di;i2+di<=nsteps;i2+=di) {
        dist=abs(steps[i1]-steps[i2]);
        if (dist>max_dist) max_dist=dist;
    }
    printf("di=%ld max_dist=%lf \n",di,max_dist);

    di=1;
    for (i1=0;i1+di<=nsteps;i1+=di) for (i2=i1+di;i2+di<=nsteps;i2+=di) {
        dist=abs(steps[i1]-steps[i2]);
        if (dist>max_dist) max_dist=dist;
    }
    printf("di=%ld max_dist=%lf \n",di,max_dist);

    return(max_dist); 
}

complex<double> complex_walk::cm_interior()
// finds center of mass of interior of loop
{
    #define MAX_NX 1000
    double min_x,max_x,min_y,max_y,x,y,x1,y1,x2,y2,dy,t,xlist[MAX_NX],
        xcm,ycm,area;
    long ngrid,i,ix,iy;
    complex<double> z,ii;
    ii=-1;
    ii=sqrt(ii);

    // find rectangle that contains loop
    min_x=1.e10; max_x=-1.e10; min_y=1.e10; max_y=-1.e10;
    for (i=0;i<=nsteps;i++) {
        x=steps[i].real();
        y=steps[i].imag();
        if (x>max_x) max_x=x;
        if (x<min_x) min_x=x;
        if (y>max_y) max_y=y;
        if (y<min_y) min_y=y;
    }

    // for array of horizontal lines find their intersections with loop
    ngrid=2000;
    dy=(max_y-min_y)/double(2*ngrid);
    xcm=0.;
    ycm=0.;
    area=0.;
    for (iy=1;iy<2*ngrid;iy+=2) {
        y=min_y + iy*dy;
        x1=steps[nsteps].real();
        y1=steps[nsteps].imag();

        ix=0;
        for (i=0;i<=nsteps;i++) {
            x2=steps[i].real();
            y2=steps[i].imag();
            if ((y1<y && y<y2) || (y1>y && y>y2)) {
                t=(y-y2)/(y1-y2);
                x=t*x1+(1.-t)*x2;
                if (ix>=MAX_NX-1) {
                    printf("MAX_NX too small in complex_walk::area()\n");
                    exit(0);
                }
                xlist[ix++]=x;
            }
            x1=x2;
            y1=y2;
        }
        qsort(xlist,ix,sizeof(double),area_cmp);
        if (ix%2==0) for (i=0;i<ix;i+=2) {
            area+=xlist[i+1]-xlist[i];
            ycm+=(xlist[i+1]-xlist[i])*y;
            xcm+=(xlist[i+1]-xlist[i])*(xlist[i+1]+xlist[i])/2.;
        }
        else printf("ODD number of intersections\n");

    }
    xcm/=area;
    ycm/=area;
//    area*=2*dy;

FILE *fptr;
fptr=fopen("temp/loop","w");
plot(fptr);
fclose(fptr);

fptr=fopen("temp/cm","w");
fprintf(fptr,"%lf %lf\n",xcm,ycm);
fclose(fptr);


    z=xcm+ii*ycm;
    return(z);


}

double complex_walk::radius_cm()
// max dist from center of mass to loop 
{
    double max=0.,temp;
    complex<double> cm;
    cm=cm_interior();
    for (long i=0;i<=nsteps;i++) {
        temp=abs(steps[i]-cm);
        if (temp>max) max=temp; 
    }
    return(max);
}

double complex_walk::approach_cm()
// min dist from center of mass to loop 
{
    double min=1.e10,temp;
    complex<double> cm;
    cm=cm_interior();
    for (long i=0;i<=nsteps;i++) {
        temp=abs(steps[i]-cm);
        if (temp<min) min=temp; 
    }
    return(min);
}

double complex_walk::radius_smallest_disc()
// radius of smallest disc enclosing a walk that is a loop
{
    printf("complex_walk::radius_smallest_disc() not implemented \n");
    exit(0); 
}

///////////////////////////////////////////////////////////////
//  checkk checkk checkk checkk checkk checkk checkk checkk  //
///////////////////////////////////////////////////////////////

