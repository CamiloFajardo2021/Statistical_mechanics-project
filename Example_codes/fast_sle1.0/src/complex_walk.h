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

// file : src/complex_walk.h 
// complex_walk is a walk with complex points 

#ifndef __complex_walk_hh
#define __complex_walk_hh

#include "local_constants.h"

class laurent;

class complex_walk {
// stores a single walk consisting of complex numbers.

// nsteps is the number of steps in the SLE. 
// In between these steps there are dsteps time intervals. 
// So the number of time intervals is dsteps*nsteps
// dt[k], k=1,2,...,dsteps*nsteps are the time intervals.
// dx[k] is change in driving function (root(kappa)*brownian motion 
// (or its approx) over that time interval.
// t[i] is the time for steps[i]. So t[i] is sum of dt[k] for 
// k=1,2,...,dsteps*k

private: 
public:
    // steps and nsteps should be private eventually
    complex<double>* steps;
    long nsteps; // number of steps in the walk  
    long max_nsteps; // max number of steps in the walk for memory allocation
    long dsteps;
    long bm_nsteps;
    long map_choice; // type of conformal map used to discretize SLE
    long sle_type; // ==1 for chordal SLE, ==2 for radial SLE, ==3 for dipolar
    long sle_domain; // domain for chordal SLE
    // 1=half plane, 2=circle, 3=strip
    long nblock;
    long blength;
    long nterms;
    double laurent_factor;
    double kappa;
    double drift;

    // SLE(kappa,rho)
    long nforce_pt; // number of force points
    double force_pt[MAX_NFORCE_PT]; // initial locations of force points
    double rho[MAX_NFORCE_PT]; // as in SLE(kappa,rho)

    long num_unzip; // used to record number of steps that were unzipped
    double *dt;
    double *dx;
    double *t;
    // parameters for conformal maps. We compute these and store them in 
    // arrays. It would be more elegant to just send dt and dx to 
    // tilted_slit_cmap() and compute alpha,xl,xr inside that function.
    // But this would probably be a lot slower. 
    double *alpha,*xl,*xr,*aa,*bb;
    complex<double> *ww;
    // sewing map: sewing_plus[i] and sewing_minus[i] are sewn together
    double *sewing_minus,*sewing_plus;

    // used for blocks
    laurent *fblock;
    double *radius;

    complex_walk();
    void allocate(long n, long ds, long mc, long bl, long ch); 
    void deallocate();
    
    void initialize(); 
    long get_nsteps();
    void set_nsteps(long n);
    void set_dsteps(long n);
    void set_max_nsteps(long n);
    complex<double> get_step(long i);
    void assign_step(long i,complex<double> z);
    void print(FILE *fptr);
    void scan(FILE *fptr);
    void plot(FILE *fptr);
    void plot_skip(FILE *fptr,long dn);
    void rescale(double c); 
    void invert();
    void read_bm(FILE *fptr);
    void write_bm(FILE *fptr);
    void sample_bm(double max_dt, double *next_dt);
    void cauchy_deterministic_dt();
    void bm_drift_deterministic_dt();
    void bm_deterministic_dt(long dt_choice,long dx_choice);
    void compute_block_ch_zip(long ib);
    void compute_block_ch_unzip(long ib);
    void compute_block_rad_zip(long ib);
    void compute_block_half_rad_zip(long ib);
    void compute_block_di_zip(long ib);
    void compute_parms_ch_zip(long istep,double ddt,double ddx);
    void compute_parms_rad_zip(long istep,double ddt,double ddx);
    void compute_parms_ch_unzip(long ipt,complex<double> z);
    void compute_parms_di_zip(long istep,double ddt,double ddx);
    void ch_zip(long ns,long ds,long mc,long bl,long nt,double lf);
    void adapt_chordal_sle_laurent(long sle_domain,complex<double> sle_z1, 
      complex<double> sle_z2,long nsteps, long ds,
      long map_choice, long bm_approx, double kappa, double drift,
      long stop_choice,
      double stop_max, double nu, double min_dt, double max_dt, double max_dx, 
      long blength,long nterms,double laurent_factor);
    void adapt_radial_sle_laurent(long sle_domain, complex<double> sle_z1, 
      complex<double> sle_z2, long nsteps, long ds,
      long map_choice, long bm_approx, double kappa, long stop_choice,
      double stop_max, double nu, double min_dt, double max_dt, double max_dx, 
      long blength,long nterms,double laurent_factor);
    void adapt_half_radial_sle_laurent(long sle_domain, complex<double> sle_z1, 
      complex<double> sle_z2, long nsteps, long ds,
      long map_choice, long bm_approx, double kappa, long stop_choice,
      double stop_max, double nu, double min_dt, double max_dt, double max_dx, 
      long blength,long nterms,double laurent_factor);
    void ch_unzip(long npt,long blength,long nterms,double laurent_factor,
      long map_choice,double max_t); 
    void adapt_dipolar_sle_laurent(long sle_domain, complex<double> sle_z1, 
      complex<double> sle_z2, long nsteps, long ds,
      long map_choice, long bm_approx, double kappa, long stop_choice,
      double stop_max, double nu, double min_dt, double max_dt, double max_dx, 
      long blength,long nterms,double laurent_factor);
    complex<double> conformal_map_ch_zip(long ns, long imag_flag, 
       complex<double> z);
    complex<double> conformal_map_ch_unzip(long ipt, double *rad, 
       complex<double> z);

    // qqqqqqqqqqqqqqqq
    void conformal_map_ch_sewing(long ipt, long npts, double *xm,double *xp);
    //qqqqqqqqqqqqqqqq
    complex<double> conformal_map_rad_zip(long ns, complex<double> z);
    complex<double> conformal_map_half_rad_zip(long ns, complex<double> z);
    complex<double> conformal_map_di_zip(long ns, long imag_flag, 
       complex<double> z);

    void ch_sewing(long npts);

    double max_radius();
    void brownian_outer_loop(long nsteps);
    void brownian_loop(long nsteps);
    void rw_loop(long nsteps);
    void circle(long ns);
    long interior(complex<double> z);
    double area();
    double diameter();
    double radius_smallest_disc();
    complex<double> cm_interior();
    double radius_cm();
    double approach_cm();

    void generate_radial_rw_continuum(long explorer_type,long domain_type, 
        double *domain_parms);

}; 

double max_distance(complex_walk *cw1,complex_walk *cw2);
double average_distance(complex_walk *cw1,complex_walk *cw2);

#endif

///////////////////////////////////////////////////////////////
//  checkk checkk checkk checkk checkk checkk checkk checkk  //
///////////////////////////////////////////////////////////////

