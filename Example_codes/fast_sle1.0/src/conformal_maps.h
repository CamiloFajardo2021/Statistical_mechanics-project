// file: src/conformal_maps.h

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

#ifndef __conformal_maps_hh
#define __conformal_maps_hh

complex<double> vertical_slit_ch_zip(double dt,double dx,complex<double> z);
complex<double> tilted_slit_ch_zip(
    double alpha, double xl, double xr,complex<double> z);
complex<double> arc_ch_zip(double a, double b, complex<double> z);

complex<double> vertical_slit_ch_unzip(double dt,double dx,
    complex<double> z);
complex<double> tilted_slit_ch_unzip(
    double alpha, double xl, double xr,complex<double> z,complex<double> iz);
complex<double> arc_ch_unzip(double a, double b, complex<double> z);

// conformal maps on real axis
double tilted_slit_ch_zip_real(double alpha, double xl, double xr,double x);
double tilted_slit_ch_unzip_real(double alpha, double xl, double xr,double x);
void tilted_slit_ch_unzip_origin(double alpha, double xl, double xr,
			       double *xm, double *xp);
void arc_ch_unzip_origin(double a, double b, double *xm, double *xp);
double arc_ch_unzip_real(double a, double b, double xm);

complex<double> half_to_strip(complex<double> z,complex<double> z1);
complex<double> disc_to_half(complex<double> z,complex<double> z1);

complex<double> vertical_slit_rad_zip(
    double dt,double dx,complex<double> z);
complex<double> tilted_slit_rad_zip(double dx,
    double alpha, double xl, double xr, complex<double> ww, complex<double> z);

complex<double> vertical_slit_half_rad_zip(
    double dt,double dx,complex<double> z);
complex<double> tilted_slit_half_rad_zip(double dx,
    double alpha, double xl, double xr, complex<double> ww, complex<double> z);

complex<double> vertical_slit_di_zip(double dt,double dx,complex<double> z);

complex<double> assign(double x, double y);
complex<double> sc_triangle_integrand(complex<double> z);
complex<double> sc_square_integrand(complex<double> z);
complex<double> sc_square_rotate_90_hplane(complex<double> z);
complex<double> sc_square_rotate_45_hplane(complex<double> z);
complex<double> sc_triangle_rotate_120_hplane(complex<double> z);
complex<double> sc_triangle_rotate_60_hplane(complex<double> z);

complex<double> sc_triangle_integrate_ps(complex<double> z);
complex<double> sc_triangle_integrate(complex<double> z);
complex<double> sc_triangle_map(complex<double> z);
complex<double> sc_triangle_map_inv(complex<double> z);

complex<double> sc_square_integrate_ps(complex<double> z);
complex<double> sc_square_integrate(complex<double> z);
complex<double> sc_square_map(complex<double> z);
complex<double> sc_square_map_inv(complex<double> z);

complex<double> scd_square_integrate_ps(complex<double> z);
complex<double> scd_square_integrate(complex<double> z);
complex<double> scd_square_map(complex<double> z);
complex<double> scd_square_map_inv(complex<double> z);

double sc_triangle_theta(double frac, double t);
double sc_square_theta(double frac, double t);

double sc_triangle_theta_test(double frac, double t);
complex<double> sc_triangle_integrate_ps_test(complex<double> z);
complex<double> sc_triangle_map_test(complex<double> z);
complex<double> sc_triangle_map_inv_test(complex<double> z);

double sc_square_theta_test(double frac, double t);
complex<double> sc_square_integrate_ps_test(complex<double> z);
complex<double> sc_square_map_test(complex<double> z);
complex<double> sc_square_map_inv_test(complex<double> z);

#endif

///////////////////////////////////////////////////////////////
//  checkk checkk checkk checkk checkk checkk checkk checkk  //
///////////////////////////////////////////////////////////////
