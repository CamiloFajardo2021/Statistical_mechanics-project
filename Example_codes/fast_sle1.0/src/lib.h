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


// file name : src/lib.h
//

#ifndef __lib_hh
#define __lib_hh

double mytime();
// returns system clock in secs 

double cauchy_rv(); // returns "standard" cauchy RV 

double normal_rv(); // returns a normal random variable with mean 0, variance 1

// least squares fit
void least_squares(long nx,double* y,double* x1,double *sd,
    double* beta,double* beta_err,double *rss); 

#endif

///////////////////////////////////////////////////////////////
//  checkk checkk checkk checkk checkk checkk checkk checkk  //
///////////////////////////////////////////////////////////////
