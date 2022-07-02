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


// file name : src/laurent.h
//

#ifndef __laurent_hh
#define __laurent_hh

#define MAX_NTERMS 30

class laurent {
private: 
public:
    long nterms;  // should be private eventually ?
    complex<double> coef[MAX_NTERMS+1];  // should be private eventually ?

    long get_nterms();
    void fprint(FILE *fptr);
    void assign(long n,complex<double> *c);

    inline laurent operator+=(laurent f1) 
    {
        if (nterms!=f1.nterms) {
            printf("ERROR in operator+\n");
            exit(1);  
        }
        for (long i=0;i<=nterms;i++) coef[i]+=f1.coef[i];
        return *this;
    } 

    inline laurent operator+(laurent f) 
    {
        *this+=f;
        return *this;
    } 

    inline laurent operator-=(laurent f1) 
    {
        if (nterms!=f1.nterms) {
            printf("ERROR in laurent::operator-\n");
            exit(1);  
        }
        for (long i=0;i<=nterms;i++) coef[i]-=f1.coef[i];
        return *this;
    } 

    inline laurent operator-(laurent f) 
    {
        *this-=f;
        return *this;
    } 

    inline laurent operator*=(laurent f1) 
    {
        laurent f;
        complex<double> temp;
        f=*this;
        for (long i1=0;i1<=nterms;i1++) coef[i1]=0.;
        for (long i1=0;i1<=nterms;i1++) for (long i2=0;i2<=nterms-i1;i2++) {
            coef[i1+i2]+=f1.coef[i1]*f.coef[i2];
        }
        return *this;
    } 

    inline laurent operator*(laurent f) 
    {
        (*this)*=f;
        return *this;
    } 

    inline laurent operator*=(complex<double> c) 
    {
        for (long i1=0;i1<=nterms;i1++) coef[i1]*=c;
        return *this;
    } 

    laurent operator=(const laurent f1) 
    {
        nterms=f1.nterms;
        for (long i=0;i<=nterms;i++) coef[i]=f1.coef[i];
        return *this;
    } 

    complex<double> evaluate(complex<double> z);
    void compose(laurent f1,laurent f2);
    void inverse(laurent f1);
    laurent shift(complex<double> c);
    void constant(long n,complex<double> c);
    void identity(long n);
    void geometric(std::complex<double>, long int);
    void approximate_frac_pow(double alpha, double c,long n);

    void approximate_vertical_slit_ch_zip(double dt,double dx,long n);
    void approximate_tilted_slit_ch_zip(double alpha,double xl,double xr,
        long n);
    void approximate_arc_ch_zip(double a,double b, long n);

    void approximate_vertical_slit_ch_unzip(double dt,double dx,long n);
    void approximate_tilted_slit_ch_unzip(double alpha,double xl,double xr,
        long n);
    void approximate_arc_ch_unzip(double a,double b, long n);
    void approximate_vertical_slit_rad_zip(double dt,double dx,long n);

    void approximate_vertical_slit_half_rad_zip(double dt,double dx,long n);
    void approximate_tilted_slit_half_rad_zip(double dx,double alpha,
        double xl,double xr,complex<double> ww,long n);
};

inline laurent operator*(laurent f1,complex<double> c) 
{
    f1*=c;
    return f1;
}

inline laurent operator*(complex<double> c,laurent f1) 
{
    f1*=c;
    return f1;
}

#endif

///////////////////////////////////////////////////////////////
//  checkk checkk checkk checkk checkk checkk checkk checkk  //
///////////////////////////////////////////////////////////////
