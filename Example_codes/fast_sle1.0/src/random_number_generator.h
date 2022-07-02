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


// file name : src/random_number_generator.h
//
// Specifies which random number generator to use. 
// Choice is determined by value of RANDOM_NUMBER_GENERATOR, set below

// RANDOM_NUMBER_GENERATOR=1 uses drand48(). This is the easy one since it 
//  is probably already installed on your system.
//  Note however that it has been declared obsolete and so may disappear.
//
// RANDOM_NUMBER_GENERATOR=2 uses sprng(). This is what I have been using.
//  sprng stands for Scalable Parallel Random Number Generators.
//  It is available at http://sprng.cs.fsu.edu/ 
//  Depending on where you install it, you may need to modify the file 
//  'Makefile'. Makefile presently contains paths to where I have 
//  sprng installed (/home/tgk/sprng/include and /home/tgk/sprng/lib).
//  If you install it somewhere other than the places where the compiler
//  usually searches you will need to modify these paths appropriately.
// 
// For comments on how to add your own random number generator, see 
//  comments below in the block for RANDOM_NUMBER_GENERATOR==3


// NB : if you set the choice to 2 without sprng installed you will 
//  get a host of indecipherable compiler errors 

#ifndef __random_number_generator_hh
#define __random_number_generator_hh

#define RANDOM_NUMBER_GENERATOR 1

#if RANDOM_NUMBER_GENERATOR==1
    #define RNG drand48

    inline void initialize_random_number_generator(long iseed)
    {
        // NB : negative iseed does not use a system generated seed
        if (iseed<0) iseed=0; 
        srand48(iseed);
    }
#endif

#if RANDOM_NUMBER_GENERATOR==2
    #define SIMPLE_SPRNG 1 // used by sprng(); line must precede following 
    #include "sprng.h" // header file for sprng() random number generator
    // Note that we have violated "no nested includes" rule, but the goal 
    // was to get all the dependencies on the random number generator into 
    // this one file 
    #define RNG sprng

    inline void initialize_random_number_generator(long iseed)
    {
        // negative iseed make system generator the seed
        if (iseed<0) iseed=make_sprng_seed();  
        init_sprng(iseed,SPRNG_DEFAULT);	
    }
#endif

#if RANDOM_NUMBER_GENERATOR==3 
    // header file for your random number generator (if needed)
    #include "your_rng_header.h" 
    // your_rng should be the name of a function prototyped as
    //     double your_rng(void);
    // We use macro substitution to avoid an extra system call here
    #define RNG your_rng

    inline void initialize_random_number_generator(long iseed)
    {
        // code to seed your random number generator
        // We don't use a macro here - seeding is only done once
        printf("RANDOM_NUMBER_GENERATOR==3 is not implemented\n");
	exit(1);
    }
#endif

#endif


///////////////////////////////////////////////////////////////
//  checkk checkk checkk checkk checkk checkk checkk checkk  //
///////////////////////////////////////////////////////////////
