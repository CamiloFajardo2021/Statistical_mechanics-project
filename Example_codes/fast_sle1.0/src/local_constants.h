// file name : src/local_constants.h
//


#ifndef __local_constants_hh
#define __local_constants_hh


#include "variable.h"
// variable.h contains constants that must be defined outside this file. 
// Usually variable.h will be created by a run_ script 
// What must be in this file depends on the model   
// 
// For nmodel==1 (SAW), it should contain the following #define statements:
// #define MAX_NPIVOT 1000 // maximum number of pivots in class walk
// #define ALGORITHM_NUMBER ?
//   ALGORITHM_NUMBER should be one of following to to specify the algorithm
//     7  square lattice, 7 nontrivial lattice symmetries 
//     2 square lattice, 2 nontrivial lattice symmetries 
//     5 hexagonal lattice, 5 nontrivial lattice symmetries 
//     11 triangular lattice, 11 nontrivial lattice symmetries 
//     1 Manhattan lattice
//     47 Cubic lattice, all 47 nontrivial lattice symmetries 
//     383 4d hypercubic lattice, all 383 nontrivial lattice symmetries 
//   The following two lines can be included in variable.h (either one or both)
//   Each of them will change slightly how the pivot algorithm is implemented.
//   They will not change the sequence of walks generated. Depending
//   on the parameters values they will make the program run slightly
//   faster or slower.
// #define USE_FIND_SEGMENT
// #define USE_THIRD

// For nmodel==6 (percolation), file should contain following 
// #define PERCOLATION_HEXAGONAL_LATTICE

// For nmodel==9 (Ising), file should contain one of 
// #define SQUARE_LATTICE
// #define ISING_TRIANGULAR_LATTICE

// kludge: if ALGORITHM_NUMBER not defined, it is not needed, but 
//     will get compiler errors from code that is never called.
#ifndef ALGORITHM_NUMBER
#define ALGORITHM_NUMBER 0
#endif

#ifdef PRIVATE
#define MAX_NRV 200   // maximum number of RV's 
#define MAX_NPFUNC 100 // maximum number of pfunction's 
#define MAX_NFDISTRIB 20   // maximum number of full distributions 
#define MAX_NFORCE_PT 100   // maximum number of force pts for SLE(kappa,rho)
#define MAX_NDOMAIN_PARMS 100 // max number of domain parms for percolation
#define MAX_NSAMPLES 10000  // max number of samples in class samples

#define MAX_CHAR 101 
// values of the exponent nu. Used to rescale things
#define NU_TWOD 0.75 
#define NU_THREED 0.588
#define NU_FOURD 0.5
// following used in routines that look for hitting points of various curves
#define HIT_TOLERANCE 1.e-8
#endif

// following depend on model and lattice
// NUM_SYM= number of nontrivial lattice symmetries, identity not counted 
//    So this is the order of the symmetry group minus 1. 
// NUM_SYM_USED= number of lattice symmetries used by pivot algorithm
//    It is not used at present.

#define NUM_SYM 1 
#define NUM_SYM_USED 1 

#if ALGORITHM_NUMBER==1
#define NUM_SYM 7 
#define NUM_SYM_USED 1 
#define SQUARE_LATTICE
#endif

#if ALGORITHM_NUMBER==2
#define NUM_SYM 7 
#define NUM_SYM_USED 2 
#define SQUARE_LATTICE
#endif

#if ALGORITHM_NUMBER==7
#define NUM_SYM 7 
#define NUM_SYM_USED 7 
#define SQUARE_LATTICE
#endif

#if ALGORITHM_NUMBER==5
#define NUM_SYM 5 
#define NUM_SYM_USED 5 
#define HEXAGONAL_TRIANGULAR_LATTICE
#define TWO_DIMENSIONS
#define COORD_NUM 3
#endif

#if ALGORITHM_NUMBER==11
#define NUM_SYM 11 
#define NUM_SYM_USED 11 
#define HEXAGONAL_TRIANGULAR_LATTICE
#define TWO_DIMENSIONS
#define COORD_NUM 6
#endif

#if ALGORITHM_NUMBER==47
#define NUM_SYM 47
#define NUM_SYM_USED 47 
#define CUBIC_LATTICE
#define THREE_DIMENSIONS
#define COORD_NUM 6
#endif

#if ALGORITHM_NUMBER==383
#define NUM_SYM 383
#define NUM_SYM_USED 383 
#define HYPER_LATTICE
#define FOUR_DIMENSIONS
#define COORD_NUM 8
#endif

#ifdef SQUARE_LATTICE
#define TWO_DIMENSIONS
#define COORD_NUM 4
#endif

#ifdef ISING_TRIANGULAR_LATTICE
#define NUM_SYM 0 
#define TWO_DIMENSIONS
#define COORD_NUM 6
#endif

#ifdef PERCOLATION_HEXAGONAL_LATTICE
#define NUM_SYM 0 
#define NUM_SYM_USED 0 
#define TWO_DIMENSIONS
#define COORD_NUM 3
#endif

#endif

///////////////////////////////////////////////////////////////
//  checkk checkk checkk checkk checkk checkk checkk checkk  //
///////////////////////////////////////////////////////////////
