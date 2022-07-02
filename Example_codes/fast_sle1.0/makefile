###########################################################################
#
#    This file is part of fast_sle (version 1.0), which implements 
#    a fast algorithm for simulating the Schramm-Loewner evolution (SLE)  
#
#    Copyright (C) 2005 Tom Kennedy
#
#    fast_sle is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    fast_sle is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with fast_sle; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#    E-mail contact: tgk@math.arizona.edu
#    Ordinary mail: 
#       Tom Kennedy
#       Mathematics Department
#       University of Arizona
#       Tucson, ZA 85721, USA
#
###########################################################################


# makefile for fast sle public code distribution

###########################################################
# classes
###########################################################

conformal_maps.o: src/conformal_maps.h src/conformal_maps.c; \
 g++ -I/home/tgk/sprng/include -L/home/tgk/sprng/lib -O3 -Wall  \
 src/conformal_maps.c -c -o conformal_maps.o

laurent.o: src/laurent.h src/laurent.c; \
 g++ -I/home/tgk/sprng/include -L/home/tgk/sprng/lib -O3 -Wall  \
 src/laurent.c -c -o laurent.o

lib.o: src/lib.h src/lib.c; \
 g++ -I/home/tgk/sprng/include -L/home/tgk/sprng/lib -O3 -Wall  \
 src/lib.c -c -o lib.o

complex_walk.o: src/complex_walk.c src/complex_walk.h conformal_maps.o\
 src/random_number_generator.h;\
 g++ -I/home/tgk/sprng/include -L/home/tgk/sprng/lib -O3 -Wall  \
 src/complex_walk.c -c -o complex_walk.o


###########################################################
# programs
###########################################################

driver_sle: src/driver_sle.c complex_walk.o conformal_maps.o laurent.o lib.o\
 src/random_number_generator.h;\
	g++ -I/home/tgk/sprng/include -L/home/tgk/sprng/lib -O3 -Wall  \
src/driver_sle.c complex_walk.o conformal_maps.o laurent.o lib.o \
-lm -o driver_sle

