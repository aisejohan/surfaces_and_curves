/*
 *	compute.h
 *
 * 	Copyright 2006 Johan de Jong
 *
 *	This file is part of Frobenius
 *
 *	Frobenius is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation; either version 2 of the License, or
 *	(at your option) any later version.
 *
 *	Frobenius is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with Frobenius; if not, write to the Free Software Foundation, 
 *	Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *									*/


struct exponents {
	unsigned int e1;
	unsigned int e2;
	unsigned int e3;
	unsigned int e4;
	unsigned int e5;
};

struct base_change {
	struct polynomial bc1;
	struct polynomial bc2;
	struct polynomial bc3;
	struct polynomial bc4;
	struct polynomial bc5;
};

struct lijst {
	struct base_change **BC;
	struct exponents **ee;
	struct polynomial **ff;
	unsigned int len;
};

/* Global variables used in other files */
extern struct lijst G;
extern struct polynomial myf;

/* Functions used in other files. */
int setup(int silent);
void deallocate_GVMnewMMold(void );
