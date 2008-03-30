/*
 *	scalar.h
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

int rr;

void setup_scalars(void);

void printmscalar(mscalar a);

void make_scalar(mscalar *a);

void free_scalar(mscalar a);

unsigned int valuation(mscalar x);

void sc_add(mscalar a, mscalar b, mscalar c);
void sc_add_variant(mscalar a, mscalar b, mscalar c);
void clean_scalar(mscalar a);

void sc_mult(mscalar a, mscalar b, mscalar c);

void sc_imult(int a, mscalar b, mscalar c);

void sc_inv(mscalar a, mscalar b);

void sc_div(mscalar a, mscalar b, mscalar c);	

void div_p(int k, mscalar a);

#define	sc_add_replace(a,b)	sc_add(a,b,b)

#define	sc_mult_replace(a,b)	sc_mult(a,b,b)

#define	sc_imult_replace(a,b)	sc_imult(a,b,b)

void sc_zero(mscalar a);

void sc_one(mscalar a);

void sc_copy(mscalar a, mscalar b);

void sc_negate(mscalar a);

void ito_sc(int a, mscalar b);

int sc_is_zero(mscalar a);
