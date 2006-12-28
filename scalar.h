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

mpz_t modulus;
mpz_t prime;
mpz_t temp;

void setup_scalars(void);

void printmscalar(mscalar a);

#ifndef PROFILER
#define make_scalar(a)		mpz_init(a)
#define free_scalar(a)		mpz_clear(a)
#else
void make_scalar(mscalar a);
void free_scalar(mscalar a);
#endif

#define valuation(x)		mpz_remove(temp,x,prime)

#define sc_add(a,b,c)		{mpz_add(c,a,b);mpz_mod(c,c,modulus);}

#define sc_mult(a,b,c)		{mpz_mul(c,a,b);mpz_mod(c,c,modulus);}

#define sc_imult(a,b,c)		{mpz_mul_si(c,b,(long)a);mpz_mod(c,c,modulus);}

#define sc_inv(a,b)		mpz_invert(b,a,modulus)

void sc_div(mscalar a, mscalar b, mscalar c);	

#define div_p(a)		mpz_divexact_ui(a,a,(unsigned long)p)

#define sc_add_replace(a,b)	{mpz_add(b,a,b);mpz_mod(b,b,modulus);}

#define sc_mult_replace(a,b)	{mpz_mul(b,a,b);mpz_mod(b,b,modulus);}

#define sc_imult_replace(a,b)	{mpz_mul_si(b,b,(long)a);mpz_mod(b,b,modulus);}

#define sc_zero(a)		mpz_set_ui(a,0)

#define sc_one(a)		mpz_set_ui(a,1)

#define sc_copy(a,b)		mpz_set(b,a)

#define sc_negate(a)		{mpz_neg(a,a);mpz_mod(a,a,modulus);}

#define ito_sc(a,b)		{mpz_set_si(b,(long)a);mpz_mod(b,b,modulus);}

#define sc_is_zero(a)		mpz_divisible_p(a,modulus)
