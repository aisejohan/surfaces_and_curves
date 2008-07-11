/*
 *	scalar.c
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

#include <stdio.h>
#include <stdlib.h>
#include "data.h"
#include "scalar.h"

extern unsigned long TT[4];
extern unsigned long SSS[4];

#ifdef KIJKEN
void empty_function(void)
{
	printf("This just for breakpoints when debugging.\n");
}
#endif

/* Only called once. */
void setup_scalars(void)
{
	printf("No set up needed.\n");
}

/*
unsigned long int valuation(mscalar x)
{
	return(mpz_remove(temp,x,prime));
}
*/

/*
void sc_add(mscalar a, mscalar b, mscalar c)
{
	mpz_add(c,a,b);
	mpz_mod(c,c,modulus);
}
*/

/*
void sc_mult(mscalar a, mscalar b, mscalar c)
{
	mpz_mul(c,a,b);
	mpz_mod(c,c,modulus);
}
*/

/*
void sc_imult(int a, mscalar b, mscalar c)
{
	mpz_mul_si(c,b,(long) a);
	mpz_mod(c,c,modulus);
}
*/

/*
void sc_inv(mscalar a, mscalar b)
{
	mpz_invert(b,a,modulus);
}
*/

/* Divides a by b. If b is not a unit then this assumes 	*
 * valuation(a) >= valuation(b), and the result is lifted	*
 * to an integer mod p^r.					*/
/* Does not destroy a and b.					*/
void sc_div(mscalar a, mscalar b, mscalar c)
{
	int e;

	e = VAL4(b);
	SET4(TT, b);
	SET4(SSS, a);
	while (e >= 64) {
		TT[0] = TT[1];
		TT[1] = TT[2];
		TT[2] = TT[3];
		TT[3] = 0;
		SSS[0] = SSS[1];
		SSS[1] = SSS[2];
		SSS[2] = SSS[3];
		SSS[3] = 0;
		e = e - 64;
	}
	DIV4(TT, e);
	DIV4(SSS, e);
	INV4(TT);
	MUL4(c, SSS, TT);
}

/* Divides a by p. The assumption is that this can be done. */
/*
void div_p(mscalar a)
{
	mpz_divexact_ui(a,a,(unsigned long) p);
}
*/

/*
void sc_add_replace(mscalar a, mscalar b)
{
	sc_add(a,b,b);
}
*/

/*
void sc_mult_replace(mscalar a, mscalar b)
{
	sc_mult(a,b,b);
}
*/

/*
void sc_imult_replace(int a, mscalar b)
{
	sc_imult(a,b,b);
}
*/

/*
void sc_zero(mscalar a)
{
	mpz_set_ui(a,(unsigned long) 0);
}
*/

/*
void sc_one(mscalar a)
{
	mpz_set_ui(a,(unsigned long) 1);
}
*/

/*
void sc_copy(mscalar a, mscalar b)
{
	mpz_set(b,a);
}
*/

/*
void sc_negate(mscalar a)
{
	mpz_neg(a,a);
	mpz_mod(a,a,modulus);
}
*/

/*
void ito_sc(int a, mscalar b)
{
	mpz_set_si(b,(long) a);
	mpz_mod(b,b,modulus);
}
*/

/*
int sc_is_zero(mscalar a)
{
	return(mpz_divisible_p(a,modulus));
}
*/
