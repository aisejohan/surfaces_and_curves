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
#include "helper.h"

static int extra;
static mpz_t *modulus;
static mpz_t temp;
static mpz_t prime;

/* Only called once. */
void setup_scalars(void)
{
	int c,i,j;

	/* Initialize extra. */
	extra=0;
	for(i=0;i<=q;i++) {
		j = (3+i)*p-1;
		c = -i-3;
		while (j > 0) {
			c += ivaluation(j);
			j--;
		}
		if (c > extra) extra = c;
	}
	printf("The invariant extra is equal to %d.\n",extra);

	mpz_init_set_ui(prime,(unsigned long) p);
	mpz_init(temp);
	
	modulus = (mpz_t *)malloc((r+extra)*sizeof(mpz_t));
	i=0;
	while (i < r + extra) {
		mpz_init(modulus[i]);
		j = r + extra - i;
		printf("Here %d.\n",j);
		mpz_ui_pow_ui(modulus[i], (unsigned long) p, (unsigned long) j);
		i++;
	}
}

#ifdef KIJKEN
void test_scalar(mscalar a)
{
	if (a->e > r) {
		printf("Valuation too big.\n");
		exit(1);
	}
	if (a->e == r) {
		if (mpz_cmp_si(a->i, 0) != 0) {
			printf("Zero incorrect.\n");
			exit(1);
		} else {
			return;
		}
	}
	if (mpz_cmp_si(a->i, 0) == 0) {
		printf("Power not >= r but coeff == 0.\n");
		exit(1);
	}
	if (mpz_remove(temp, a->i, prime)) {
		printf("Coefficient pos valuation!\n");
		exit(1);
	}
	if (mpz_cmp(a->i, modulus[extra + a->e]) > 0) {
		printf("a->i not reduced.\n");
		exit(1);
	}
	if (mpz_cmp_si(a->i, 0) < 0) {
		printf("a->i negative.\n");
		exit(1);
	}
}
#endif

void printmscalar(mscalar a)
{
#ifdef KIJKEN
	test_scalar(a);
#endif
	mpz_t s,t;

	mpz_init(s);
	mpz_init(t);

	if (a->e == r) {
		printf("0");
		return;
	}

	if (a->e) printf("p^%d * ", a->e);

	mpz_set(s, modulus[extra + a->e]);
	mpz_cdiv_q_ui(t, s, 2);
	if (mpz_cmp(a->i, t)>0) {
		mpz_sub(t, a->i, s);
		mpz_out_str(stdout, (int) 10, t);
	} else {
		mpz_out_str(stdout, (int) 10, a->i);
	}

	mpz_clear(s);
	mpz_clear(t);
	return;
}

void make_scalar(mscalar a)
{
	a->e = r;
	mpz_init_set_ui(a->i, 0);
}

void free_scalar(mscalar a)
{
	mpz_clear(a->i);
}

int valuation(mscalar x)
{
#ifdef KIJKEN
	test_scalar(x);
	if (x->e >= r) {printf("Valuation of zero!");exit(1);}
#endif
	return(x->e);
}

void sc_add(mscalar a, mscalar b, mscalar c)
{
#ifdef KIJKEN
	test_scalar(a);
	test_scalar(b);
	test_scalar(c);
#endif
	if (a->e < b->e) {
		mpz_mul(temp, b->i, modulus[extra + r - b->e + a->e]);
		mpz_add(temp, temp, a->i);
		c->e = a->e;
		mpz_mod(c->i, temp, modulus[extra + c->e]);
		return;
	}
	if (a->e > b->e) {
		mpz_mul(temp, a->i, modulus[extra + r - a->e + b->e]);
		mpz_add(temp, temp, b->i);
		c->e = b->e;
		mpz_mod(c->i, temp, modulus[extra + c->e]);
		return;
	}
	c->e = a->e;
	mpz_add(temp, a->i, b->i);
	c->e += mpz_remove(temp, temp, prime);
	if (c->e < r) {
		mpz_mod(c->i, temp, modulus[extra + c->e]);
		return;
	}
	c->e = r;
	mpz_set_ui(c->i, (unsigned long) 0);
	return;
}

void sc_mult(mscalar a, mscalar b, mscalar c)
{
#ifdef KIJKEN
	test_scalar(a);
	test_scalar(b);
	test_scalar(c);
#endif
	c->e = a->e + b->e;
	if (c->e < r) {
		mpz_mul(temp, a->i, b->i);
		mpz_mod(c->i, temp, modulus[extra + c->e]);
		return;
	}
	c->e = r;
	mpz_set_ui(c->i, (unsigned long) 0);
	return;
}

void sc_imult(int a, mscalar b, mscalar c)
{
#ifdef KIJKEN
	test_scalar(b);
	test_scalar(c);
#endif
	if (a) {
		mpz_mul_si(temp, b->i, (long) a);
		c->e = b->e + mpz_remove(temp, temp, prime);
		if (c->e < r) {
			mpz_mod(c->i, temp, modulus[extra + c->e]);
			return;
		}
	}
	c->e = r;
	mpz_set_ui(c->i, (unsigned long) 0);
	return;
}

void sc_inv(mscalar a, mscalar b)
{
#ifdef KIJKEN
	test_scalar(a);
	test_scalar(b);
	if (a->e != 0) {
		printf("Not a unit!");
		exit(1);
	}
#endif
	b->e = 0;
	mpz_invert(b->i, a->i, modulus[extra]);
}

/* Divides a by b. If b is not a unit then this assumes 		*
 * valuation(a) >= valuation(b), and the result is lifted		*
 * to an integer mod p^r. It is also assumed that a is not zero.	*/
/* Does not destroy a and b.						*/
void sc_div(mscalar a, mscalar b, mscalar c)
{
#ifdef KIJKEN
	test_scalar(a);
	test_scalar(b);
	test_scalar(c);
	if (a->e < b->e) {printf("Not divisible in sc_div.\n");exit(1);}
	if (a->e - b->e >= r) {printf("Zero divided by something.\n");exit(1);}
#endif
	mpz_invert(temp, b->i , modulus[extra + b->e]);
	mpz_mul(temp, temp, a->i);
	c->e = a->e - b->e;
	mpz_mod(c->i, temp, modulus[extra + c->e]);
	return;
}

/* Divides a by p. */
void div_p(mscalar a)
{
#ifdef KIJKEN
	test_scalar(a);
#endif
	a->e--;
	return;
}

void sc_add_replace(mscalar a, mscalar b)
{
	sc_add(a,b,b);
}


void sc_mult_replace(mscalar a, mscalar b)
{
	sc_mult(a,b,b);
}


void sc_imult_replace(int a, mscalar b)
{
	sc_imult(a,b,b);
}


void sc_zero(mscalar a)
{
#ifdef KIJKEN
	test_scalar(a);
#endif
	a->e = r;
	mpz_set_ui(a->i, (unsigned long) 0);
}


void sc_one(mscalar a)
{
#ifdef KIJKEN
	test_scalar(a);
#endif
	a->e = 0;
	mpz_set_ui(a->i, (unsigned long) 1);
}

void sc_copy(mscalar a, mscalar b)
{
#ifdef KIJKEN
	test_scalar(a);
	test_scalar(b);
#endif
	b->e = a->e;
	mpz_set(b->i, a->i);
}

void sc_negate(mscalar a)
{
#ifdef KIJKEN
	test_scalar(a);
#endif
	mpz_neg(a->i, a->i);
	mpz_mod(a->i, a->i, modulus[extra + a->e]);
}

void ito_sc(int a, mscalar b)
{
#ifdef KIJKEN
	test_scalar(b);
#endif
	if (a) {
		mpz_set_si(temp, (long) a);
		b->e = mpz_remove(temp, temp, prime);
		if (b->e < r) {
			mpz_mod(b->i, temp, modulus[extra + b->e]);
			return;
		}
	}
	b->e = r;
	mpz_set_ui(b->i, 0);
}

int sc_is_zero(mscalar a)
{
#ifdef KIJKEN
	test_scalar(a);
#endif
	return((a->e == r));
}
