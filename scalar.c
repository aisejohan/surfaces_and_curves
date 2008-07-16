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

static mpz_t *modulus;
static mpz_t temp;
static mpz_t prime;

/* Only called once. */
void setup_scalars(void)
{
	int c, i, j, extra;

	/* Initialize extra. */
	extra = 0;
	for (i = 0; i <= q; i++) {
		j = (3 + i)*p - 1;
		c = -i - 3;
		while (j > 0) {
			c += ivaluation(j);
			j--;
		}
		if (c > extra) extra = c;
	}
	rr = r + extra;
	printf("The invariant extra is equal to %d and rr is %d.\n", extra, rr);

	mpz_init_set_ui(prime, p);
	mpz_init(temp);
	
	modulus = (mpz_t *)malloc(rr*sizeof(mpz_t));
	i = 0;
	while (i < rr) {
		mpz_init(modulus[i]);
		j = rr - i;
		mpz_ui_pow_ui(modulus[i], p, j);
		i++;
	}
	return;
}

void close_scalars(void )
{
	int i;

	i = rr - 1;
	while (i >= 0) {
		mpz_clear(modulus[i]);
		i--;
	}
	free(modulus);
	mpz_clear(temp);
	mpz_clear(prime);
}

#ifdef KIJKEN
void test_scalar(mscalar a)
{
	if (a->e > rr) {
		printf("Valuation too big.\n");
		exit(1);
	}
	if (a->e == rr) {
		if (mpz_cmp_si(a->i, 0) != 0) {
			printf("Zero incorrect.\n");
			exit(1);
		} else {
			return;
		}
	}
	if (a->e < 0) {
		printf("Negative power!");
		exit(1);
	}
	if (mpz_cmp_si(a->i, 0) == 0) {
		printf("Power not >= rr but coeff == 0.\n");
		exit(1);
	}
	if (mpz_remove(temp, a->i, prime)) {
		printf("Coefficient pos valuation!\n");
		exit(1);
	}
	if (mpz_cmp(a->i, modulus[a->e]) > 0) {
		printf("a->i not reduced.\n");
		exit(1);
	}
	if (mpz_cmp_si(a->i, 0) < 0) {
		printf("a->i negative.\n");
		exit(1);
	}
	return;
}
#endif

void printmscalar(mscalar a)
{
#ifdef KIJKEN
	test_scalar(a);
#endif
	mpz_t s,t;


	if (a->e == rr) {
		printf("0");
		return;
	}
	mpz_init(s);
	mpz_init(t);

	if (a->e) printf("p^%d * ", a->e);

	mpz_set(s, modulus[a->e]);
	mpz_cdiv_q_ui(t, s, 2);
	if (mpz_cmp(a->i, t)>0) {
		mpz_sub(t, a->i, s);
		mpz_out_str(stdout, 10, t);
	} else {
		mpz_out_str(stdout, 10, a->i);
	}

	mpz_clear(s);
	mpz_clear(t);
	return;
}

#ifdef PROFILER
void make_scalar(mscalar a)
{
	mpz_init(a->i);
	return;
}

void free_scalar(mscalar a)
{
	mpz_clear(a->i);
	return;
}
#endif


#ifdef KIJKEN
unsigned int valuation(mscalar x)
{
	test_scalar(x);
	if (x->e == rr) {printf("Valuation of zero!");exit(1);}
	return(x->e);
}
#endif

static inline unsigned int my_p_remove(mpz_t b)
{
	unsigned int e=0;
	mpz_t x;

	if (mpz_sgn(b) == 0) return(rr);

	mpz_init(x);

	while (mpz_tdiv_q_ui(x, b, p) == 0) {
		mpz_set(b, x);
		e++;
	}

	mpz_clear(x);
	
	return(e);
}

void sc_add(mscalar a, mscalar b, mscalar c)
{
#ifdef KIJKEN
	test_scalar(a);
	test_scalar(b);
#endif
	if (a->e < b->e) {
		mpz_mul(temp, b->i, modulus[rr - b->e + a->e]);
		mpz_add(temp, temp, a->i);
		c->e = a->e;
		mpz_mod(c->i, temp, modulus[c->e]);
		return;
	}
	if (a->e > b->e) {
		mpz_mul(temp, a->i, modulus[rr - a->e + b->e]);
		mpz_add(temp, temp, b->i);
		c->e = b->e;
		mpz_mod(c->i, temp, modulus[c->e]);
		return;
	}
	c->e = a->e;
	mpz_add(temp, a->i, b->i);
	c->e += my_p_remove(temp);
	if (c->e < rr) {
		mpz_mod(c->i, temp, modulus[c->e]);
		return;
	}
	c->e = rr;
	mpz_set_ui(c->i, 0);
	return;
}

void sc_add_variant(mscalar a, mscalar b, mscalar c)
{
	if (a->e < b->e) {
		if (b->e < rr) {
			mpz_mul(temp, b->i, modulus[rr - b->e + a->e]);
			mpz_add(c->i, temp, a->i);
			c->e = a->e;
			return;
		}
		c->e = a->e;
		mpz_set(c->i, a->i);
		return;
	}
	if (a->e > b->e) {
		if (a->e < rr) {
			mpz_mul(temp, a->i, modulus[rr - a->e + b->e]);
			mpz_add(c->i, temp, b->i);
			c->e = b->e;
			return;
		}
		c->e = b->e;
		mpz_set(c->i, b->i);
		return;
	}
	c->e = a->e;
	mpz_add(c->i, a->i, b->i);
	return;
}

void clean_scalar(mscalar a)
{
	if (a->e >= rr) {
		a->e = rr;
		mpz_set_ui(a->i, 0);
		return;
	}
	mpz_mod(temp, a->i, modulus[a->e]);
	a->e += my_p_remove(temp);
	if (a->e < rr) {
		mpz_mod(a->i, temp, modulus[a->e]);
		return;
	}
	a->e = rr;
	mpz_set_ui(a->i, 0);
	return;
}

void sc_mult(mscalar a, mscalar b, mscalar c)
{
#ifdef KIJKEN
	test_scalar(a);
	test_scalar(b);
#endif
	c->e = a->e + b->e;
	if (c->e < rr) {
		mpz_mul(temp, a->i, b->i);
		mpz_mod(c->i, temp, modulus[c->e]);
		return;
	}
	c->e = rr;
	mpz_set_ui(c->i, 0);
	return;
}

void sc_imult(int a, mscalar b, mscalar c)
{
#ifdef KIJKEN
	test_scalar(b);
#endif

	mpz_mul_si(temp, b->i, a);
	c->e = b->e + my_p_remove(temp);
	if (c->e < rr) {
		mpz_mod(c->i, temp, modulus[c->e]);
		return;
	}
	c->e = rr;
	mpz_set_ui(c->i, 0);
	return;
}

void sc_inv(mscalar a, mscalar b)
{
#ifdef KIJKEN
	test_scalar(a);
	if (a->e != 0) {
		printf("Not a unit!");
		exit(1);
	}
#endif

	b->e = 0;
	mpz_invert(b->i, a->i, modulus[0]);
	return;
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
	if (a->e < b->e) {
		printf("Not divisible in sc_div.\n");
		exit(1);
	}
	if (a->e == rr) {
		printf("Zero divided by something.\n");
		exit(1);
	}
#endif

	mpz_invert(temp, b->i, modulus[b->e]);
	mpz_mul(temp, temp, a->i);
	c->e = a->e - b->e;
	mpz_mod(c->i, temp, modulus[c->e]);
	return;
}

/* Divides a by p^k where k is an integer. */
void div_p(int k, mscalar a)
{
#ifdef KIJKEN
	test_scalar(a);
	if (a->e - k < 0) {
		printf("Negative power.");
		exit(1);
	}
#endif

	if (a->e < rr) {
		a->e = a->e - k;
	}
	return;
}

void sc_zero(mscalar a)
{
	a->e = rr;
	mpz_set_ui(a->i, 0);
	return;
}


void sc_one(mscalar a)
{
	a->e = 0;
	mpz_set_ui(a->i, 1);
	return;
}

#ifdef KIJKEN
void sc_copy(mscalar a, mscalar b)
{
	test_scalar(a);
	b->e = a->e;
	mpz_set(b->i, a->i);
	return;
}
#endif

void sc_negate(mscalar a)
{
#ifdef KIJKEN
	test_scalar(a);
#endif

	mpz_neg(a->i, a->i);
	mpz_mod(a->i, a->i, modulus[a->e]);
	return;
}

void ito_sc(int a, mscalar b)
{
	mpz_set_si(temp, a);
	b->e = my_p_remove(temp);
	if (b->e < rr) {
		mpz_mod(b->i, temp, modulus[b->e]);
		return;
	}
	b->e = rr;
	mpz_set_ui(b->i, 0);
	return;
}


#ifdef KIJKEN
int sc_is_zero(mscalar a)
{
	test_scalar(a);
	return((a->e == rr));
}
#endif
