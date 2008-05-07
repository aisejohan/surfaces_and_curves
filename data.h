/*
 *	data.h
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

/* #define KIJKEN */

#define KLEINER	1
#define GELIJK	-1
#define GROTER	0

/* There are two orderings. LEX and REVLEX. Chose one by
 * commenting out one of the following two lines. */
#define LEX_ORDER
/* #define REVLEX_ORDER */

/* There are three strategies for reducing the polynomials
 * encountered during the computation. OLD, NEW and MIXED.
 * Setting MIXED_GROBNER is a good choice on average, but
 * for certain choices of degrees etc the combinations
 * 	LEX_ORDER + NEW_GROBNER
 * or
 * 	REVLEX_ORDER + OLD_GROBNER
 * can be much faster. */
#define NEW_GROBNER

#define d1	9
#define d2	13
#define d3	14
#define d	126
#define p	3
#define r	100		/* Exponent. */
#define q	20		/* Largest power of Delta. */

#define maxlength	512

/* This type will be used for our scalars. */
#include <gmp.h>
typedef struct scalar {
	int e;
	mpz_t i;
} mscalar[1];

struct term {
	mscalar c;
	unsigned int n1;
	unsigned int n2;
	unsigned int n3;
	struct term *next;
};

struct polynomial {
	unsigned int degree;
	struct term *leading;
};
