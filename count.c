/*
 *	helper.c
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

#include <stddef.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "data.h"
#include "scalar.h"
#include "pol.h"
#include "helper.h"
#include "compute.h"



unsigned int count_points(int degree)
{
	int n, i, m, j, u, v, x, y, z, w, nr;
	int *a, *b, *c, *e, *powers;
	struct term *t;

	if (degree != 1) {
		printf("Counting points over extension fields not defined!\n");
		return 0;
	}

	n = number_terms(myf);
	a = (int *)malloc(n*(sizeof(int)));
	b = (int *)malloc(n*(sizeof(int)));
	c = (int *)malloc(n*(sizeof(int)));
	e = (int *)malloc(n*(sizeof(int)));
	t = myf.leading;
	i = 0;
	m = 0;
	while (t) {
		if (t->n1 > m) { m = t->n1; }
		if (t->n2 > m) { m = t->n2; }
		if (t->n3 > m) { m = t->n3; }
		if (t->n4 > m) { m = t->n4; }
		mpz_mod(temp,t->c,prime);
		a[i] = mpz_get_si(temp);
		i = i + 1;
		t = t->next;
	}

	/* Compute powers mod p */
	powers = (int *)malloc(p*(m + 1)*sizeof(int));
	i = 0;
	while (i < p) {
		powers[i] = 1;
		i = i + 1;
	}
	j = 1;
	while (j <= m) {
		i = 0;
		while (i < p) {
			powers[i + j*p] = (i * powers[i + (j - 1)*p]) % p;
			i = i + 1;
		}
		j = j + 1;
	}

	nr = 0;
	x = 0;
	while (x < p) {
		t = myf.leading;
		i = 0;
		while (t) {
			b[i] = (a[i] * powers[x + p*t->n1]) % p;
			i = i + 1;
			t = t->next;
		}
		y = 0;
		while (y < p) {
			t = myf.leading;
			i = 0;
			while (t) {
				c[i] = (b[i] * powers[y + p*t->n2]) % p;
				i = i + 1;
				t = t->next;
			}
			z = 0;
			while (z < p) {
				t = myf.leading;
				i = 0;
				while (t) {
					e[i] = (c[i] * powers[z + p*t->n3]) % p;
					i = i + 1;
					t = t->next;
				}
				w = 0;
				while (w < p) {
					u = 0;
					t = myf.leading;
					i = 0;
					while (t) {
						v = (e[i] * powers[w + p*t->n4]) % p;
						u = (u + v) % p;
						i = i + 1;
						t = t->next;
					}
					if (u == 0) {
						nr = nr + 1;
					}
					w = w + 1;
				}
				z = z + 1;
			}
			y = y + 1;
		}
		x = x + 1;
	}
	free(a);
	free(b);
	free(c);
	free(e);
	free(powers);
	return nr;
}
