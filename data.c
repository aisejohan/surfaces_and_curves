/*
 *	data.c
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
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "data.h"
#include "scalar.h"
#include "pol.h"
#include "helper.h"

int set_data(void )
{
	int c,i,j,extra,success;

	printf("Please enter the prime you want to use:\t\t ");
	success = 0;
	while (!success) {
		scanf("%d",&p);
		if (
			(p != 2) &&
			(p != 3) &&
			(p != 5) &&
			(p != 7) &&
			(p != 11) &&
			(p != 13) &&
			(p != 17) &&
			(p != 19)) {
				printf("Not a prime or too large (>19)."
				" Try again:\t ");
		} else {
				success = 1;
		}
	}
	printf("Please enter the three weights you want to use\n");
	printf("(enter three integers separated by spaces):\t ");
	success = 0;
	while (!success) {
		scanf("%d %d %d", &d1, &d2, &d3);
		if (
			(d1 > 500) ||
			(d2 > 500) ||
			(d3 > 500)) {
				printf("Weights too large (>500)."
				" Try again:\t\t ");
		} else {
			success = 1;
		}
	}
	printf("Please enter the degree of the curve:\t ");
	success = 0;
	while (!success) {
		scanf("%d", &d);
		if (
			(d > 2000) ||
			(count_sum(d) == 0) ||
			(hilbert(d - d1 - d2 - d3) == 0)) {
				printf("Degree too large (>2000),"
				" or no curve\nof that degree,"
				" or g = 0. Try again:\t\t ");
		} else {
			success = 1;
		}
	}
	printf("Please enter the number of terms in the\nexpansion");
	printf(" (try a small number first, like 5):\t ");
	success = 0;
	while (!success) {
		scanf("%d", &q);
		if (
			(q > 200)) {
				printf("Too large (>200). Try again:\t\t\t ");
		} else {
			success = 1;
		}
	}
	extra = 0;
	for (i = 0; i <= q; i++) {
		j = (2 + i)*p - 1;
		c = -i - 2;
		while (j > 0) {
			c += ivaluation(j);
			j--;
		}
		if (c > extra) extra = c;
	}
	r = extra + q + 10; /* For good luck. */
	if (p == 2) r = r + 32; /* More good luck for the even prime. */
	printf("Doing computations modulo %d^%d.\n", p, r);
	printf("Everything OK? (1=yes and 0=no):\t\t ");
	scanf("%d", &success);
	return(success);
}

int set_random(void )
{
	int success,i;

	printf("Please enter 0 to enter coefficients for the\nhypersurface"
	" and 1 for a random one:\t\t ");
	success = 0;
	while (!success) {
		scanf("%d", &i);
		if (
			(i != 0) &&
			(i != 1)) {
				printf("Reply of %d not understood."
				" Try again:\t\t ", i);
		} else {
			success = 1;
		}
	}
	return(i);
}

