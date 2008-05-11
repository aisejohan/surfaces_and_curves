/*
 *	basis.c
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
#include <stdlib.h>
#include <unistd.h>

#include "data.h"
#include "scalar.h"
#include "pol.h"
#include "grobner.h"
#include "helper.h"
#include "compute.h"
#include "delta.h"
#include "reduce.h"
#include "char_p_0.h"
#include "basis.h"

int compute_frobenius_matrix(void )
{
	int i, j, k, extra;
	int c;
	mscalar cc;
	int *gap;
	int e;
	int blen1, blen2;
	struct term **basis1, **basis2;
	int glen1, glen2;
	struct term **gens1, **gens2;
	mscalar **matrix, **fmatrix, **Fmatrix;
	mscalar *column;
	struct term *aaterm;
	struct polynomial Delta;
	struct polynomial T;
	struct polynomial **aa, **bb, **dd, **hh;
	struct polynomial ***hhh;
	struct polynomial ***fbasis;
	T.leading = NULL;
	Delta.leading = NULL;
	aaterm = NULL;
	setup_scalars();
	make_scalar(cc);
	
	gap = find_gap();
	blen1 = char_0(d - d1 - d2 - d3, gap);
	blen2 = char_0(2*d - d1 - d2 - d3, gap);
	glen1 = char_p(d - d1 - d2 - d3);
	glen2 = char_p(2*d - d1 - d2 - d3);
	printf("Lengths: %d %d %d %d.\n", blen1, glen1, blen2, glen2);
	basis1 = char_0_basis(d - d1 - d2 - d3, blen1, gap);
	basis2 = char_0_basis(2*d - d1 - d2 - d3, blen2, gap);
	gens1 = char_p_generators(d - d1 - d2 - d3, glen1);
	gens2 = char_p_generators(2*d - d1 - d2 - d3, glen2);
	printf("For %d = 2*d-d1-d2-d3 you get %d in char 0.\n",
						2*d - d1 - d2 - d3, blen2);
	print_terms(basis2, blen2);
	printf("For %d = d-d1-d2-d3 you get %d in char 0.\n",
						d - d1 - d2 - d3, blen1);
	print_terms(basis1, blen1);

	e = 0;
	matrix = gens_to_basis(blen1, basis1, blen2, basis2,
					glen1, gens1, glen2, gens2, &e);
	e -= clean_matrix(glen1 + glen2, blen1 + blen2, matrix);
	printf("Extra powers of p used %d.\n", e);

	/* Initialize fmatrix */
	fmatrix = (mscalar **)malloc((blen1 + blen2)*sizeof(mscalar *));
	if (!fmatrix) {
		perror("Malloc failed!");
		exit(1);
	}
	
	/* Initialize fbasis. */
	fbasis = (struct polynomial ***)
		malloc((blen1 + blen2)*sizeof(struct polynomial **));
	if (!fbasis) {
		perror("Malloc failed!");
		exit(1);
	}

	/* Initialize hhh. */
	hhh = (struct polynomial ***)
		malloc(2*sizeof(struct polynomial **));
	if (!hhh) {
		perror("Malloc failed!");
		exit(1);
	}

	/* Initialize extra. */
	extra = rr - r;
	printf("The invariant extra is equal to %d.\n", extra);

	/* Initialize bb which is going to be equal to
	 * 	p^i Delta^i p^3 (x1...x3)^(p-1)
	 * at various stages. */
	T.degree = (p - 1)*(d1 + d2 + d3);
	make_term(&T.leading);
	sc_one(T.leading->c);
	for (k = 1; k <= extra + 2 + 0; k++) { /* Note extra powers of p. */
		sc_imult_replace(p, T.leading->c);
	}
	T.leading->n1 = p - 1;
	T.leading->n2 = p - 1;
	T.leading->n3 = p - 1;
	T.leading->next = NULL;
	bb = split_up(&T);

	/* Highest degree and term is first basis element. 	*
	 * This is the case i=0,j=2 of expansion in the file	*
	 * short_explanation.					*
	 * Note (i+j-1 choose j-1) is 1 in this case. 		*/
	sc_one(cc);
	hhh[0] = copy_pol_star(cc, bb);
	for (i = 0; i + 1 <= blen2; i++) {
		T.degree = p*(2*d - d1 - d2 - d3);
		make_term(&T.leading);
		sc_one(T.leading->c);
		T.leading->n1 = p*basis2[i]->n1;
		T.leading->n2 = p*basis2[i]->n2;
		T.leading->n3 = p*basis2[i]->n3;
		T.leading->next = NULL;
		fbasis[i] = split_up(&T);
	}
	
	/* This is the case i=0,j=1 of expansion in the file	*
	 * short_explanation.					*
	 * Note (i+j-1 choose j-1) is 1 in this case. 		*/
	sc_one(cc);
	hhh[1] = copy_pol_star(cc, bb);
	for (i = 0; i + 1 <= blen1; i++) {
		T.degree = p*(1*d - d1 - d2 - d3);
		make_term(&T.leading);
		sc_one(T.leading->c);
		T.leading->n1 = p*basis1[i]->n1;
		T.leading->n2 = p*basis1[i]->n2;
		T.leading->n3 = p*basis1[i]->n3;
		T.leading->next = NULL;
		fbasis[blen2 + i] = split_up(&T);
	}

	/* This actually computes p*Delta */
	Delta = compute_delta();
	dd = split_up(&Delta);
	for (i = 1; i <= q; i++) {
		/* Compute next version of bb which is
		 * 	p^i Delta^i p^3 (x1...x4)^(p-1) 
		 * in split form. */
		printf("Start computing %d^%d Delta^%d"
			" (x1x2x3x4)^%d... ", p, 3 + i, i, p - 1);
		fflush(stdout);
		hh = mult_split(dd, bb);
		free_star(bb);
		free(bb);
		bb = hh;
		hh = NULL;
		printf("Done.\n");

		/* Highest degree and term is first basis element. 	*
		 * This is the case j=2,i=i of the file			*
		 * short_explanation.					*/
		sc_one(cc);
		c = i + 1;
		sc_imult_replace(c, cc);
		/* Note (i+j-1 choose j-1) is (i+1) in this case. */
		aa = copy_pol_star(cc, bb);
		merge_add_split(&(hhh[0]), aa);	

		/* This is the case j=1,i=i of the file			*
		 * short_explanation.					*/
		sc_one(cc);
		/* Note (i+j-1 choose j-1) is 1 in this case. */
		aa = copy_pol_star(cc, bb);
		merge_add_split(&(hhh[1]), aa);
	}
	free_star(bb);
	free(bb);
	free_star(dd);
	free(dd);

	printf("Start computing aa.\n");
	for (j = 0; j + 1 <= blen2 + blen1; j++) {
		if (j + 1 <= blen2) c = 0;
		else c = 1;
		printf("%d", j+1);
		fflush(stdout);
		hh = mult_split(fbasis[j], hhh[c]);
		printf(" ");
		fflush(stdout);
		aa = all_the_way_split(hh);
		column = coefficients(aa, glen1, gens1, glen2, gens2);
		if (!column) {
			printf("FIXME!\n");
			exit(1);
		} else {
			fmatrix[j] = column;
			column = NULL;
		}
	}
	printf("\n");

	extra -= clean_matrix(blen1 + blen2, glen1 + glen2, fmatrix);
	
	Fmatrix = prod_matrix(blen1 + blen2, glen1 + glen2,
					blen1 + blen2, fmatrix, matrix);

	k = extra + e + 1;
	k -= clean_matrix(blen1 + blen2, blen1 + blen2, Fmatrix);
	print_matrix(blen1 + blen2, blen1 + blen2, Fmatrix);

	if (k == 0) {
		printf("This should be the matrix of frobenius!\n");
	} else {
		printf("This matrix times %d^(-%d)"
		" should be the matrix of frobenius.\n", p, k);
	}

	/************************************************
	 * Neurotic freeing continues even now.		*
	 * The reason for this is that it makes 	*
	 * it easier to detect memory leaks.		*
	 ************************************************/
	free_star(hhh[0]);
	free(hhh[0]);
	free_star(hhh[1]);
	free(hhh[1]);
	free(hhh);
	free(gap);
	free_list_terms(basis1, blen1);
	free_list_terms(basis2, blen2);
	free_list_terms(gens1, glen1);
	free_list_terms(gens2, glen2);

	/* Free fbasis. */
	for (i = 0; i + 1 <= blen1 + blen2; i++) {
		free_star(fbasis[i]);
		free(fbasis[i]);
	}
	free(fbasis);
	/* Free matrices */
	free_matrix(glen1 + glen2, blen1 + blen2, matrix);
	free_matrix(blen1 + blen2, glen1 + glen2, fmatrix);
	free_matrix(blen1 + blen2, blen1 + blen2, Fmatrix);
	free_scalar(cc);
	/********************************************************
	 * End Neurotic freeing. 				*
	 ********************************************************/

	return(0);
}
