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

int main(void )
{
	int i, j, k, retry, extra;
	int c;
	mscalar cc, cc1, cc2;
	int *gap;
	int e;
	int blen1, blen2, blen3;
	struct term **basis1, **basis2, **basis3;
	int glen1, glen2, glen3;
	struct term **gens1, **gens2, **gens3;
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
	make_scalar(cc);
	make_scalar(cc1);
	make_scalar(cc2);
	
#ifdef KIJKEN
	printf("Debug is set! To unset do not define KIJKEN.\n");
#endif
	if (d % 2 != 0) {
		printf("Error: degree should be even.\n");
		exit(1);
	}

	if (d < 2*d1+2*d2+2*d3) {
		printf("Error: degree should be >= 2(d1+d2+d3).\n");
		exit(1);
	}

	if (p == 2) {
		printf("Error: p cannot be 2.\n");
		exit(1);
	}

	/* Setup the scalars. */
	setup_scalars();

	/* Seed the randomness. */
	set_seed(0);

	retry = 1;
	while (retry == 1) {
		while (retry == 1) {
			retry = setup(0);
		}

		gap = find_gap();
		blen1 = char_0(d/2 - d1 - d2 - d3, gap);
		blen2 = char_0(3*d/2 - d1 - d2 - d3, gap);
		blen3 = char_0(5*d/2 - d1 - d2 - d3, gap);
		glen1 = char_p(d/2 - d1 - d2 - d3);
		glen2 = char_p(3*d/2 - d1 - d2 - d3);
		glen3 = char_p(5*d/2 - d1 - d2 - d3);
		printf("Lengths: %d %d %d %d %d %d.\n",
			blen1, glen1, blen2, glen2, blen3, glen3);
		basis1 = char_0_basis(d/2 - d1 - d2 - d3, blen1, gap);
		basis2 = char_0_basis(3*d/2 - d1 - d2 - d3, blen2, gap);
		basis3 = char_0_basis(5*d/2 - d1 - d2 - d3, blen3, gap);
		gens1 = char_p_generators(d/2 - d1 - d2 - d3, glen1);
		gens2 = char_p_generators(3*d/2 - d1 - d2 - d3, glen2);
		gens3 = char_p_generators(5*d/2 - d1 - d2 - d3, glen3);
		printf("For %d = 5*d/2-d1-d2-d3 you get %d in char 0.\n",
					5*d/2 - d1 - d2 - d3, blen3);
		print_terms(basis3, blen3);
		printf("For %d = 3*d/2-d1-d2-d3 you get %d in char 0.\n",
					3*d/2 - d1 - d2 - d3, blen2);
		print_terms(basis2, blen2);
		printf("For %d = d/2-d1-d2-d3 you get %d in char 0.\n",
					d/2 - d1 - d2 - d3, blen1);
		print_terms(basis1, blen1);

		if (p == 3) {
			e = __extra(5*d/2 - d1 - d2 - d3, gap);
		} else {
			e = 0;
		}
		matrix = gens_to_basis(blen1, basis1, blen2, basis2,
			blen3, basis3, glen1, gens1, glen2, gens2,
			glen3, gens3, &e);
		e -= clean_matrix(glen1 + glen2 + glen3,
				blen1 + blen2 + blen3, matrix);
		printf("Extra powers of p used %d.\n", e);
	}

	/* Initialize fmatrix */
	fmatrix = (mscalar **)malloc((blen1+blen2+blen3)*sizeof(mscalar *));
	if (!fmatrix) {
		perror("Malloc failed!");
		exit(1);
	}
	
	/* Initialize fbasis. */
	fbasis = (struct polynomial ***)
		malloc((blen1 + blen2 + blen3)*sizeof(struct polynomial **));
	if (!fbasis) {
		perror("Malloc failed!");
		exit(1);
	}

	/* Initialize hhh. */
	hhh = (struct polynomial ***)
		malloc(3*sizeof(struct polynomial **));
	if (!hhh) {
		perror("Malloc failed!");
		exit(1);
	}

	/* Initialize extra. */
	extra = 0;
	for (i = 0; i <= q; i++) {
		j = (5 + 2*i)*p - 2;
		c = -i - 2;
		while (j > 0) {
			c += ivaluation(j);
			j -= 2;
		}
		if (c > extra) extra = c;
	}
	printf("The invariant extra is equal to %d.\n", extra);

	/* Initialize bb which is going to be equal to
	 * 	p^i Delta^i p^2 (x1...x3)^(p-1)
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
	 * This is the case i=0,j=5/2 of expansion in the file	*
	 * short_explanation.					*
	 * Note (i+j-1 choose j-1) is 1 in this case. 		*/
	sc_one(cc);
	hhh[0] = copy_pol_star(cc, bb);
	for (i = 0; i + 1 <= blen3; i++) {
		T.degree = p*(5*d/2 - d1 - d2 - d3);
		make_term(&T.leading);
		sc_one(T.leading->c);
		T.leading->n1 = p*basis3[i]->n1;
		T.leading->n2 = p*basis3[i]->n2;
		T.leading->n3 = p*basis3[i]->n3;
		T.leading->next = NULL;
		fbasis[i] = split_up(&T);
	}
	
	/* This is the case i=0,j=2 of expansion in the file	*
	 * short_explanation.					*
	 * Note (i+j-1 choose j-1) is 1 in this case. 		*/
	sc_one(cc);
	hhh[1] = copy_pol_star(cc, bb);
	for (i = 0; i + 1 <= blen2; i++) {
		T.degree = p*(3*d/2 - d1 - d2 - d3);
		make_term(&T.leading);
		sc_one(T.leading->c);
		T.leading->n1 = p*basis2[i]->n1;
		T.leading->n2 = p*basis2[i]->n2;
		T.leading->n3 = p*basis2[i]->n3;
		T.leading->next = NULL;
		fbasis[blen3 + i] = split_up(&T);
	}

	/* This is the case i=0,j=1 of expansion in the file	*
	 * short_explanation.					*
	 * Note (i+j-1 choose j-1) is 1 in this case. 		*/
	sc_one(cc);
	hhh[2] = copy_pol_star(cc, bb);
	for (i = 0; i + 1 <= blen1; i++) {
		T.degree = p*(1*d/2 - d1 - d2 - d3);
		make_term(&T.leading);
		sc_one(T.leading->c);
		T.leading->n1 = p*basis1[i]->n1;
		T.leading->n2 = p*basis1[i]->n2;
		T.leading->n3 = p*basis1[i]->n3;
		T.leading->next = NULL;
		fbasis[blen3 + blen2 + i] = split_up(&T);
	}

	/* This actually computes p*Delta */
	Delta = compute_delta();
	dd = split_up(&Delta);
	for (i = 1; i <= q; i++) {
		/* Compute next version of bb which is
		 * 	p^i Delta^i p^2 (x1...x4)^(p-1) 
		 * in split form. */
		printf("Start computing %d^%d Delta^%d"
			" (x1x2x3)^%d... ", p, 2 + i, i, p - 1);
		fflush(stdout);
		hh = mult_split(dd, bb);
		free_star(bb);
		free(bb);
		bb = hh;
		hh = NULL;
		printf("Done.\n");

		/* Highest degree and term is first basis element. 	*
		 * This is the case j=5/2,i=i of the file		*
		 * short_explanation.					*/
		sc_one(cc1);
		sc_one(cc2);
		sc_one(cc);
		c = i;
		while (c > 0) {
			k = 2*c;
			sc_imult_replace(k, cc1);
			k = 5 + 2*i - 2*c;
			sc_imult_replace(k, cc2);
			c--;
		}
		sc_div(cc2, cc1, cc);
		/* Note the above equals (i+j-1 choose j-1) in this case. */
		aa = copy_pol_star(cc, bb);
		merge_add_split(&(hhh[0]), aa);

		/* This is the case j=3/2,i=i of the file		*
		 * short_explanation.					*/
		sc_one(cc1);
		sc_one(cc2);
		c = i;
		while (c > 0) {
			sc_imult_replace(2*c, cc1);
			sc_imult_replace(3 + 2*i - 2*c, cc2);
			c--;
		}
		sc_div(cc2, cc1, cc);
		/* Note the above equals (i+j-1 choose j-1) in this case. */
		aa = copy_pol_star(cc, bb);
		merge_add_split(&(hhh[1]), aa);

		/* This is the case j=1,i=i of the file			*
		 * short_explanation.					*/
		sc_one(cc1);
		sc_one(cc2);
		c = i;
		while (c > 0) {
			sc_imult_replace(2*c, cc1);
			sc_imult_replace(1 + 2*i - 2*c, cc2);
			c--;
		}
		sc_div(cc2, cc1, cc);
		/* Note the abve equals (i+j-1 choose j-1) in this case. */
		aa = copy_pol_star(cc, bb);
		merge_add_split(&(hhh[2]), aa);
	}
	free_star(bb);
	free(bb);
	free_star(dd);
	free(dd);

	printf("Start computing aa.\n");
	for (j = 0; j + 1 <= blen3 + blen2 + blen1; j++) {
		if (j + 1 <= blen3) c = 0;
		else if (j + 1 <= blen3 + blen2) c = 1;
		else c = 2;
		printf("%d", j + 1);
		fflush(stdout);
		hh = mult_split(fbasis[j], hhh[c]);
		printf(" ");
		fflush(stdout);
		aa = all_the_way_split(hh);
		column = coefficients(aa, glen1, gens1, glen2, gens2,
								glen3, gens3);
		if (!column) {
			printf("FIXME!\n");
			exit(1);
		} else {
			fmatrix[j] = column;
			column = NULL;
		}
	}
	printf("\n");

	extra -= clean_matrix(blen1 + blen2 + blen3,
				glen1 + glen2 + glen3, fmatrix);
	
	Fmatrix = prod_matrix(blen1 + blen2 + blen3, glen1 + glen2 + glen3,
				blen1 + blen2 + blen3, fmatrix, matrix);

	k = extra + e;
	k -= clean_matrix(blen1 + blen2 + blen3,
					blen1 + blen2 + blen3, Fmatrix);
	print_matrix(blen1 + blen2 + blen3, blen1 + blen2 + blen3, Fmatrix);

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
	free_star(hhh[2]);
	free(hhh[2]);
	free(hhh);
	/* Free G and myf. */
	free_tail(myf.leading);
	for (i = 0; i + 1 <= G.len; i++) {
		free_tail(G.BC[i]->bc1.leading);
		free_tail(G.BC[i]->bc2.leading);
		free_tail(G.BC[i]->bc3.leading);
		free_tail(G.BC[i]->bc4.leading);
		free_tail(G.ff[i]->leading);
	}
	free(gap);
	free_list_terms(basis1, blen1);
	free_list_terms(basis2, blen2);
	free_list_terms(basis3, blen3);
	free_list_terms(gens1, glen1);
	free_list_terms(gens2, glen2);
	free_list_terms(gens3, glen3);

	/* Free fbasis. */
	for (i = 0; i + 1 <= blen1 + blen2 + blen3; i++) {
		free_star(fbasis[i]);
		free(fbasis[i]);
	}
	free(fbasis);
	/* Free matrices */
	free_matrix(glen1 + glen2 + glen3, blen1 + blen2 + blen3, matrix);
	free_matrix(blen1 + blen2 + blen3, glen1 + glen2 + glen3, fmatrix);
	free_matrix(blen1 + blen2 + blen3, blen1 + blen2 + blen3, Fmatrix);
	free_scalar(cc);
	free_scalar(cc1);
	free_scalar(cc2);
	/********************************************************
	 * End Neurotic freeing. 				*
	 ********************************************************/

	return(0);
}
