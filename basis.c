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

/* External variables. */
int blen1,blen2;
struct term **basis1,**basis2;
mscalar **fmatrix;


/* Takes the coefficients and frees aa. */
static void add_coefficients(struct polynomial **aa, int column)
{
	int row;
	struct term *aaterm;
	row = 0;
	aaterm = aa[0]->leading;
	while(aaterm) {
		if(
		(aaterm->n1 == basis2[row]->n1) &&
		(aaterm->n2 == basis2[row]->n2) &&
		(aaterm->n3 == basis2[row]->n3)) {
			sc_add_replace(aaterm->c,fmatrix[row][column]);
			row++;
			aaterm = aaterm->next;
		} else {
			row++;
		};
	};
	free_tail(aa[0]->leading);
	free(aa[0]);
	row = blen2;
	aaterm = aa[1]->leading;
	while(aaterm) {
		if(
		(aaterm->n1 == basis1[row-blen2]->n1) &&
		(aaterm->n2 == basis1[row-blen2]->n2) &&
		(aaterm->n3 == basis1[row-blen2]->n3)) {
			sc_add_replace(aaterm->c,fmatrix[row][column]);
			row++;
			aaterm = aaterm->next;
		} else {
			row++;
		};
	};
	free_tail(aa[1]->leading);
	free(aa[1]);
	free(aa);
	return;
}

static void print_fmatrix(void)
{
	int i,j;
	printf("[");
	for(i=0;i+1<=blen1+blen2;i++) {
		for(j=0;j+1<=blen1+blen2;j++) {
			printmscalar(fmatrix[i][j]);
			if(j+1 < blen1+blen2) printf(",");
		};
		if(i+1 < blen1+blen2) printf(";\\\n");
	};
	printf("]\n");
	return;
}

int main() 
{
	int i,j,k,retry,extra;
	int c;
	mscalar cc;
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
	
#ifdef KIJKEN
	printf("Debug is set! To unset do not define KIJKEN.\n");
#endif
	/* Setup the scalars. */
	setup_scalars();

	/* Seed the randomness. */
	set_seed(0);

	retry = 1;
	while(retry == 1) {
		while(retry == 1) {
			retry = setup();
		};

		if(d>=d1+d2+d3) {
			blen1=check_flatness(d-d1-d2-d3);
			printf("For %d = d-d1-d2-d3 you get %d\n",
					d-d1-d2-d3,blen1);
			if(blen1<=0) {
				retry = 1;
				sleep(10);
				/* Free up G and myf. */
				free_tail(myf.leading);
				for(i=0;i+1<=G.len;i++) {
					free_tail(G.BC[i]->bc1.leading);
					free_tail(G.BC[i]->bc2.leading);
					free_tail(G.BC[i]->bc3.leading);
					free_tail(G.BC[i]->bc4.leading);
					free_tail(G.ff[i]->leading);
				};
				for(i=0;i+1<=maxlength;i++) {
					free(G.BC[i]);
					free(G.ff[i]);
					free(G.ee[i]);
				};
				free(G.BC);
				free(G.ff);
				free(G.ee);
			} else {
				basis1 = find_basis(d-d1-d2-d3,blen1);
				for(i=0;i+1<=blen1;i++) {
					T.degree = d-d1-d2-d3;
					T.leading = basis1[i];
					print_pol(T);
					T.leading = NULL;
				};
				printf("\n");
			};
		};
		if((retry == 0) && (2*d>=d1+d2+d3)) {
			blen2=check_flatness(2*d-d1-d2-d3);
			printf("For %d = 2*d-d1-d2-d3 you get %d\n",
					2*d-d1-d2-d3,blen2);
			if(blen2<=0) {
				retry = 1;
				/* Free up G and myf. */
				free_tail(myf.leading);
				for(i=0;i+1<=G.len;i++) {
					free_tail(G.BC[i]->bc1.leading);
					free_tail(G.BC[i]->bc2.leading);
					free_tail(G.BC[i]->bc3.leading);
					free_tail(G.BC[i]->bc4.leading);
					free_tail(G.ff[i]->leading);
				};
				for(i=0;i+1<=maxlength;i++) {
					free(G.BC[i]);
					free(G.ff[i]);
					free(G.ee[i]);
				};
				free(G.BC);
				free(G.ff);
				free(G.ee);
				/* Free up basis1. */
				for(i=0;i+1<=blen1;i++) {
					free_term(basis1[i]);
				};
				free(basis1);
			} else {
				basis2 = find_basis(2*d-d1-d2-d3,blen2);
				for(i=0;i+1<=blen2;i++) {
					T.degree = 2*d-d1-d2-d3;
					T.leading = basis2[i];
					print_pol(T);
					T.leading = NULL;
				};
				printf("\n");
			};
		};
	};
	
	if(d < d1+d2+d3) {
		printf("Degree too small and p_g=0!\n");
		exit(0);
	};

	/* Initialize fmatrix */
	fmatrix = (mscalar **)malloc((blen1+blen2)*sizeof(mscalar *));
	if(!fmatrix) {
		perror("Malloc failed!");
		exit(1);
	};
	for(i=0;i+1<=blen1+blen2;i++) {
		fmatrix[i] = (mscalar *)
			malloc((blen1+blen2)*sizeof(mscalar));
		if(!fmatrix[i]) {
			perror("Malloc failed!");
			exit(1);
		};
	};
	for(i=0;i+1<=blen1+blen2;i++) {
		for(j=0;j+1<=blen1+blen2;j++) {
			make_scalar(fmatrix[i][j]);
			sc_zero(fmatrix[i][j]);
		};
	};

	/* Initialize fbasis. */
	fbasis = (struct polynomial ***)
		malloc((blen1+blen2)*sizeof(struct polynomial **));
	if(!fbasis) {
		perror("Malloc failed!");
		exit(1);
	};

	/* Initialize hhh. */
	hhh = (struct polynomial ***)
		malloc(2*sizeof(struct polynomial **));
	if(!hhh) {
		perror("Malloc failed!");
		exit(1);
	};

	/* Initialize extra. */
	extra=0;
	for(i=0;i<=q;i++) {
		j=(2+i)*p-1;
		c=-i-2;
		while (j > 0) {
			c += ivaluation(j);
			j--;
		}
		if (c > extra) extra = c;
	}
	printf("The invariant extra is equal to %d.\n",extra);

	/* Initialize bb which is going to be equal to
	 * 	p^i Delta^i p^3 (x1...x3)^(p-1)
	 * at various stages. */
	T.degree = (p-1)*(d1+d2+d3);
	make_term(&T.leading);
	sc_one(T.leading->c);
	for(k=1;k<=extra+2+0;k++) { /* Note extra powers of p. */
		sc_imult_replace(p,T.leading->c);
	}
	T.leading->n1 = p-1;
	T.leading->n2 = p-1;
	T.leading->n3 = p-1;
	T.leading->next = NULL;
	bb = split_up(&T);

	/* Highest degree and term is first basis element. 	*
	 * This is the case i=0,j=2 of expansion in the file	*
	 * short_explanation.					*
	 * Note (i+j-1 choose j-1) is 1 in this case. 		*/
	sc_one(cc);
	hhh[0] = copy_pol_star(cc,bb);
	for(i=0;i+1<=blen2;i++) {
		T.degree = p*(2*d-d1-d2-d3);
		make_term(&T.leading);
		sc_one(T.leading->c);
		T.leading->n1 = p*basis2[i]->n1;
		T.leading->n2 = p*basis2[i]->n2;
		T.leading->n3 = p*basis2[i]->n3;
		T.leading->next = NULL;
		fbasis[i] = split_up(&T);
	};
	
	/* This is the case i=0,j=1 of expansion in the file	*
	 * short_explanation.					*
	 * Note (i+j-1 choose j-1) is 1 in this case. 		*/
	sc_one(cc);
	hhh[1] = copy_pol_star(cc,bb);
	for(i=0;i+1<=blen1;i++) {
		T.degree = p*(1*d-d1-d2-d3);
		make_term(&T.leading);
		sc_one(T.leading->c);
		T.leading->n1 = p*basis1[i]->n1;
		T.leading->n2 = p*basis1[i]->n2;
		T.leading->n3 = p*basis1[i]->n3;
		T.leading->next = NULL;
		fbasis[blen2+i] = split_up(&T);
	};

	/* This actually computes p*Delta */
	Delta = compute_delta();
	dd = split_up(&Delta);
	for(i=1;i<=q;i++) {
		/* Compute next version of bb which is
		 * 	p^i Delta^i p^3 (x1...x4)^(p-1) 
		 * in split form. */
		printf("Start computing %d^%d Delta^%d"
			" (x1x2x3x4)^%d... ",p,3+i,i,p-1);
		fflush(stdout);
		hh = mult_split(dd,bb);
		free_star(bb);
		free(bb);
		bb = hh;
		hh = NULL;
		printf("Done.\n");

		/* Highest degree and term is first basis element. 	*
		 * This is the case j=2,i=i of the file			*
		 * short_explanation.					*/
		sc_one(cc);
		c=i+1;
		sc_imult_replace(c,cc);
		/* Note (i+j-1 choose j-1) is (i+1) in this case. */
		aa = copy_pol_star(cc,bb);
		merge_add_split(&(hhh[0]),aa);	


		/* This is the case j=1,i=i of the file			*
		 * short_explanation.					*/
		sc_one(cc);
		/* Note (i+j-1 choose j-1) is 1 in this case. */
		aa = copy_pol_star(cc,bb);
		merge_add_split(&(hhh[1]),aa);
	}

	printf("Start computing aa.\n");
	for(j=0;j+1<=blen2+blen1;j++) {
		if (j+1 <= blen2) c=0;
		else c=1;
		printf("%d",j+1); fflush(stdout);
		hh = mult_split(fbasis[j],hhh[c]);
		printf(" "); fflush(stdout);
		aa = all_the_way_split(hh);
		add_coefficients(aa,j);
	}
	printf("\n");

	k=extra+1;
	for(i=0;i+1<=blen1+blen2;i++) {
		for(j=0;j+1<=blen1+blen2;j++) {
			c=valuation(fmatrix[i][j]);
			if (!sc_is_zero(fmatrix[i][j]) && (c < k)) k = c;
		}
	}
	for(i=0;i+1<=blen1+blen2;i++) {
		for(j=0;j+1<=blen1+blen2;j++) {
			for(c=1;c<=k;c++) div_p(fmatrix[i][j]);
		}
	}

	print_fmatrix();
	if (k == extra+1 ) {
		printf("This should be the matrix of frobenius!\n");
	} else {
		printf("This matrix times %d^(-%d)"
		" should be the matrix of frobenius.\n",p,extra+1-k);
	}

	/************************************************
	 * Neurotic freeing continues even now.		*
	 * The reason for this is that it makes 	*
	 * it easier to detect memory leaks.		*
	 ************************************************/
	free_star(bb); free(bb);
	free_star(dd); free(dd);
	free_star(hhh[0]); free(hhh[0]);
	free_star(hhh[1]); free(hhh[1]);
	free(hhh);
	/* Free G and myf. */
	free_tail(myf.leading);
	for(i=0;i+1<=G.len;i++) {
		free_tail(G.BC[i]->bc1.leading);
		free_tail(G.BC[i]->bc2.leading);
		free_tail(G.BC[i]->bc3.leading);
		free_tail(G.BC[i]->bc4.leading);
		free_tail(G.ff[i]->leading);
	};
	for(i=0;i+1<=maxlength;i++) {
		free(G.BC[i]);
		free(G.ff[i]);
		free(G.ee[i]);
	};
	free(G.BC);
	free(G.ff);
	free(G.ee);
	/* Free basis1. */
	for(i=0;i+1<=blen1;i++) {
		free_term(basis1[i]);
	};
	free(basis1);
	/* Free basis2. */
	for(i=0;i+1<=blen2;i++) {
		free_term(basis2[i]);
	};
	free(basis2);
	/* Free fbasis. */
	for(i=0;i+1<=blen1+blen2;i++) {
		k = 1 + fbasis[i][0]->degree/d;
		for(j=0;j+1<=k;j++) {
			free_tail(fbasis[i][j]->leading);
			free(fbasis[i][j]);
		};
		free(fbasis[i]);
	};
	free(fbasis);
	/* Free fmatrix */
	for(i=0;i+1<=blen1+blen2;i++) {
		for(j=0;j+1<=blen1+blen2;j++) {
			free_scalar(fmatrix[i][j]);
		};
	};
	for(i=0;i+1<=blen1+blen2;i++) {
		free(fmatrix[i]);
	};
	free(fmatrix);
	free_scalar(cc);
	/********************************************************
	 * End Neurotic freeing. 				*
	 ********************************************************/


	return(0);
}
