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
void add_coefficients(struct polynomial **aa, int column)
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

void print_fmatrix(void)
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
};	

int main() 
{
	int i,j,k,retry,p_pow,precision=0;
	int c;
	struct term *aaterm;
	struct polynomial Delta;
	struct polynomial T;
	struct polynomial **aa, **bb, **dd, **hh;
	struct polynomial ***fbasis;
	T.leading = NULL;
	Delta.leading = NULL;
	aaterm = NULL;
	
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
	
	 /* This is the case i=0,j=2 of expansion in the file	*
	  * short_explanation.					*/
	for(i=0;i+1<=blen2;i++) {
		T.degree = 2*(p*d)-d1-d2-d3;
		make_term(&T.leading);
		sc_one(T.leading->c);
		T.leading->n1 = p*basis2[i]->n1 + p - 1;
		T.leading->n2 = p*basis2[i]->n2 + p - 1;
		T.leading->n3 = p*basis2[i]->n3 + p - 1;
		T.leading->next = NULL;
		fbasis[i] = split_up(&T);
		bb = copy_pol_star(fbasis[i]);
		aa = all_the_way_split(bb);
		/* Now aa[0] still has to be			*
		 * 	multiplied by p^2, and			*
		 * 	divided by p-parts of 2p-1,2p-2,...,3.	*
		 * For p>2 times p for 2 times 4		*/
		if(p>2) {
			c = p;
			times_int(c,aa[0]);
		} else {
			times_int(4,aa[0]);
		};
		/* For aa[1] we					*
		 * 	multiply by p^2				*
		 * 	divide by p-parts of 2p-1,...,2		*
		 * For all p we multiply by p.			*/
		c = p;
		times_int(c,aa[1]);
		add_coefficients(aa,i);
	};

	 /* This is the case i=0,j=1 of expansion in the file	*
	 * short_explanation.					*/
	for(i=0;i+1<=blen1;i++) {
		T.degree = p*d-d1-d2-d3;
		make_term(&T.leading);
		sc_one(T.leading->c);
		T.leading->n1 = p*basis1[i]->n1 + p - 1;
		T.leading->n2 = p*basis1[i]->n2 + p - 1;
		T.leading->n3 = p*basis1[i]->n3 + p - 1;
		T.leading->next = NULL;
		fbasis[blen2+i] = split_up(&T);
		bb = copy_pol_star(fbasis[blen2+i]);
		aa = all_the_way_split(bb);
		/* Now aa[0] still has to be			*
		 * 	multiplied by p^2, and			*
		 * 	divided by p-parts of p-1,p-2,...,3.	*
		 * For all p times p^2.				*/
		c = p*p;
		times_int(c,aa[0]);
		/* For aa[1] we					*
		 * 	multiply by p^2				*
		 * 	divide by p-parts of p-1,...,2		*
		 * For all p we get p^2				*/
		times_int(c,aa[1]);
		add_coefficients(aa,blen2+i);
	};

	Delta = compute_delta();
	dd = split_up(&Delta);
	bb = (struct polynomial **)malloc(sizeof(struct polynomial *));
	if(!bb) {
		perror("Malloc failed!");
		exit(1);
	};
	bb[0] = NULL;
	make_pol(&bb[0]);
	bb[0]->degree = 0;
	make_term(&bb[0]->leading);
	sc_one(bb[0]->leading->c);
	bb[0]->leading->n1 = 0;
	bb[0]->leading->n2 = 0;
	bb[0]->leading->n3 = 0;
	bb[0]->leading->next = NULL;
	for(i=1;i<=q;i++) {
		/* Compute Delta^i in split form. */
		printf("Starting computing Delta^%d... ",i); fflush(stdout);
		hh = mult_split(dd,bb);
		for(j=0;j<=(i-1)*p;j++) {
			free_tail(bb[j]->leading);
			free(bb[j]);
		};
		free(bb);
		bb = hh;
		hh = NULL;
		printf("Done.\n");

		/* Highest degree and term is first basis element. 	*
		 * This is the case j=2,i=i of the file			*
		 * short_explanation.					*/
		for(j=0;j+1<=blen2;j++) {
			printf("Starting computing hh... "); fflush(stdout);
			hh = mult_split(fbasis[j],bb);
			printf("Done.\n");
			printf("Starting computing aa... "); fflush(stdout);
			aa = all_the_way_split(hh);
			printf("Done.\n");
			/* Now aa[0] still has to be			*
			 * 	multiplied by p^2			*
			 * 	multiplied by p^i			*
			 * 	multiplied by (2-1+i choose i) 		*
			 * 	and divided by the p-parts of 		*
			 * 		(2+i)p-1,...,3.			*/
			p_pow = 2+i;
			c = i+1;
			while(c % p == 0) {
				c = c/p;
				p_pow++;
			};
			times_int(c,aa[0]);
			for(k=(2+i)*p-1;k>=3;k--) {
				p_pow -= ivaluation(k);
			};
			if(p_pow >= 0) {
				for(k=1;k<=p_pow;k++) {
					times_int(p,aa[0]);
				};
			} else {
				if(p_pow < precision) precision = p_pow;
				for(k=1;k<=-p_pow;k++) {
					aaterm = aa[0]->leading;
					while(aaterm) {
						if(valuation(aaterm->c) > 0) {
							div_p(aaterm->c);
							aaterm = aaterm->next;
						} else {
							printf("FIXME!\n");
							exit(1);
						};
					};
				};
			};
			/* For aa[1] and aa[2] we			*
			 * 	multiply by p^2				*
			 * 	multiply by p^i				*
			 * 	multiply by (2-1+i choose i)		*
			 * 	and divde by p-parts of 		*
			 * 		(2+i)p-1,...,2.			*/
			times_int(c,aa[1]);
			p_pow -= ivaluation(2);
			if(p_pow >= 0) {
				for(k=1;k<=p_pow;k++) {
					times_int(p,aa[1]);
					times_int(p,aa[2]);
				};
			} else {
				if(p_pow < precision) precision = p_pow;
				for(k=1;k<=-p_pow;k++) {
					aaterm =aa[1]->leading;
					while(aaterm) {
						if(valuation(aaterm->c) > 0) {
							div_p(aaterm->c);
							aaterm = aaterm->next;
						} else {
							printf("FIXME!\n");
							exit(1);
						};
					};
					aaterm =aa[2]->leading;
					while(aaterm) {
						if(valuation(aaterm->c) > 0) {
							div_p(aaterm->c);
							aaterm = aaterm->next;
						} else {
							printf("FIXME!\n");
							exit(1);
						};
					};
				};
			};
			add_coefficients(aa,j);
		};

		/* This is the case j=1,i=i of the file			*
		 * short_explanation.					*/
		for(j=0;j+1<=blen1;j++) {
			printf("Starting computing hh... "); fflush(stdout);
			hh = mult_split(fbasis[blen2+j],bb);
			printf("Done.\n");
			printf("Starting computing aa... "); fflush(stdout);
			aa = all_the_way_split(hh);
			printf("Done.\n");
			/* Now aa[0] still has to be			*
			 * 	multiplied by p^2			*
			 * 	multiplied by p^i			*
			 * 	multiplied by (1-1+i choose i)=1	*
			 * 	and divided by the p-parts of 		*
			 * 	   (1+i)p-1,...,3.			*/
			p_pow = 2+i;
			c = 1;
			times_int(c,aa[0]);
			for(k=(1+i)*p-1;k>=3;k--) {
				p_pow -= ivaluation(k);
			};
			if(p_pow >= 0) {
				for(k=1;k<=p_pow;k++) {
					times_int(p,aa[0]);
				};
			} else {
				if(p_pow < precision) precision = p_pow;
				for(k=1;k<=-p_pow;k++) {
					aaterm = aa[0]->leading;
					while(aaterm) {
						if(valuation(aaterm->c) > 0) {
							div_p(aaterm->c);
							aaterm = aaterm->next;
						} else {
							printf("FIXME!\n");
							exit(1);
						};
					};
				};
			};
			/* For aa[1] and aa[2] we 			*
			 * 	multiply by p^3				*
			 * 	multiply by p^i				*
			 * 	multiply by (1-1+i choose i)=1		*
			 * 	divide by p-parts of 			*
			 * 		(1+i)p-1,...,2.			*/
			times_int(c,aa[1]);
			p_pow -= ivaluation(2);
			if(p_pow >= 0) {
				for(k=1;k<=p_pow;k++) {
					times_int(p,aa[1]);
				};
			} else {
				if(p_pow < precision) precision = p_pow;
				for(k=1;k<=-p_pow;k++) {
					aaterm = aa[1]->leading;
					while(aaterm) {
						if(valuation(aaterm->c) > 0) {
							div_p(aaterm->c);
							aaterm = aaterm->next;
						} else {
							printf("FIXME!\n");
							exit(1);
						};
					};
					aaterm = aa[2]->leading;
					while(aaterm) {
						if(valuation(aaterm->c) > 0) {
							div_p(aaterm->c);
							aaterm = aaterm->next;
						} else {
							printf("FIXME!\n");
							exit(1);
						};
					};
				};
			};
			add_coefficients(aa,blen2+j);
		};

	
	};
	
	print_fmatrix();
	printf("The variable precision is %d.\n",precision);

	/************************************************
	 * Neurotic freeing continues even now.		*
	 * The reason for this is that it makes 	*
	 * it easier to detect memory leaks.		*/
	for(j=0;j<=q*p;j++) {
		free_tail(bb[j]->leading);
		free(bb[j]);
	};
	free(bb);
	for(j=0;j<=p;j++) {
		free_tail(dd[j]->leading);
		free(dd[j]);
	}
	free(dd);
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
	/********************************************************
	 * End Neurotic freeing. 				*
	 ********************************************************/
		

	return(0);
};
