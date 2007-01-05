/*
 *	reduce.c
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
#include "pol.h"
#include "helper.h"
#include "grobner.h"
#include "delta.h"
#include "compute.h"

/* Here f is replaced by its grobner reduction.			*
 * The degree of f should be jd-s for some j.			*
 * The resulting polynomial of degree (j-1)d-s is returned.	*/
struct polynomial one_step_down(struct polynomial *f)
{
	int i,j,k;
	mscalar c;
	struct polynomial T;
	struct polynomial **aa;
	struct base_change fBC;
	make_scalar(c);
	
	if(((f->degree+d1+d2+d3+d4) % d != 0) || (f->degree+d1+d2+d3+d4 == d)) {
		printf("Incorrect degree. Stop.");
		exit(1);
	};
	j = (f->degree + d1 + d2 + d3 + d4) / d;
	
	fBC.bc1.leading = NULL;
	fBC.bc1.degree = f->degree - (d - d1);
	fBC.bc2.leading = NULL;
	fBC.bc2.degree = f->degree - (d - d2);
	fBC.bc3.leading = NULL;
	fBC.bc3.degree = f->degree - (d - d3);
	fBC.bc4.leading = NULL;
	fBC.bc4.degree = f->degree - (d - d4);
	fBC.bc5.leading = NULL;
	fBC.bc5.degree = f->degree - d;
	
	aa = gen_division(f,G.len,G.ff);
	
	for(i=0; i+1 <= G.len; i++) {
if(aa[i]->leading) {
		times_int(-1,aa[i]); /* Sign! */
		if((G.BC[i]->bc1.leading) && (fBC.bc1.leading)) {
			T = pol_mult(*aa[i],G.BC[i]->bc1);
			merge_add(&(fBC.bc1),T);
		} else if (G.BC[i]->bc1.leading) {
			fBC.bc1 = pol_mult(*aa[i],G.BC[i]->bc1);
		};
		if((G.BC[i]->bc2.leading) && (fBC.bc2.leading)) {
			T = pol_mult(*aa[i],G.BC[i]->bc2);
			merge_add(&(fBC.bc2),T);
		} else if (G.BC[i]->bc2.leading) {
			fBC.bc2 = pol_mult(*aa[i],G.BC[i]->bc2);
		};
		if((G.BC[i]->bc3.leading) && (fBC.bc3.leading)) {
			T = pol_mult(*aa[i],G.BC[i]->bc3);
			merge_add(&(fBC.bc3),T);
		} else if (G.BC[i]->bc3.leading) {
			fBC.bc3 = pol_mult(*aa[i],G.BC[i]->bc3);
		};
		if((G.BC[i]->bc4.leading) && (fBC.bc4.leading)) {
			T = pol_mult(*aa[i],G.BC[i]->bc4);
			merge_add(&(fBC.bc4),T);
		} else if (G.BC[i]->bc4.leading) {
			fBC.bc4 = pol_mult(*aa[i],G.BC[i]->bc4);
		};
		if((G.BC[i]->bc5.leading) && (fBC.bc5.leading)) {
			T = pol_mult(*aa[i],G.BC[i]->bc5);
			merge_add(&(fBC.bc5),T);
		} else if (G.BC[i]->bc5.leading) {
			fBC.bc5 = pol_mult(*aa[i],G.BC[i]->bc5);
		};
};
	};
	
	/* Free aa up. */
	for(i=0;i+1<=G.len;i++) {
		free_tail(aa[i]->leading);
		free(aa[i]);
	};
	free(aa);

	/* Derivatives. */
	rep_deriv(&(fBC.bc1),1);
	rep_deriv(&(fBC.bc2),2);
	rep_deriv(&(fBC.bc3),3);
	rep_deriv(&(fBC.bc4),4);
	
	/* Divide fBC.bci by p-primary part of j-1. 	*
	 * and multiply fBC.bc5 by p-part of j-1.	*/
	k = 1;
	i = j-1;
	/* Note that j is not 1, so i is not 0.		*/
	while(i % p == 0) {
		i = i/p;
		k = p*k;
	};
	/* c becomes the inverse of i */
	ito_sc(i,c);
	sc_inv(c,c);

	times_scalar(c,&(fBC.bc1));
	times_scalar(c,&(fBC.bc2));
	times_scalar(c,&(fBC.bc3));
	times_scalar(c,&(fBC.bc4));

	ito_sc(k,c);
	times_scalar(c,&(fBC.bc5));

	/* Adding up to get the result. */	
	merge_add(&(fBC.bc5), fBC.bc4);
	merge_add(&(fBC.bc5), fBC.bc3);
	merge_add(&(fBC.bc5), fBC.bc2);
	merge_add(&(fBC.bc5), fBC.bc1);

	free_scalar(c);
	return(fBC.bc5);
};

#if 0
/* This returns the complete reduction and destroys f.		*/
struct polynomial **all_the_way(struct polynomial *f)
{
	struct polynomial tmp;
	struct polynomial *g;
	struct polynomial **uit;
	tmp.leading=NULL;
	
	while(f->degree + d1+d2+d3+d4 > 3*d) {
		tmp = one_step_down(f);
#ifdef KIJKEN
		if(f->leading) {
			printf("Error: f did not reduce to zero!");
			exit(1);
		};
#endif
		f->degree = tmp.degree;
		f->leading = tmp.leading;
		tmp.degree = 0;
		tmp.leading = NULL;
	};
	uit=(struct polynomial **)malloc(3*sizeof(struct polynomial *));
	if(!uit) {
		perror("Malloc failed!");
		exit(1);
	};
	uit[0] = NULL;
	uit[1] = NULL;
	uit[2] = NULL;
	make_pol(&uit[0]);
	make_pol(&uit[1]);
	make_pol(&uit[2]);
	
	uit[0]->degree = 3*d - (d1+d2+d3+d4);
	uit[1]->degree = 2*d - (d1+d2+d3+d4);
	uit[2]->degree = d - (d1+d2+d3+d4);

	if(f->degree + d1+d2+d3+d4 == 3*d) {
		g = NULL;
		make_pol(&g);
		*g = one_step_down(f);
		uit[0] = f;
		*uit[2] = one_step_down(g);
		uit[1] = g;
	} else if(f->degree + d1+d2+d3+d4 == 2*d) {
		*uit[2] = one_step_down(f);
		uit[1] = f;
	} else if(f->degree + d1+d2+d3+d4 == d) {
		uit[2] = f;
	} else {
		printf("Wrong degree again!");
		exit(1);
	};
	
	return(uit);
};
#endif

/* Deals with a split up polynomial and reduces all the		*
 * way down. Destroys bb. 					*
 * ALTERNATIVE VERSION.						*/
struct polynomial **all_the_way_split(struct polynomial **bb)
{
	int i,j,ii,jj,k,c,tel;
	struct polynomial T;
	struct polynomial **aa,**cc;
	
	/* bb[0] has degree (jj)d-s */
	jj = 1 + bb[0]->degree/d;
	/* This means we have bb[0],...,bb[jj-1] */
	
	aa = (struct polynomial **)malloc(3*sizeof(struct polynomial *));
	if(!aa) {
		perror("Malloc failed!");
		exit(1);
	};
	aa[0] = NULL;
	aa[1] = NULL;
	aa[2] = NULL;
	make_pol(&aa[0]);
	make_pol(&aa[1]);
	make_pol(&aa[2]);

	aa[0]->degree = 3*d-d1-d2-d3-d4;
	aa[1]->degree = 2*d-d1-d2-d3-d4;
	aa[2]->degree = d-d1-d2-d3-d4;
	
	for(ii=0;ii+1<=jj;ii++) {
		
		/* bb[ii] has degree 	*
		 * (jj-ii)d-s = jd-s,	*
		 * so j=jj-ii 		*/
		j = (bb[ii]->degree+d1+d2+d3+d4)/d;

		/* In the one_step_down function (inside all_the_way)	*
		 * we multiply by the correct factor up to factors of	*
		 * p, so we have to correct for that here. We do this	*
		 * stupidly, mimicking what happens in one_down so we 	*
		 * don't make an error.					*
		 * Note that k is not cumulative (reset back to 1)	*/
		if(j>1) {
			/* This will have degree (j-1)d - s	*/
			T = one_step_down(bb[ii]);
			cc = split_up(&T);
			i = j-1;
			k = 1;
			while(i % p == 0) {
				i = i/p;
				k = p*k;
			};
			c = k;
			for(tel=1;tel+ii+1<=jj;tel++){
				times_int(c,bb[ii+tel]);
				merge_add(bb[ii+tel],*cc[tel-1]);
				free(cc[tel-1]);
			};
			free(cc);
		};
	};

	merge_add(aa[0],*bb[jj-3]);
	free(bb[jj-3]);
	merge_add(aa[1],*bb[jj-2]);
	free(bb[jj-2]);
	merge_add(aa[2],*bb[jj-1]);
	free(bb[jj-1]);
	free(bb);
	
	return(aa);
};
