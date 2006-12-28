/*
 *	pol.c
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

/* Makes a term.				*/
void make_term(struct term **mon)
{
#ifdef KIJKEN
	if (*mon) {
		printf("Term already exists!\n");
		exit(1);
	};
#endif
	*mon = (struct term *) malloc (sizeof (struct term));
	if(!*mon) {
		perror("Malloc failed.");
		exit(1);
	};
	make_scalar((*mon)->c);
	(*mon)->next=NULL;
};

void make_pol(struct polynomial **f)
{
#ifdef KIJKEN
	if (*f) {
		printf("Polynomial already exists!\n");
		exit(1);
	};
#endif
	*f = (struct polynomial *) malloc (sizeof (struct polynomial));
	if(!*f) {
		perror("Malloc failed.");
		exit(1);
	};
	(*f)->leading=NULL;
};

/* Frees a term. */
void free_term(struct term *mon)
{
#ifdef KIJKEN
	if(!mon) {
		printf("free_term: term does not exist!\n");
		exit(1);
	};
#endif
	free_scalar(mon->c);
	free(mon);
	mon=NULL;
};

/* Copies data not pointer.						*/
/* sc_copy does allocation.						*/
void copy_term(struct term *mon1, struct term *mon2)
{
	sc_copy(mon1->c, mon2->c);
	mon2->n1 = mon1->n1;
	mon2->n2 = mon1->n2;
	mon2->n3 = mon1->n3;
	mon2->n4 = mon1->n4;
};

/* This frees memory starting with mon.					*
 * It does not set mon=NULL. This means you will/might			*
 * get a segfault if you use memory management wrong.			*/
void free_tail(struct term *mon)
{
	struct term *tmp;
	while(mon) {
		tmp = mon->next;
		free_term(mon);
		mon = tmp;
	};
	return;
};

/* Remove terms with zero coefficient from putative pol.	*
 * The result is a polynomial of the same degree as the input	*
 * which may have ``degree'' NOT zero! This is usefull in case	*
 * where we directly after using this add this to another pol	*
 * of the same degree as the input.				*/
void clean_pol(struct polynomial *pol)
{
	struct term *tmp;
	struct term **ptrterm;

	ptrterm = &(pol->leading);
	while(*ptrterm) {
		if(sc_is_zero((*ptrterm)->c)) {
			tmp = *ptrterm;
			*ptrterm = (*ptrterm)->next;
			free_term(tmp);
		} else {
			ptrterm = &((*ptrterm)->next);
		};
	};
	return;
};

/* Create data structure of the same length.		*/
void make_tail(struct term *mon, struct term **ptrterm)
{
	while(mon) {
		make_term(ptrterm);
		ptrterm = &((*ptrterm)->next);
		mon = mon->next;
	};
	*ptrterm = NULL;
};


/* Copies tail of pol.					*/
void copy_tail(struct term *mon, struct term **ptrterm)
{
	while(mon) {
		make_term(ptrterm);
		copy_term(mon,*ptrterm);
		ptrterm = &((*ptrterm)->next);
		mon = mon->next;
	};
};


/* Copies a polynomial.					*/
struct polynomial copy_pol(struct polynomial f)
{
	struct polynomial uit;
	uit.leading = NULL;
	
	uit.degree = f.degree;
	copy_tail(f.leading,&(uit.leading));
	return(uit);
};

/* This function assumes terms of the same degree. 	*
 * It compares the monomials not the coefficients.	*
 * Returns						*
 * 		GELIJK if equal				*
 * 		KLEINER if mon1 < mon2			*
 * 		GROTER if mon1 > mon2			*
 * 							*/
int kleiner(struct term *mon1, struct term *mon2)
{
	if(mon1->n1 != mon2->n1) return((mon1->n1 < mon2->n1));
	if(mon1->n2 != mon2->n2) return((mon1->n2 < mon2->n2));
	if(mon1->n3 != mon2->n3) return((mon1->n3 < mon2->n3));
	if(mon1->n4 != mon2->n4) return((mon1->n4 < mon2->n4));
	return(-1);
};

/* This function tests for divisibility of terms.	*/
int deelbaar(struct term *mon1, struct term *mon2)
{
	return(((mon1->n1 <= mon2->n1) && (mon1->n2 <= mon2->n2) && (mon1->n3 <= mon2->n3) && (mon1->n4 <= mon2->n4) && (valuation(mon1->c) <= valuation(mon2->c))));
};

/* Prints a polynomial. 				*/
void print_pol(struct polynomial f)
{
	struct term *fterm;

	fterm = f.leading;
	if(!fterm) {
		printf("0\n");
		return;
	};
	while(fterm) {
		if(!sc_is_zero(fterm->c)) {
			printmscalar(fterm->c);
			if(fterm->n1) printf("* x^%d ",fterm->n1);
			if(fterm->n2) printf("* y^%d ",fterm->n2);
			if(fterm->n3) printf("* z^%d ",fterm->n3);
			if(fterm->n4) printf("* w^%d ",fterm->n4);
			if(fterm->next) printf("+\n");
			else printf("\n");
		} else {
			printf("ZERO COEFFICIENT! FIXME!\n");
			exit(1);
		};
		fterm = fterm->next;
	};
	return;
};

/* This function adds polynomials of the same degree	*
 * into a new polynomial structure: 			*
 * The pols MUST have the same degree.			*
 * The pols may be zero.				*
 * The pols MUST NOT have terms which are zero.		*/
struct polynomial pol_add(struct polynomial f, struct polynomial g)
{
	mscalar c;
	int vergelijk;
	struct polynomial uit;
	struct term *fterm, *gterm;
	struct term **ptrterm;
	uit.leading = NULL;
	make_scalar(c);

	if(f.degree != g.degree) {
		printf("Can't add these!\n");
		exit(1);
	};
	uit.degree = f.degree;
	ptrterm = &(uit.leading); /* uit.leading will be set later.*/
	fterm=f.leading;
	gterm=g.leading;
	while (1) {
		if(!fterm) {
			*ptrterm = NULL;
			copy_tail(gterm,ptrterm);
			free_scalar(c);
			return(uit);
		};
		if(!gterm) {
			copy_tail(fterm,ptrterm);
			free_scalar(c);
			return(uit);
		};

		vergelijk=kleiner(fterm,gterm);
		if (vergelijk == GROTER) {
			make_term(ptrterm);
			copy_term(fterm,*ptrterm);
			ptrterm = &((*ptrterm)->next);
			fterm = fterm->next;
		} else if (vergelijk == KLEINER) {
			make_term(ptrterm);
			copy_term(gterm,*ptrterm);
			ptrterm = &((*ptrterm)->next);
			gterm = gterm->next;
		} else {
			/* vergelijk == GELIJK */
			sc_add(fterm->c,gterm->c,c);
			if(sc_is_zero(c)) {
				fterm=fterm->next;
				gterm=gterm->next;
			} else {
				make_term(ptrterm);
				sc_copy(c,(*ptrterm)->c);
				(*ptrterm)->n1 = fterm->n1;
				(*ptrterm)->n2 = fterm->n2;
				(*ptrterm)->n3 = fterm->n3;
				(*ptrterm)->n4 = fterm->n4;
				ptrterm = &((*ptrterm)->next);
				fterm = fterm->next;
				gterm = gterm->next;
			};
		};
	};
	printf("Fall through! Cannot happen.");
	exit(1);
};


/* Same as above but replace f by (f+g).		*
 * Here g may have terms that are zero but not f.	*/
void rep_pol_add(struct polynomial *f, struct polynomial g)
{
	int vergelijk;
	struct term *fterm, *gterm;
	struct term **ptrterm;

	if(f->degree != g.degree) {
		printf("Can't add these!\n");
		printf("Degree f is %d and degree g is %d\n",
				f->degree,g.degree);
		print_pol(*f);
		printf("\n");
		print_pol(g);
		exit(1);
	};
	ptrterm = &(f->leading);
	fterm=f->leading;
	*ptrterm = NULL;
	gterm=g.leading;
	while (1) {
		if(!fterm) {
			/* if(!gterm) return; */
			/* copy_tail with check for zero */
			while(gterm)  {
				if(!sc_is_zero(gterm->c)) {
					make_term(ptrterm);
					copy_term(gterm,*ptrterm);
					ptrterm = &((*ptrterm)->next);
				};
				gterm = gterm->next;
			};
			return;
		};

		if(!gterm) {
			*ptrterm = fterm;
			return;
		};

		vergelijk=kleiner(fterm,gterm);
		if (vergelijk == GROTER) {
			*ptrterm = fterm;
			ptrterm = &(fterm->next);
			fterm = fterm->next;
			*ptrterm = NULL;
		} else if (vergelijk == KLEINER) {
			/* Check for zero in g. */
			if(!sc_is_zero(gterm->c)) {
				make_term(ptrterm);
				copy_term(gterm,*ptrterm);
				ptrterm = &((*ptrterm)->next);
			};
			gterm = gterm->next;
		} else {
			/* vergelijk == GELIJK */
			sc_add_replace(gterm->c,fterm->c);
			if(sc_is_zero(fterm->c)) {
				/* Here we use *ptrterm as temp	*
				 * storage. A little ugly.	*/
				*ptrterm = fterm->next;
				free_term(fterm);
				fterm = *ptrterm;
				/* Now we put *ptrterm back to NULL. */
				*ptrterm = NULL;
				gterm=gterm->next;
			} else {
				*ptrterm = fterm;
				ptrterm = &(fterm->next);
				fterm = fterm->next;
				*ptrterm = NULL;
				gterm = gterm->next;
			};
		};
	};
	printf("Fall through! Cannot happen.");
	exit(1);
};

/* Replaces f by the product (possibly null).		*/
void times_int(int c, struct polynomial *f)
{
	struct term *tmp;
	struct term **ptrterm;

	ptrterm = &(f->leading);
	while(*ptrterm) {
		sc_imult_replace(c,(*ptrterm)->c);
		if(sc_is_zero((*ptrterm)->c)) {
			tmp = *ptrterm;
			*ptrterm = (*ptrterm)->next;
			free_term(tmp);
		} else {
			ptrterm = &((*ptrterm)->next);
		};
	};
	return;
};


/* Replaces f by the product (possibly null).		*/
void times_scalar(mscalar c, struct polynomial *f)
{
	struct term *tmp;
	struct term **ptrterm;

	ptrterm = &(f->leading);
	while(*ptrterm) {
		sc_mult_replace(c,(*ptrterm)->c);
		if(sc_is_zero((*ptrterm)->c)) {
			tmp = *ptrterm;
			*ptrterm = (*ptrterm)->next;
			free_term(tmp);
		} else {
			ptrterm = &((*ptrterm)->next);
		};
	};
	return;
};


/* This function assumes the data structures f and g	*
 * have the same length, and f is initialized, but	*
 * g need not be. The result t*f is put into g, but	*
 * it may not be a polynomial since some terms may be	*
 * zero.						*/
void times_term(struct term t, struct polynomial f, struct polynomial *g)
{
	struct term *fterm, *gterm;

	g->degree = f.degree + d1*t.n1 + d2*t.n2 + d3*t.n3 + d4*t.n4;
	gterm = g->leading;
	fterm = f.leading;
	while(fterm) {
		sc_mult(t.c, fterm->c, gterm->c);
		gterm->n1 = t.n1 + fterm->n1;
		gterm->n2 = t.n2 + fterm->n2;
		gterm->n3 = t.n3 + fterm->n3;
		gterm->n4 = t.n4 + fterm->n4;
		fterm = fterm->next;
		gterm = gterm->next;
	};
	return;
};


/* Same as above but it creates and outputs the product. 	*
 * Same caveats as above.					*/
struct polynomial
make_times_term(struct term t, struct polynomial f)
{
	struct term *fterm;
	struct term **ptrterm;
	struct polynomial uit;
	uit.leading = NULL;

	uit.degree = f.degree + d1*t.n1 + d2*t.n2 + d3*t.n3 + d4*t.n4;
	ptrterm = &(uit.leading);
	fterm = f.leading;
	while(fterm) {
		make_term(ptrterm);
		sc_mult(t.c, fterm->c, (*ptrterm)->c);
		(*ptrterm)->n1 = t.n1 + fterm->n1;
		(*ptrterm)->n2 = t.n2 + fterm->n2;
		(*ptrterm)->n3 = t.n3 + fterm->n3;
		(*ptrterm)->n4 = t.n4 + fterm->n4;
		fterm = fterm->next;
		ptrterm = &((*ptrterm)->next);
	};
	return(uit);
};

struct polynomial pol_mult(struct polynomial f, struct polynomial g)
{
	struct polynomial uit, tmppol;
	struct term *fterm;
	tmppol.leading = NULL;
	uit.leading = NULL;

	if((!f.leading) || (!g.leading)) {
		uit.degree = f.degree + g.degree;
		return(uit);
	};
	make_tail(g.leading, &(tmppol.leading));
	fterm = f.leading;
	while((!uit.leading) && (fterm)) {
		uit = make_times_term(*fterm, g);
		clean_pol(&uit);
		fterm = fterm->next;
	};
	while(fterm) {
		times_term(*fterm, g, &tmppol);
		rep_pol_add(&uit,tmppol);
		fterm = fterm->next;
	};
	free_tail(tmppol.leading);
	return(uit);
};


#ifdef DRAADJES

/* This code is a quick implementation of using threads to multiply	*
 * polynomials. Not used currently.					*/

/* Small function to execute thread. */
void *pol_mult_usu(struct polynomial **ab)
{
	*ab[2] = pol_mult(*ab[0],*ab[1]);
	pthread_exit(NULL);
}

/* Compute the product of two polynomials by doing a pair of threads. */
struct polynomial pol_mult_new(struct polynomial f, struct polynomial g)
{
	int nf=0;
	struct polynomial f1,f2,uit;
	struct polynomial **fg1,**fg2;
	struct term *fterm;
	pthread_t th1, th2;
	uit.leading = NULL;

	if((!f.leading) || (!g.leading)) {
		uit.degree = f.degree + g.degree;
		return(uit);
	};

	/* First count the number of terms in f. */
	fterm=f.leading;
	while(fterm) {
		nf++;
		fterm = fterm->next;
	};

	/* Two threads, so split f in two equal length pieces. */
	nf = nf/2;
	f1.degree = f.degree;
	f1.leading = f.leading;
	f2.degree = f.degree;
	fterm = f.leading;
	while(nf>0) {
		nf--;
		fterm = fterm->next;
	};
	f2.leading = fterm->next;
	fterm->next = NULL;

	/* Prepare the data structures. */
	fg1 =(struct polynomial **) malloc(3*sizeof(struct polynomial *));
	fg1[0] = &f1;
	fg1[1] = &g;
	fg1[2] = &uit;
	fg1[2]->leading = NULL;
	fg2 =(struct polynomial **) malloc(3*sizeof(struct polynomial *));
	fg2[0] = &f2;
	fg2[1] = &g;
	fg2[2] = (struct polynomial *) malloc(sizeof(struct polynomial));
	fg2[2]->leading = NULL;

	/* Create threads. */
	pthread_create(&th1,NULL,(void *)pol_mult_usu,fg1);
	pthread_create(&th2,NULL,(void *)pol_mult_usu,fg2);

	/* Wait till done. */
	pthread_join(th1,NULL);
	pthread_join(th2,NULL);

	/* Add results back together. */
	rep_pol_add(fg1[2],*fg2[2]);

	/* Free second one. */
	free_tail((fg2[2])->leading);
	free(fg2[2]);
	free(fg1);
	free(fg2);
	
	/* Put f back together in one piece. */
	fterm->next = f2.leading;

	return(uit);
};
#endif
