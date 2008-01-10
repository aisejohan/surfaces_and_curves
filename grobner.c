/*
 *	grobner.c
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

#include <stdlib.h>
#include <stdio.h>

#include "data.h"
#include "scalar.h"
#include "pol.h"

/* Generic routines... */

/* This function tests for divisibility of terms.	*/
static int deelbaar(struct term *mon1, struct term *mon2)
{
	return(((mon1->n1 <= mon2->n1) &&
		(mon1->n2 <= mon2->n2) &&
		(mon1->n3 <= mon2->n3) &&
		(valuation(mon1->c) <= valuation(mon2->c))));
}


/****************************************************************
 * Really insane and screwed up. Removes all monomials it can.	*
 * The remainder ends up in pp.					*
 * The coefficients aa[i] will be returned, so that		*
 *								*
 *	(input pp) + sum_{i=0,ss-1} aa[i] * vh[i]		*
 *								*
 *		=						*
 *								*
 *	(output pp)						*
 *								*
 * **************************************************************/
#ifdef OLD_GROBNER
struct polynomial ** 
gen_division(struct polynomial *pp, unsigned int ss, struct polynomial **vh)
{
	struct polynomial tmp[ss];
	struct polynomial vh_rest[ss];
	struct polynomial **aa;
	struct polynomial *ppp;
	struct term *aaterm[ss];
	struct term **ptrterm;
	struct term *pppterm;
	struct term mon;
	unsigned int i, dividing;

	make_scalar(mon.c);
	ppp = NULL;
	make_pol(&ppp);
	aa = (struct polynomial **)malloc(ss*sizeof(struct polynomial *));
	if(!aa) {
		perror("Malloc failed!");
		exit(1);
	};
	for(i=0;i+1<=ss;i++) {
		tmp[i].leading = NULL;
		aaterm[i] = NULL;
		aa[i] = NULL;
		make_pol(&aa[i]);
		aa[i]->degree = (pp->degree > vh[i]->degree) ?
			(pp->degree - vh[i]->degree) : 0;
	};

	/* Copy pp into ppp. */
	ppp->degree = pp->degree;
	ppp->leading = pp->leading;
	/* Set pp equal to ``zero'' */
	pp->leading = NULL;
	ptrterm = &pp->leading;

	while (ppp->leading) {
		i = 0;
		dividing = 1;
		while (i+1 <= ss && dividing) {
			if (deelbaar(vh[i]->leading, ppp->leading)) {
				/* No sign in front of pppterm->c */
				sc_div(ppp->leading->c, vh[i]->leading->c,
					mon.c);
				/* Change sign mon.c */
				sc_negate(mon.c);
				mon.n1 = ppp->leading->n1 - vh[i]->leading->n1;
				mon.n2 = ppp->leading->n2 - vh[i]->leading->n2;
				mon.n3 = ppp->leading->n3 - vh[i]->leading->n3;

				pppterm = ppp->leading;
				ppp->leading = ppp->leading->next;
				free_term(pppterm);

				if (aaterm[i]) {
					times_term(mon, vh_rest[i], &(tmp[i]));
					make_term(&aaterm[i]->next);
					copy_term(&mon, aaterm[i]->next);
					aaterm[i] = aaterm[i]->next;
				} else {
					vh_rest[i].degree = vh[i]->degree;
					vh_rest[i].leading = 
						vh[i]->leading->next;
					tmp[i] = make_times_term(mon,
						vh_rest[i]);
					make_term(&aa[i]->leading);
					copy_term(&mon, aa[i]->leading);
					aaterm[i] = aa[i]->leading;
				};

				rep_pol_add(ppp, tmp[i]);

				dividing = 0;
			} else {
				i=i+1;
			};
		};
		/* dividing == 1 means that we cannot get rid of the leading
		 * term. So we put it back in pp. */
		if(dividing) {
			*ptrterm = ppp->leading;
			ptrterm = &((*ptrterm)->next);
			/* Move on to the next one. */
			ppp->leading = ppp->leading->next;
			/* Terminate pp. */
			*ptrterm = NULL;
		};
	};
	for(i=0;i+1<=ss;i++) {
		free_tail(tmp[i].leading);
	};
	free(ppp);
	free_scalar(mon.c);
	return(aa);
}
#endif

#ifdef MIXED_GROBNER
struct polynomial ** 
gen_division(struct polynomial *pp, unsigned int ss, struct polynomial **vh)
{
	struct polynomial tmp[ss];
	struct polynomial save_the_spot, uit;
	struct term test;
	struct polynomial vh_rest[ss];
	struct polynomial **aa;
	struct polynomial *ppp;
	struct term *aaterm[ss];
	struct term **ptrterm;
	struct term *pppterm;
	struct term mon;
	unsigned int i, dividing;

	make_scalar(mon.c);
	ppp = NULL;
	make_pol(&ppp);
	aa = (struct polynomial **)malloc(ss*sizeof(struct polynomial *));
	if(!aa) {
		perror("Malloc failed!");
		exit(1);
	};
	for(i=0;i+1<=ss;i++) {
		tmp[i].leading = NULL;
		aaterm[i] = NULL;
		aa[i] = NULL;
		make_pol(&aa[i]);
		aa[i]->degree = (pp->degree > vh[i]->degree) ?
			(pp->degree - vh[i]->degree) : 0;
	};

	/* Copy pp into ppp. */
	ppp->degree = pp->degree;
	ppp->leading = pp->leading;
	/* Set pp equal to ``zero'' */
	pp->leading = NULL;
	ptrterm = &pp->leading;

	while (ppp->leading) {
		i = 0;
		dividing = 1;
		while (i+1 <= ss && dividing) {
			if (deelbaar(vh[i]->leading, ppp->leading)) {
				/* No sign in front of pppterm->c */
				sc_div(ppp->leading->c, vh[i]->leading->c,
					mon.c);
				/* Change sign mon.c */
				sc_negate(mon.c);
				mon.n1 = ppp->leading->n1 - vh[i]->leading->n1;
				mon.n2 = ppp->leading->n2 - vh[i]->leading->n2;
				mon.n3 = ppp->leading->n3 - vh[i]->leading->n3;

				pppterm = ppp->leading;
				ppp->leading = ppp->leading->next;
				free_term(pppterm);

				if (aaterm[i]) {
					if (ss > 1) /* Old part. */
						times_term(mon, vh_rest[i],
							&(tmp[i]));
					make_term(&aaterm[i]->next);
					copy_term(&mon, aaterm[i]->next);
					aaterm[i] = aaterm[i]->next;
				} else {
					vh_rest[i].degree = vh[i]->degree;
					vh_rest[i].leading = 
						vh[i]->leading->next;
					if (ss > 1) /* Old part. */
						tmp[i] = make_times_term(mon,
							vh_rest[i]);
					make_term(&aa[i]->leading);
					copy_term(&mon, aa[i]->leading);
					aaterm[i] = aa[i]->leading;
				};

				if (ss == 1) {
					/* New part. This part seems to
					 * work better with LEX_ORDER
					 * than with REVLEX_ORDER. */
					save_the_spot.degree = aa[i]->degree;
					save_the_spot.leading = aaterm[i];
					test.n1 = mon.n1 +
						vh_rest[i].leading->n1;
					test.n2 = mon.n2 +
						vh_rest[i].leading->n2;
					test.n3 = mon.n3 +
						vh_rest[i].leading->n3;

					while ((ppp->leading) &&
					deelbaar(vh[i]->leading, ppp->leading)
					&& (GROTER == kleiner(ppp->leading,
					&test))) {
						/* No sign in front of
						 * pppterm->c */
						sc_div(ppp->leading->c,
						vh[i]->leading->c, mon.c);
						/* Change sign mon.c */
						sc_negate(mon.c);
						mon.n1 = ppp->leading->n1 -
							vh[i]->leading->n1;
						mon.n2 = ppp->leading->n2 -
							vh[i]->leading->n2;
						mon.n3 = ppp->leading->n3 -
							vh[i]->leading->n3;
						pppterm = ppp->leading;
						ppp->leading = 
							ppp->leading->next;
						free_term(pppterm);
						make_term(&aaterm[i]->next);
						copy_term(&mon,
							aaterm[i]->next);
						aaterm[i] = aaterm[i]->next;
					}
					uit = pol_mult(save_the_spot,
						vh_rest[i]);
					merge_add(ppp, uit);
					/* End new part. */
				} else {
					/* Old part. */
					rep_pol_add(ppp, tmp[i]);
					/* End old part. */
				}

				dividing = 0;
			} else {
				i=i+1;
			};
		};
		/* dividing == 1 means that we cannot get rid of the leading
		 * term. So we put it back in pp. */
		if(dividing) {
			*ptrterm = ppp->leading;
			ptrterm = &((*ptrterm)->next);
			/* Move on to the next one. */
			ppp->leading = ppp->leading->next;
			/* Terminate pp. */
			*ptrterm = NULL;
		};
	};
	if (ss > 1) {
		for(i=0;i+1<=ss;i++) {
			free_tail(tmp[i].leading);
		};
	}
	free(ppp);
	free_scalar(mon.c);
	return(aa);
}
#endif

#ifdef NEW_GROBNER
struct polynomial ** 
gen_division(struct polynomial *pp, unsigned int ss, struct polynomial **vh)
{
	struct polynomial save_the_spot, uit;
	struct term test;
	struct polynomial vh_rest[ss];
	struct polynomial **aa;
	struct polynomial *ppp;
	struct term *aaterm[ss];
	struct term **ptrterm;
	struct term *pppterm;
	struct term mon;
	unsigned int i, dividing;

	make_scalar(mon.c);
	ppp = NULL;
	make_pol(&ppp);
	aa = (struct polynomial **)malloc(ss*sizeof(struct polynomial *));
	if(!aa) {
		perror("Malloc failed!");
		exit(1);
	};
	for(i=0;i+1<=ss;i++) {
		aaterm[i] = NULL;
		aa[i] = NULL;
		make_pol(&aa[i]);
		aa[i]->degree = (pp->degree > vh[i]->degree) ?
			(pp->degree - vh[i]->degree) : 0;
	};

	/* Copy pp into ppp. */
	ppp->degree = pp->degree;
	ppp->leading = pp->leading;
	/* Set pp equal to ``zero'' */
	pp->leading = NULL;
	ptrterm = &pp->leading;

	while (ppp->leading) {
		i = 0;
		dividing = 1;
		while (i+1 <= ss && dividing) {
			if (deelbaar(vh[i]->leading, ppp->leading)) {
				/* No sign in front of pppterm->c */
				sc_div(ppp->leading->c, vh[i]->leading->c,
					mon.c);
				/* Change sign mon.c */
				sc_negate(mon.c);
				mon.n1 = ppp->leading->n1 - vh[i]->leading->n1;
				mon.n2 = ppp->leading->n2 - vh[i]->leading->n2;
				mon.n3 = ppp->leading->n3 - vh[i]->leading->n3;

				pppterm = ppp->leading;
				ppp->leading = ppp->leading->next;
				free_term(pppterm);

				if (aaterm[i]) {
					make_term(&aaterm[i]->next);
					copy_term(&mon, aaterm[i]->next);
					aaterm[i] = aaterm[i]->next;
				} else {
					vh_rest[i].degree = vh[i]->degree;
					vh_rest[i].leading = 
						vh[i]->leading->next;
					make_term(&aa[i]->leading);
					copy_term(&mon, aa[i]->leading);
					aaterm[i] = aa[i]->leading;
				};

				/* New part. This part seems to
				 * work better with LEX_ORDER
				 * than with REVLEX_ORDER. */
				save_the_spot.degree = aa[i]->degree;
				save_the_spot.leading = aaterm[i];

				if (vh_rest[i].leading) {
					test.n1 = mon.n1 + 
						vh_rest[i].leading->n1;
					test.n2 = mon.n2 +
						vh_rest[i].leading->n2;
					test.n3 = mon.n3 +
						vh_rest[i].leading->n3;
				}

				while ((ppp->leading) &&
				deelbaar(vh[i]->leading, ppp->leading)
				&& ((!vh_rest[i].leading) ||
				(GROTER == kleiner(ppp->leading, &test)))) {
					/* No sign in front of
					 * pppterm->c */
					sc_div(ppp->leading->c,
						vh[i]->leading->c, mon.c);
					/* Change sign mon.c */
					sc_negate(mon.c);
					mon.n1 = ppp->leading->n1 -
						vh[i]->leading->n1;
					mon.n2 = ppp->leading->n2 -
						vh[i]->leading->n2;
					mon.n3 = ppp->leading->n3 -
						vh[i]->leading->n3;
					pppterm = ppp->leading;
					ppp->leading = ppp->leading->next;
					free_term(pppterm);
					make_term(&aaterm[i]->next);
					copy_term(&mon, aaterm[i]->next);
					aaterm[i] = aaterm[i]->next;
				}
				uit = pol_mult(save_the_spot, vh_rest[i]);
				merge_add(ppp, uit);

				dividing = 0;
			} else {
				i=i+1;
			};
		};
		/* dividing == 1 means that we cannot get rid of the leading
		 * term. So we put it back in pp. */
		if(dividing) {
			*ptrterm = ppp->leading;
			ptrterm = &((*ptrterm)->next);
			/* Move on to the next one. */
			ppp->leading = ppp->leading->next;
			/* Terminate pp. */
			*ptrterm = NULL;
		};
	};
	free(ppp);
	free_scalar(mon.c);
	return(aa);
}
#endif

/* ppp does not get changed */
unsigned int
zero_on_division(struct polynomial ppp, unsigned int ss, struct polynomial **vh)
{
	struct polynomial tmp[ss];
	struct term mon;
	struct polynomial pp;
	unsigned int i, dividing;
	make_scalar(mon.c);
	pp.leading = NULL;
	for(i=0;i+1<=ss;i++) {
		tmp[i].leading = NULL;
	};

	pp.degree = ppp.degree;
	copy_tail(ppp.leading,&(pp.leading));
	while(pp.leading) {
		i = 0;
		dividing = 1;
		while((i+1<=ss) && dividing) {
			if(deelbaar(vh[i]->leading, pp.leading)) {
				/* No sign in front of ppterm->c */
				sc_div(pp.leading->c, vh[i]->leading->c,mon.c);
				/* Change sign mon.c */
				sc_negate(mon.c);
				mon.n1 = pp.leading->n1 - vh[i]->leading->n1;
				mon.n2 = pp.leading->n2 - vh[i]->leading->n2;
				mon.n3 = pp.leading->n3 - vh[i]->leading->n3;
				if(tmp[i].leading) {
					times_term(mon, *(vh[i]), &(tmp[i]));
				} else {
					tmp[i] = make_times_term(mon,
						*(vh[i]));
				};
				rep_pol_add(&pp, tmp[i]);
				dividing = 0;
			} else {
				i=i+1;
			};
		};
		if(dividing) {
			for(i=0;i+1<=ss;i++) {
				free_tail(tmp[i].leading);
			};
			free_tail(pp.leading);
			free_scalar(mon.c);
			return(0);
		};
	};
	for(i=0;i+1<=ss;i++) {
		free_tail(tmp[i].leading);
	};
	free_tail(pp.leading);
	free_scalar(mon.c);
	return(1);
}
