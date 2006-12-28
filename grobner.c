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
struct polynomial ** 
gen_division(struct polynomial *pp, unsigned int ss, struct polynomial **vh)
{
	struct polynomial tmp[ss];
	struct polynomial **aa;
	struct polynomial *ppp;
	struct term *aaterm[ss];
	struct term **ptrterm;
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
				mon.n4 = ppp->leading->n4 - vh[i]->leading->n4;

				if(tmp[i].leading) {
					times_term(mon, *(vh[i]), &(tmp[i]));
					make_term(&aaterm[i]->next);
					copy_term(&mon, aaterm[i]->next);
					aaterm[i] = aaterm[i]->next;
				} else {
					tmp[i] = make_times_term(mon,*(vh[i]));
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
};

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
				mon.n4 = pp.leading->n4 - vh[i]->leading->n4;
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
};
