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

unsigned int ivaluation(int a)
{
	unsigned int uit = 0;
	while (!(a % p)) {
		a = a / p;
		uit++;
	}
	return(uit);
}

void set_seed(unsigned int zaadje)
{
	int fd, uit;
	unsigned int willekeurig;
	if (zaadje) {
		srand(zaadje);
		return;
	}
	fd = open("/dev/urandom", O_RDONLY);
	if (fd < 0) {
		printf("Unable to open /dev/urandom. So seed=666.\n");
		willekeurig = 666;
	} else {
		uit = read(fd, &willekeurig, sizeof(willekeurig));
		if (uit <= 0) {
			printf("Failure reading /dev/urandom. So Seed=666.\n");
			willekeurig = 666;
		}
	}
	srand(willekeurig);
	uit = close(fd);
	return;
}

unsigned int count_sum(unsigned int degree)
{
	unsigned int count, a1, a2, a3;
	count = 0;
	for(a1=0;(d1*a1 <= degree);a1++) {
		for(a2=0;(d1*a1+d2*a2 <= degree);a2++) {
			for(a3=0;(d1*a1+d2*a2+d3*a3 <= degree);a3++) {
				if((degree - (a1*d1+a2*d2+a3*d3)) % d4 == 0) {
					count++;
				};
			};
		};
	};
	return(count);
}

void print_sum(unsigned int degree)
{
	unsigned int a1, a2, a3;
	for (a1 = 0; (d1*a1 <= degree); a1++) {
	  for (a2 = 0; (d1*a1 + d2*a2 <= degree); a2++) {
	    for (a3 = 0; (d1*a1 + d2*a2 + d3*a3 <= degree); a3++) {
	      if ((degree - (a1*d1 + a2*d2 + a3*d3)) % d4 == 0) {
		printf("[%d, %d, %d, %d]\n", a1, a2, a3,
				(degree - (a1*d1 + a2*d2 + a3*d3))/d4);
	      }
	    }
	  }
	}
	return;
}

/* Makes a random polynomial of degree degree.		*
 * The result may be the zero polynomial!		*/
struct polynomial make_random(unsigned int degree, int print)
{
	unsigned int a1, a2, a3, a4;
	int c;
	struct polynomial uit;
	struct term *uitterm;
	struct term **ptrterm;
	uitterm = NULL;
	uit.degree = degree;
	uit.leading = NULL;

	if (!count_sum(degree)) {
		printf("No monomials of degree %d! Stop.\n", degree);
		exit(1);
	}
#ifdef INPUT_F
 #ifndef OUTPUT_LIST
	printf("\n");
	printf("Please input coefficients below.\n");
 #endif
#endif
	for (a1 = 0; (d1*a1 <= degree); a1++) {
	  for (a2 = 0; (d1*a1 + d2*a2 <= degree); a2++) {
	    for (a3 = 0; (d1*a1 + d2*a2 + d3*a3 <= degree); a3++) {
	      if ((degree - (a1*d1 + a2*d2 + a3*d3)) % d4 == 0) {
		a4 = (degree - (a1*d1 + a2*d2 + a3*d3))/d4;
#ifdef INPUT_F
		/* Dummy input at first. */
		c = 1;
#else
		/* Stupid lift. */
		c = rand() % p;
		/* Change for compatibility with previous version.	*
		 * This does not make any difference to the pol mod p.	*/
		if (c < -(p - 1)/2) c += p; /* OK, this never happens. */
		if (c > (p - 1)/2) c -= p; /* This does happen. */
#endif
		/* Create the new term to be put in. */
		make_term(&uitterm);
		uitterm->n1 = a1;
		uitterm->n2 = a2;
		uitterm->n3 = a3;
		uitterm->n4 = a4;
		ito_sc(c,uitterm->c);
		ptrterm = &(uit.leading);
		while ((*ptrterm) && (kleiner(uitterm, *ptrterm) == KLEINER)) {
			ptrterm = &((*ptrterm)->next);
		}
		uitterm->next = *ptrterm;
		*ptrterm = uitterm;
		uitterm=NULL;
	      }
	    }
	  }
	}
	uitterm = uit.leading;
	while (uitterm) {
		a1 = uitterm->n1;
		a2 = uitterm->n2;
		a3 = uitterm->n3;
		a4 = uitterm->n4;
		c = 0;
		printf("Coefficient of   ");
		if (a1) {
			printf("x^%d",a1);
			c++;
		}
		if ((a1) && (a2+a3+a4)) {
			printf(" * ");
			c++;
		}
		if (a2) {
			printf("y^%d", a2);
			c++;
		}
		if ((a2) && (a3 + a4)) {
			printf(" * ");
			c++;
		}
		if (a3) {
			printf("z^%d", a3);
			c++;
		}
		if ((a3) && (a4)) {
			printf(" * ");
			c++;
		}
		if (a4) {
			printf("w^%d", a4);
			c++;
		}
		while (8 - c) {
			printf("   ");
			c++;
		}
		printf("= ");
#ifndef INPUT_F
		if (print) printmscalar(uitterm->c);
 		printf("\n");
#else
 #ifdef OUTPUT_LIST
 		printf("\n");
 #else
		scanf("%d", &c);
 #endif
		ito_sc(c, uitterm->c);
#endif
		uitterm = uitterm->next;
	}
	clean_pol(&uit);
#ifdef OUTPUT_LIST
	exit(0);
#endif
	return(uit);
}

/* Counts the number of terms and verifies f is a polynomial.	*
 * Allows f to be the zero polynomial, warns if degree not 0.	*/
unsigned int number_terms(struct polynomial f)
{
	int count = 0;
	struct term *fterm;

	if(!f.leading) {
		if(!f.degree) {
			printf("Zero but not degree 0.");
		}
		return(0);
	}
	fterm = f.leading;
	while(fterm) {
		count++;
		fterm = fterm->next;
	}
	return(count);
}

/* Frobenius operation by raising varaibles tot the	*
 * pth power.						*/
struct polynomial frobenius(struct polynomial f)
{
	struct polynomial uit;
	struct term *fterm;
	struct term **ptrterm;
	uit.leading = NULL;
	ptrterm = NULL;

	uit.degree = p * f.degree;
	ptrterm = &(uit.leading);
	fterm = f.leading;
	while (fterm) {
		make_term(ptrterm);
		sc_copy(fterm->c, (*ptrterm)->c);
		(*ptrterm)->n1 = p*fterm->n1;
		(*ptrterm)->n2 = p*fterm->n2;
		(*ptrterm)->n3 = p*fterm->n3;
		(*ptrterm)->n4 = p*fterm->n4;
		ptrterm = &((*ptrterm)->next);
		fterm = fterm->next;
	}
	return(uit);
}

/* Ouputs the derivative. 					*
 * The result is nonsense if the degree of f is too low.	*/
struct polynomial deriv(struct polynomial f, unsigned int i)
{
	mscalar c;
	struct polynomial uit;
	struct term *fterm;
	struct term **ptrterm;
	uit.leading = NULL;
	ptrterm = NULL;
	make_scalar(c);

	fterm = f.leading;
	switch (i) {
		case 1: 
		uit.degree = (f.degree > d1) ? (f.degree - d1) : 0;
		ptrterm = &(uit.leading);
		while (fterm) {
			sc_imult(fterm->n1, fterm->c,c);
			if (!sc_is_zero(c)) {
				make_term(ptrterm);
				sc_copy(c, (*ptrterm)->c);
				(*ptrterm)->n1 = fterm->n1 - 1;
				(*ptrterm)->n2 = fterm->n2;
				(*ptrterm)->n3 = fterm->n3;
				(*ptrterm)->n4 = fterm->n4;
				ptrterm = &((*ptrterm)->next);
			}
			fterm = fterm->next;
		}
		free_scalar(c);
		return(uit);

		case 2: 
		uit.degree = (f.degree > d2) ? (f.degree - d2) : 0;
		ptrterm = &(uit.leading);
		while (fterm) {
			sc_imult(fterm->n2, fterm->c, c);
			if (!sc_is_zero(c)) {
				make_term(ptrterm);
				sc_copy(c, (*ptrterm)->c);
				(*ptrterm)->n1 = fterm->n1;
				(*ptrterm)->n2 = fterm->n2 - 1;
				(*ptrterm)->n3 = fterm->n3;
				(*ptrterm)->n4 = fterm->n4;
				ptrterm = &((*ptrterm)->next);
			}
			fterm = fterm->next;
		}
		free_scalar(c);
		return(uit);

		case 3:
		uit.degree = (f.degree > d3) ? (f.degree - d3) : 0;
		ptrterm = &(uit.leading);
		while (fterm) {
			sc_imult(fterm->n3, fterm->c, c);
			if (!sc_is_zero(c)) {
				make_term(ptrterm);
				sc_copy(c, (*ptrterm)->c);
				(*ptrterm)->n1 = fterm->n1;
				(*ptrterm)->n2 = fterm->n2;
				(*ptrterm)->n3 = fterm->n3 - 1;
				(*ptrterm)->n4 = fterm->n4;
				ptrterm = &((*ptrterm)->next);
			}
			fterm = fterm->next;
		}
		free_scalar(c);
		return(uit);
		
		case 4:
		uit.degree = (f.degree > d4) ? (f.degree - d4) : 0;
		ptrterm = &(uit.leading);
		while (fterm) {
			sc_imult(fterm->n4, fterm->c, c);
			if (!sc_is_zero(c)) {
				make_term(ptrterm);
				sc_copy(c, (*ptrterm)->c);
				(*ptrterm)->n1 = fterm->n1;
				(*ptrterm)->n2 = fterm->n2;
				(*ptrterm)->n3 = fterm->n3;
				(*ptrterm)->n4 = fterm->n4 - 1;
				ptrterm = &((*ptrterm)->next);
			}
			fterm = fterm->next;
		}
		free_scalar(c);
		return(uit);
		
		default:
		printf("Wrong again honey!");
		exit(1);
	}
	exit(1);
}

/* Same as above. Replaces f by the derivative. 		*
 * The result is nonsense if the degree of f is too low.	*/
void rep_deriv(struct polynomial *f, unsigned int i)
{
	struct term *fterm;
	struct term **ptrterm;

	fterm = f->leading;
	ptrterm = &(f->leading);
	switch (i) {
		case 1: 
		f->degree = (f->degree > d1) ? (f->degree - d1) : 0;
		while (fterm) {
			sc_imult_replace(fterm->n1, fterm->c);
			if (!sc_is_zero(fterm->c)) {
				fterm->n1 = fterm->n1 - 1;
				*ptrterm = fterm;
				ptrterm = &(fterm->next);
				fterm = fterm->next;
			} else {
				*ptrterm = fterm->next;
				free_term(fterm);
				fterm = *ptrterm;
			}
		}
		return;

		case 2: 
		f->degree = (f->degree > d2) ? (f->degree - d2) : 0;
		while (fterm) {
			sc_imult_replace(fterm->n2, fterm->c);
			if (!sc_is_zero(fterm->c)) {
				fterm->n2 = fterm->n2 - 1;
				*ptrterm = fterm;
				ptrterm = &(fterm->next);
				fterm = fterm->next;
			} else {
				*ptrterm = fterm->next;
				free_term(fterm);
				fterm = *ptrterm;
			}
		}
		return;

		case 3: 
		f->degree = (f->degree > d3) ? (f->degree - d3) : 0;
		while (fterm) {
			sc_imult_replace(fterm->n3, fterm->c);
			if (!sc_is_zero(fterm->c)) {
				fterm->n3 = fterm->n3 - 1;
				*ptrterm = fterm;
				ptrterm = &(fterm->next);
				fterm = fterm->next;
			} else {
				*ptrterm = fterm->next;
				free_term(fterm);
				fterm = *ptrterm;
			}
		}
		return;

		case 4:
		f->degree = (f->degree > d4) ? (f->degree - d4) : 0;
		while (fterm) {
			sc_imult_replace(fterm->n4, fterm->c);
			if (!sc_is_zero(fterm->c)) {
				fterm->n4 = fterm->n4 - 1;
				*ptrterm = fterm;
				ptrterm = &(fterm->next);
				fterm = fterm->next;
			} else {
				*ptrterm = fterm->next;
				free_term(fterm);
				fterm = *ptrterm;
			}
		}
		return;

		default:
		printf("Wrong again honey!");
		exit(1);
	};
	exit(1);
}
