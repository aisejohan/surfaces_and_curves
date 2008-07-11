/*
 *	scalar.h
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

#define DIV4(a,b)	div4_inline_asm_C(a,b)
#define ADD4(a,b,c)	add4_inline_asm_C(a,b,c)
#define MUL4(a,b,c)	mul4_inline_asm_C(a,b,c)
#define NEG4(a)		neg4_inline_asm_C(a)
#define VAL4(a)		val4(a)
#define	SET4(a,b)	{a[0]=b[0];a[1]=b[1];a[2]=b[2];a[3]=b[3];}
#define	EQUAL4(a,b)	((a[0]==b[0])&&(a[1]==b[1])&&(a[2]==b[2])&&(a[3]==b[3]))
#define ISZERO4(a)	((a[0]==0)&&(a[1]==0)&&(a[2]==0)&&(a[3]==0))
#define INV4(a)		inv4_C(a)
#define	I_TO_4(a,k)	i_to_4_C(a,k)
#define	PRINT4(a)	printf("%lu + %lu*2^64 + %lu*2^128 + %lu*2^192",\
							a[0],a[1],a[2],a[3])

static inline int equal4_C(unsigned long *A, unsigned long *B)
{
	return ((A[0] == B[0]) && (A[1] == B[1]) && (A[2] == B[2])
							&& (A[3] == B[3]));
}

static inline int is_zero4_C(unsigned long *A)
{
	return ((A[0] == 0) && (A[1] == 0) && (A[2] == 0) && (A[3] == 0));
}

static inline void set4_C(unsigned long *A, unsigned long *B)
{
	A[0] = B[0];
	A[1] = B[1];
	A[2] = B[2];
	A[3] = B[3];
}

static inline void chop_C(unsigned long *couple, unsigned long a)
{
	couple[0] = (a << 32) >> 32;
	couple[1] = a >> 32;
}

static inline void add_C(unsigned long *triple, unsigned long a, unsigned long b)
{
	unsigned long temp;
	unsigned long ac[2], bc[2], tempc[2];

	chop_C(ac, a);
	chop_C(bc, b);
	temp = ac[0] + bc[0];
	chop_C(tempc, temp);
	triple[0] = tempc[0];
	temp = tempc[1] + ac[1] + bc[1];
	chop_C(tempc, temp);
	triple[1] = tempc[0];
	triple[2] = tempc[1];
}

static inline void add4_C(unsigned long *S, unsigned long *A, unsigned long *B)
{
	unsigned long couple[2];
	unsigned long triple[3];
	unsigned long s[8];

	add_C(triple, A[0], B[0]);
	s[0] = triple[0];
	s[1] = triple[1];
	s[2] = triple[2];
	add_C(triple, A[1], B[1]);
	s[2] += triple[0];
	s[3] = triple[1];
	s[4] = triple[2];
	add_C(triple, A[2], B[2]);
	s[4] += triple[0];
	s[5] = triple[1];
	s[6] = triple[2];
	add_C(triple, A[3], B[3]);
	s[6] += triple[0];
	s[7] = triple[1];
	chop_C(couple, s[2]);
	s[2] = couple[0];
	s[3] += couple[1];
	chop_C(couple, s[3]);
	s[3] = couple[0];
	s[4] += couple[1];
	chop_C(couple, s[4]);
	s[4] = couple[0];
	s[5] += couple[1];
	chop_C(couple, s[5]);
	s[5] = couple[0];
	s[6] += couple[1];
	chop_C(couple, s[6]);
	s[6] = couple[0];
	s[7] += couple[1];
	S[0] = s[0] + (s[1] << 32);
	S[1] = s[2] + (s[3] << 32);
	S[2] = s[4] + (s[5] << 32);
	S[3] = s[6] + (s[7] << 32);
}

static inline void mul4_C(unsigned long *S, unsigned long *A, unsigned long *B)
{
	unsigned long t1, t2;
	unsigned long ac[2], bc[2], sc[2];
	unsigned long triple[3];
	unsigned long s[8];

	chop_C(ac, A[0]);
	chop_C(bc, B[0]);
	t1 = ac[0] * bc[0];
	chop_C(sc, t1);
	s[0] = sc[0];
	s[1] = sc[1];
	t1 = ac[0] * bc[1];
	t2 = ac[1] * bc[0];
	add_C(triple, t1, t2);
	s[1] += triple[0];
	s[2] = triple[1];
	s[3] = triple[2];
	t1 = ac[1] * bc[1];
	chop_C(sc, t1);
	s[2] += sc[0];
	s[3] += sc[1];

/*	chop_C(ac, A[0]); */
	chop_C(bc, B[1]);
	t1 = ac[0] * bc[0];
	chop_C(sc, t1);
	s[2] += sc[0];
	s[3] += sc[1];
	t1 = ac[0] * bc[1];
	t2 = ac[1] * bc[0];
	add_C(triple, t1, t2);
	s[3] += triple[0];
	s[4] = triple[1];
	s[5] = triple[2];
	t1 = ac[1] * bc[1];
	chop_C(sc, t1);
	s[4] += sc[0];
	s[5] += sc[1];

/*	chop_C(ac, A[0]); */
	chop_C(bc, B[2]);
	t1 = ac[0] * bc[0];
	chop_C(sc, t1);
	s[4] += sc[0];
	s[5] += sc[1];
	t1 = ac[0] * bc[1];
	t2 = ac[1] * bc[0];
	add_C(triple, t1, t2);
	s[5] += triple[0];
	s[6] = triple[1];
	s[7] = triple[2];
	t1 = ac[1] * bc[1];
	chop_C(sc, t1);
	s[6] += sc[0];
	s[7] += sc[1];

/*	chop_C(ac, A[0]); */
/*	chop_C(bc, B[3]); */
	t1 = A[0] * B[3];
	chop_C(sc, t1);
	s[6] += sc[0];
	s[7] += sc[1];

	chop_C(ac, A[1]);
	chop_C(bc, B[0]);
	t1 = ac[0] * bc[0];
	chop_C(sc, t1);
	s[2] += sc[0];
	s[3] += sc[1];
	t1 = ac[0] * bc[1];
	t2 = ac[1] * bc[0];
	add_C(triple, t1, t2);
	s[3] += triple[0];
	s[4] += triple[1];
	s[5] += triple[2];
	t1 = ac[1] * bc[1];
	chop_C(sc, t1);
	s[4] += sc[0];
	s[5] += sc[1];

/*	chop_C(ac, A[1]); */
	chop_C(bc, B[1]);
	t1 = ac[0] * bc[0];
	chop_C(sc, t1);
	s[4] += sc[0];
	s[5] += sc[1];
	t1 = ac[0] * bc[1];
	t2 = ac[1] * bc[0];
	add_C(triple, t1, t2);
	s[5] += triple[0];
	s[6] += triple[1];
	s[7] += triple[2];
	t1 = ac[1] * bc[1];
	chop_C(sc, t1);
	s[6] += sc[0];
	s[7] += sc[1];

/*	chop_C(ac, A[1]); */
/*	chop_C(bc, B[2]); */
	t1 = A[1] * B[2];
	chop_C(sc, t1);
	s[6] += sc[0];
	s[7] += sc[1];

	chop_C(ac, A[2]);
	chop_C(bc, B[0]);
	t1 = ac[0] * bc[0];
	chop_C(sc, t1);
	s[4] += sc[0];
	s[5] += sc[1];
	t1 = ac[0] * bc[1];
	t2 = ac[1] * bc[0];
	add_C(triple, t1, t2);
	s[5] += triple[0];
	s[6] += triple[1];
	s[7] += triple[2];
	t1 = ac[1] * bc[1];
	chop_C(sc, t1);
	s[6] += sc[0];
	s[7] += sc[1];

/*	chop_C(ac, A[2]); */
/*	chop_C(bc, B[1]); */
	t1 = A[2] * B[1];
	chop_C(sc, t1);
	s[6] += sc[0];
	s[7] += sc[1];

/*	chop_C(ac, A[3]); */
/*	chop_C(bc, B[0]); */
	t1 = A[3] * B[0];
	chop_C(sc, t1);
	s[6] += sc[0];
	s[7] += sc[1];

	chop_C(sc, s[0]);
	s[0] = sc[0];
	s[1] += sc[1];
	chop_C(sc, s[1]);
	s[1] = sc[0];
	s[2] += sc[1];
	chop_C(sc, s[2]);
	s[2] = sc[0];
	s[3] += sc[1];
	chop_C(sc, s[3]);
	s[3] = sc[0];
	s[4] += sc[1];
	chop_C(sc, s[4]);
	s[4] = sc[0];
	s[5] += sc[1];
	chop_C(sc, s[5]);
	s[5] = sc[0];
	s[6] += sc[1];
	chop_C(sc, s[6]);
	s[6] = sc[0];
	s[7] += sc[1];
	chop_C(sc, s[7]);
	s[7] = sc[0];

	S[0] = s[0] + (s[1] << 32);
	S[1] = s[2] + (s[3] << 32);
	S[2] = s[4] + (s[5] << 32);
	S[3] = s[6] + (s[7] << 32);
}

static inline void neg4_C(unsigned long *A)
{
	unsigned long B[4], C[4];

	C[0] = A[0] ^ ((unsigned long) -1);
	C[1] = A[1] ^ ((unsigned long) -1);
	C[2] = A[2] ^ ((unsigned long) -1);
	C[3] = A[3] ^ ((unsigned long) -1);
	B[0] = 1;
	B[1] = 0;
	B[2] = 0;
	B[3] = 0;
	add4_C(A, B, C);
}

static inline void i_to_4_C(unsigned long *A, int k)
{
	A[1] = 0;
	A[2] = 0;
	A[3] = 0;
	if (k >= 0) {
		A[0] = k;
	} else {
		A[0] = -k;
		neg4_C(A);
	}
}

static inline void div4_C(unsigned long *A, unsigned int k)
{
	if (k == 0) return;
	A[0] = (A[0] >> k) | (A[1] << (64 - k));
	A[1] = (A[1] >> k) | (A[2] << (64 - k));
	A[2] = (A[2] >> k) | (A[3] << (64 - k));
	A[3] = A[3] >> k;
}


static inline int val4_C(unsigned long *A)
{
	unsigned long a, b;

        if (A[0] == 0) {
		if (A[1] == 0) {
			if (A[2] == 0) {
				if (A[3] == 0) {
					return(256);
				} else {
					a = 192;
					b = A[3];
				}
			} else {
				a = 128;
				b = A[2];
			}
		} else {
			a = 64;
			b = A[1];
		}
	} else {
		a = 0;
		b = A[0];
	}
        __asm__ ("bsfq %1,%0" : "=r" (b) : "r" (b));
        return((int) a + b);
}

static inline void rand4_C(unsigned long *A)
{
	A[0] = (unsigned long) rand() + ((unsigned long) rand() << 1) + ((unsigned long) rand() << 32) + ((unsigned long) rand() << 33);
	A[1] = (unsigned long) rand() + ((unsigned long) rand() << 1) + ((unsigned long) rand() << 32) + ((unsigned long) rand() << 33);
	A[2] = (unsigned long) rand() + ((unsigned long) rand() << 1) + ((unsigned long) rand() << 32) + ((unsigned long) rand() << 33);
	A[3] = (unsigned long) rand() + ((unsigned long) rand() << 1) + ((unsigned long) rand() << 32) + ((unsigned long) rand() << 33);
}

static inline void inv4_C(unsigned long *A)
{
	int e;
	unsigned long B[4], C[4], D[4];

	SET4(C, A);
	A[0]--; /* No overflow as A is unit 1+2... */
	e = VAL4(A); /* C is 1+2^e x */
	NEG4(A);
	A[0]++; /* No overflow as A[0] is even. Now A = 1 - 2^e x */
	MUL4(D, C, A);
	SET4(C, D);
	e = 2*e;

	while (e < 256) {
		SET4(B, C); /* What is left over */
		B[0]--;
		NEG4(B);
		B[0]++; /* Now B = 1 - 2^{2e} y */ 
		MUL4(D, C, B);
		SET4(C, D); /* C is updated. */
		MUL4(D, A, B);
		SET4(A, D); /* A is updated. */
		e = 2*e;
	}
}
