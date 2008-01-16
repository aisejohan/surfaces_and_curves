#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

/* really dumb. */
#define MAX	200
#define CUTOFF	100

int storage=0;
int **lijst;
int len;
int pairs[6][4] = {{0,1,2,3},{0,2,1,3},{0,3,1,2},{1,2,0,3},{1,3,0,2},{2,3,0,1}};

void my_alloc(void )
{
	int i;

	if (!storage) {
		lijst = (int **) malloc(100*sizeof(int *));
		storage = 100;
		for(i=0; i+1 <= 100;i++) {
			lijst[i] = (int *) malloc(3*sizeof(int));
		}
		return;
	}
	lijst = (int **) realloc(lijst, 2*storage*sizeof(int *));
	for(i=storage; i+1 <= 2*storage;i++) {
		lijst[i] = (int *) malloc(3*sizeof(int));
	}
	storage = 2*storage;
	return;
}
		

/* Returns the number of monomials in degree degree, and sets lijst
 * equal to the list of them.*/
int count_sum(int d1, int d2, int d3, int degree, int save)
{
	int len,a1,a2;
	
	len=0;
	a1 = 0;
	while (d1*a1 <= degree) {
		a2 = 0;
		while (d1*a1+d2*a2 <= degree) {
			if((degree - (a1*d1+a2*d2)) % d3 == 0) {
				if (save) {
					if (len == storage) my_alloc();
					lijst[len][0] = a1;
					lijst[len][1] = a2;
					lijst[len][2] = (degree - 
						(a1*d1+a2*d2)) / d3;
					}
					len++;
			}
			a2++;
		}
		a1++;
	}
	return(len);
}

/* This prints out all the multiindices whose associated monomials have degree degree. */
void print_sum(int len)
{
	int i;

	for(i=0;i+1<=len;i++) {
		printf("[%d, %d, %d]\n",
			lijst[i][0],
			lijst[i][1],
			lijst[i][2]);
	}
	return;
}

int gcd(int a, int b)
{
	int t;

	while (b > 0) {
		t = a % b;
		a = b;
		b = t;
	}
	return a;
}

/* This tests whether three out of four are divisible by the same prime.*
 * If so then the weighted projective space has a codimension 1 locus	*
 * of quotient singularities.						*/
int well_formed(int d1,  int d2,  int d3)
{
	if (gcd(d1,d2) > 1) return(0);
	if (gcd(d1,d3) > 1) return(0);
	if (gcd(d2,d3) > 1) return(0);
	return(1);
}

/* This tests whether there can be a quasi-smooth surface in the linear
 * system. The rule is that for each i either x_i^power can occur or that
 * x_i^power x_j should occur.... etc. See paper by Iano-Fletcher. 
 * ASSUMES: degree is bigger than d_i for all i. */
int d_suitable(int len)
{
	int i, j, k, success, sum;
	int test_one[3];
	int test_two[3];

	for (j=0;j<=2;j++) {
		test_two[j] = 0;
		test_one[j] = 0;
	}	

	for (i=0;i+1<=len;i++) {
		sum = lijst[i][0] + lijst[i][1] + lijst[i][2];
		for(j=0;j<=2;j++) {
			if (lijst[i][j] == 0) {
				test_two[j] = 1;
				k=0;
				while (k <= 2) {
					if ((k != j) && lijst[i][j] == 1)
						test_one[3-j-k] = 1;
					k++;
				}
			}
			if (sum == lijst[i][j]) test_one[j] = 1;
		}
	}
	success=0;
	for (j=0;j<=2;j++) {
		if (test_two[j] == 1) success++;
		if (test_one[j] == 1) success++;
	}
	if (success == 6) return 1;
	return 0;
}

/* This computes the dimension of the automorphism group. *
 * Note that we are computing the dimension of the cone.  */
int dim_aut(int d1, int d2, int d3)
{
	int totaal;
	totaal = 0;
	totaal = totaal + count_sum(d1,d2,d3,d1,0);
	totaal = totaal + count_sum(d1,d2,d3,d2,0);
	totaal = totaal + count_sum(d1,d2,d3,d3,0);
	return(totaal);
}

int main()
{
	int count,totaal,d1,d2,d3,degree,g;
	totaal=0;
/*
	count = count_sum(5,6,11,60,1);
	printf("Here is %d.\n",count);
	print_sum(count);
	d1 = well_formed(5,6,11);
	printf("Here is %d.\n",d1);
	d2 = d_suitable(count);
	printf("Here is %d.\n",d2);
	printf("%d \n ", gcd(1,2));
	exit(0);
*/
	printf("d1 d2 d3 d #f dim_aut g.\n");
	d1 = 1;
	while (d1 <= MAX) {
	  d2 = d1;
	  while (d2 <= MAX) {
	    d3 = d2;
	    while (d3 <= MAX) {
		degree = d1+d2+d3;
		while (degree <= MAX) {
			g = count_sum(d1,d2,d3,degree-d1-d2-d3,0);
if (g < 11) {
			count = count_sum(d1, d2, d3, degree,1);
			if (d_suitable(count)) {
				printf("%d %d %d %d %d ",
					d1,d2,d3,degree,count);
				count = dim_aut(d1,d2,d3);
				printf("%d ",count);
				printf("%d\n",g);
			}
}
			degree++;
	      };
	      d3++;
	    };
	    d2++;
	  };
	  d1++;
	};

	return(0);
}
