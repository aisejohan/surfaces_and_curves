#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

/* really dumb. */
#define MAX	500
#define MAX_g	11
#define MAX_pg	11

int storage=0;
int **lijst;

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

/* This tests whether there can be a quasi-smooth curve in the linear
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
					if ((k != j) && lijst[i][k] == 1)
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

int hilbert_function(int d1, int d2, int d3, int d, int i)
{
	int phi;

	if (i < 0) return(0);

	phi = count_sum(d1,d2,d3,i,0);
	if (i >= d-d1) phi -= count_sum(d1,d2,d3,i-d+d1,0);
	if (i >= d-d2) phi -= count_sum(d1,d2,d3,i-d+d2,0);
	if (i >= d-d3) phi -= count_sum(d1,d2,d3,i-d+d3,0); 
	if (i >= 2*d-d1-d2) phi += count_sum(d1,d2,d3,i-2*d+d1+d2,0);
	if (i >= 2*d-d1-d3) phi += count_sum(d1,d2,d3,i-2*d+d1+d3,0);
	if (i >= 2*d-d2-d3) phi += count_sum(d1,d2,d3,i-2*d+d2+d3,0);
	if (i >= 3*d-d1-d2-d3) phi -= count_sum(d1,d2,d3,i-3*d+d1+d2+d3,0);
	return(phi);
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

void list_curves()
{
	int count,d1,d2,d3,degree,g;

	printf("d1 d2 d3 d #f dim_aut g.\n");
	d1 = 1;
	while (d1 <= MAX) {
	  d2 = d1;
	  while (d2 <= MAX) {
	    d3 = d2;
	    while (d3 <= MAX) {
	      if (well_formed(d1,d2,d3)) {
		degree = d1+d2+d3;
		while (degree <= MAX) {
		  g = count_sum(d1, d2, d3, degree-d1-d2-d3, 0);
		  if (g < MAX_g) {
		    count = count_sum(d1,d2,d3,degree,1);
		    if (d_suitable(count)) {
		        printf("%d %d %d %d %d ",d1,d2,d3,degree,count);
		        count = dim_aut(d1,d2,d3);
		        printf("%d ",count);
		        printf("%d\n",g);
		    }
		  }
		  degree++;
	        }
	      }
	      d3++;
	    }
	    d2++;
	  }
	  d1++;
	}

	exit(0);
}

void list_double_covers()
{
	int count,d1,d2,d3,degree,pg,i,h11;

	printf("d1 d2 d3 d #f dim_aut i pg h11.\n");
	d1 = 1;
	while (d1 <= MAX) {
	  d2 = d1;
	  while (d2 <= MAX) {
	    d3 = d2;
	    while (d3 <= MAX) {
	      if (well_formed(d1,d2,d3)) {
		degree = 2*d1+2*d2+2*d3;
		while (degree <= MAX) {
		  i = (degree/2)-d1-d2-d3;
		  pg = hilbert_function(d1,d2,d3,degree,i);
		  if (pg < MAX_pg) {
		    count = count_sum(d1, d2, d3, degree,1);
		    if (d_suitable(count)) {
		      printf("%d %d %d %d %d ",d1,d2,d3,degree,count);
		      count = dim_aut(d1,d2,d3);
		      printf("%d ",count);
		      printf("%d ",i);
		      printf("%d ",pg);
		      h11 = hilbert_function(d1,d2,d3,degree,i+degree);
		      printf("%d\n",h11);
		    }
		  }
		  degree = degree + 2;
		}
	      }
	      d3++;
	    }
	    d2++;
	  }
	  d1++;
	}

	exit(0);
}

int main(void )
{
	list_double_covers();

	return(0);
}
