myN = 20;
myepsilon = 0.0001;

findsmall(f,p,n) = 
{
	local(d,c,i,polygon,j,e,g);
	if(isprime(p),,error("Second argument should be a prime."));
	d = poldegree(f);
	g = x^d;
	polygon = newtonpoly(f,p);
	for(i=0,d-1,\
		c = polcoeff(f,i);\
		e = floor(sum(j=1,d-i,polygon[d+1-j]));\
		c = c/p^e;\
		c = lift(Mod(c,p^n));\
		if((c>(p^n)/2),c=c-p^n);\
		g = g + c*p^e*x^i;\
	);
	return(g);
}

check_symmetry(f,p) =
{
	local(d,c);
	if(isprime(p),,error("Second argument should be a prime."));
	d = poldegree(f);
	c = polcoeff(f,0)/p^d;
	if((c^2 - 1),error("Not a Weil polynomial."));
	if((c == 1),\
		if(x^d*p^(-d)*subst(f,x,p^2/x)-f,\
			print("Not symmetric.")\
		)\
	);
	if((c == -1),\
		if(x^d*p^(-d)*subst(f,x,p^2/x)+f,\
			print("Not symmetric.")\
		)\
	);
}

/* Finds Weil polynomial of "weight" 2. */
findweil(f,p,initial) =
{
	local(n,g,lijst,success);
	if(isprime(p),,error("Second argument should be a prime."));
	if(initial,,initial=1+floor(log(poldegree(f)*p)/log(p)));
	n=myN+1;
	while(n>=initial,\
		n = n - 1;\
		g = findsmall(f,p,n);\
		lijst = abs(polroots(g));\
		success = 1;\
		for(i=1,matsize(lijst)[2],\
			if(((p-myepsilon > lijst[i]) || (lijst[i] > p+myepsilon)),\
				success = 0\
			)\
		);\
		if(success,\
			check_symmetry(g,p);\
			return(g))\
	);
	print("No luck this time!");
}
