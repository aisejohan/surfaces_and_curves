myN = 20;
myepsilon = 0.0001;

findsmall(f,p,n) = 
{
	local(d,c,i,polygon,j,e,g);
	if(isprime(p),,error("Second argument should be a prime."));
	d = poldegree(f);
	g = x^d;
	polygon = newtonpoly(f,p);
	for(i=0,d-1,
		c = polcoeff(f,i);
		e = floor(sum(j=1,d-i,polygon[d+1-j]));
		c = c/p^e;
		c = lift(Mod(c,p^n));
		if((c>(p^n)/2),c=c-p^n);
		g = g + c*p^e*x^i;
	);
	return(g);
}

check_symmetry(f,p,w) =
{
	local(d,c);
	if(isprime(p),,error("Second argument should be a prime."));
	d = poldegree(f);
	c = polcoeff(f,0)/p^(w*d/2);
	if((c^2 - 1),
		print("Not a Weil polynomial.");
		return(0)
	);
	if((c == 1),
		if(x^d*p^(-w*d/2)*subst(f,x,p^w/x)-f,
			print("Not symmetric.");
			return(0)
		,
			return(1)
		)
	);
	if((c == -1),
		if(x^d*p^(-w*d/2)*subst(f,x,p^w/x)+f,
			print("Not symmetric.");
			return(0)
		,
			return(1)
		)
	);
}

/* Finds Weil polynomial. */
findweil(f,p,initial) =
{
	local(d,polygon,w,n,g,lijst,success);
	if(isprime(p),,error("Second argument should be a prime."));
	if(initial,,initial=1+floor(log(poldegree(f)*p)/log(p)));
	d = poldegree(f);
	polygon = newtonpoly(f,p);
	w = 2*sum(i=1,d,polygon[i])/d;
	print1("The weight is ",w,".");
	n=myN+initial+1;
	while(n>=initial,
		n = n - 1;
		g = findsmall(f,p,n);
		lijst = abs(polroots(g));
		success = 1;
		for(i=1,matsize(lijst)[1],
			if(((p^(w/2)-myepsilon > lijst[i]) ||
				(lijst[i] > p^(w/2)+myepsilon)),
				success = 0
			)
		);
		if(success,
			print1(" Something found... ");
			if(check_symmetry(g,p,w),
				print(" Success!");
				return(g)
			)
		)
	);
	print("No luck this time!");
}
