
/*Set series precision to 100 terms*/
\ps 100

/*Set the digit precision to 100; this number is set so high particularly for grabbing Taylor series; as high precision is needed.*/
 
\p 100

/*
This function constructs the approximate Abel function
The variable z is the main variable we care about; this grows like tetration, so expect overflow errors for large arguments.
The variable l is the multiplier of the approximate Abel function
the variable n is the depth of iteration required

The variable v is a marker as to whether you want to grab taylor series or not.
For instance, writing beta_function(Pi+z,1,100,z), will grab the taylor series with respect to z about the point Pi.
If left blank, then the function will run pointwise as t_COMPLEX, and grabbing taylor series will produce an error.

n should be set to 100 or higher for better precision, but produces fairly good accuracy for about n=15 
The functional equation this satisfies is exp(beta_function(z,l,n))/(1+exp(-l*z)) = beta_function(z+1,l,n); and this program approaches the solution for n to infinity

****Added an additional if statement to avoid overflow errors.
****It catches whether exp(out) will over flow.
*/

beta_function(z,l,n,{v=0}) =
{
	my(out = 0);
	if(v==0,
		for(i=0,n-1,
			if(abs(out)<= 1E6, 
				out = exp(out)/(exp(l*(n-i-z)) +1),
				return(out);
			)
		);
		out,
		Ser(out,v);
		for(i=0,n-1,
			if(abs(polcoef(out,0,v))<= 1E6, 
				out = exp(out)/(exp(l*(n-i-z)) +1),
				return(out);
			)
		);
		out;
	)
	
}



/*
This function is the error term between the approximate Abel function and the actual Abel function
The variable z is the main variable we care about
The variable l is the multiplier
The variable n is the depth of iteration inherited from beta_function

The variable v is a marker as to whether you want to grab taylor series or not.
For instance, writing tau(Pi+z,1,100,z), will grab the taylor series with respect to z about the point Pi.
If left blank, then the function will run pointwise as t_COMPLEX, and grabbing taylor series will produce an error.

n can be set about 100, still; but 15 or 20 is sufficient for small precision. 
Making 1E6 larger produces better accuracy (about 1E8 is the limit); too large and it'll overflow; too small and it'll diverge
*/


tau(z,l,n, {v=0})={
	if(v==0,
		if(real(beta_function(z,l,n)) <= 1E6, 
			log(1 + tau(z+1,l,n)/beta_function(z+1,l,n)) - log(1+exp(-l*z)),
			-log(1+exp(-l*z))
		),
		if(real(polcoef(beta_function(z,l,n,v),0,v)) <= 1E6,
			log(1 + tau(z+1,l,n,v)/beta_function(z+1,l,n,v)) - log(1+exp(-l*z)),
			-log(1+exp(-l*z))
		)
	);
}

/*
This is the actual Abel function
The variable z is the main variable we care about
The variable l is the multiplier
The variable n is the depth of iteration inherited from beta_function

The variable v is a marker as to whether you want to grab taylor series or not.
For instance, writing Abl(Pi+z,1,100,z), will grab the taylor series with respect to z about the point Pi.
If left blank, then the function will run pointwise as t_COMPLEX, and grabbing taylor series will produce an error.

The functional equation this satisfies is exp(Abl(z,l,n)) = Abl(z+1,l,n); and this function approaches that solution for n,k to infinity
*/

Abl(z,l,n,{v=0}) = {
	beta_function(z,l,n,v) + tau(z,l,n,v)
}

/*
This is the pasted together beta_function. This function combines all the multipliers to find a new asymptotic solution.
n should be set to about 100 or higher; this will diverge fairly fast for large values of z, as it grows like tetration; and is more volatile than the previous asymptotic solutions.

v flags, again, whether you want to grab the taylor series or not.
*/

beta(z,n,{v=0}) = {
	beta_function(z,1/sqrt(1+z),n,v);
}

/*
This is the logarithmic solution to Tetration. This code will work fine away from the real line and on the real-line precisely.
This code will diverge in a neighborhood of the real-line, and produce incorrect results. To handle the neighborhood of the real-line, use the taylor series approach.

This means, to get the value of say Tet(1+0.8*I,100); it is better to get Tet(1+z,100,z) and expand the Taylor Series output.
*/

Tet(z,n, {v=0}) ={
	if(v==0,
		if(real(beta(z,n)) <= 1E6,
			log(Tet(z+1,n)),
			beta(z,n)
		),
		if(real(polcoef(beta(z,n,v),0,v)) <= 1E6,
			log(Tet(z+1,n,v)),
			beta(z,n,v)
		)
	);
}

/*
This is the normalized tetration function; which we call the super-exponential. The normalization constant is accurate enough for high-precision.
*/

Sexp(z,n,{v=0}) = {
	Tet(z+1.969637739698065306544624079350257708852542229771084623924562193889980567396859073585886113019833431,n,v);
}
}

/*
This function will produce 50 terms of the Taylor series about a point A. The value n is the depth of iteration inherited from beta.
*/

TAYLOR_SERIES(A,n) = {
	my(ser = vector(100,i,0));
	my(out = Tet(A+v,n,v));
	for(i=1,100, ser[i] = polcoef(out,i-1,v));
	return(ser);
}

/*
This sums the first 50 terms of the Taylor series about A. 
The variable C is an array of Taylor coefficients; these can be grabbed with TAYLOR_SERIES.
*/

SUM_TAYLOR(z,A,C) = {
	sum(j=1,100, C[j]*(z-A)^(j-1));
}

