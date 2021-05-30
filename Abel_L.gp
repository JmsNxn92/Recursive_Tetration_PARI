
/*Set series precision to 50 terms*/
\ps 50 

/*Set the digit precision to 500; this number is set so high particularly for grabbing Taylor series; as high precision is needed.*/
 
\p 500 

/*
This function constructs the approximate Abel function
The variable z is the main variable we care about; this grows like tetration, so expect overflow errors for large arguments.
The variable l is the multiplier of the approximate Abel function
the variable n is the depth of iteration required
n should be set to 100 or higher for better precision, but produces fairly good accuracy for about n=15 
The functional equation this satisfies is exp(beta_function(z,l,n))/(1+exp(-l*z)) = beta_function(z+1,l,n); and this program approaches the solution for n to infinity
*/

beta_function(z,l,n) =
{
	my(out = 0);
	for(i=0,n-1,
		out = exp(out)/(exp(l*(n-i-z)) +1));
	out;
}

/*
This function is the error term between the approximate Abel function and the actual Abel function
The variable z is the main variable we care about
The variable l is the multiplier
The variable n is the depth of iteration inherited from beta_function
n can be set about 100, still; but 15 or 20 is sufficient for small precision. 
Making 0.0001 smaller produces better accuracy (about 0.00000001 is the limit); too small and it'll overflow; too large and it'll diverge
*/


tau(z,l,n)={
	if(1/real(beta_function(z,l,n)) <= 0.0001, 
		-log(1+exp(-l*z)),
		log(1 + tau(z+1,l,n)/beta_function(z+1,l,n)) - log(1+exp(-l*z))
	)
}

/*
This is the actual Abel function
The variable z is the main variable we care about
The variable l is the multiplier
The variable n is the depth of iteration inherited from beta_function
The functional equation this satisfies is exp(Abl(z,l,n)) = Abl(z+1,l,n); and this function approaches that solution for n to infinity and 0.0001 set smaller
*/

Abl(z,l,n) = {
	beta_function(z,l,n) + tau(z,l,n)
}

/*
This is the pasted together beta_function. This function combines all the multipliers to find a new asymptotic solution.
n should be set to about 100 or higher; this will diverge fairly fast for large values of z, as it grows like tetration; and is more volatile than the previous asymptotic solutions.
*/

beta(z,n) = {
	beta_function(z,1/sqrt(1+z),n);
}

/*
This is the logarithmic solution to Tetration. This code will work fine away from the real line and on the real-line precisely.
This code will diverge in a neighborhood of the real-line, and produce incorrect results. To handle the neighborhood of the real-line, use the taylor series approach.
*/

Tet(z,n) ={
	if(1/real(beta(z,n)) <= 0.0000001,
		beta(z,n),
		log(Tet(z+1,n))
	)
}

/*
This function estimates the depth of recursion needed to produce an accurate Taylor series about A
Setting 0.0001 smaller produces better accuracy; capping out at about 0.00000001.
*/

Tet_GRAB_k(A,n) ={
	my(k=0);
	while( 1/real(beta(A+k,n)) >= 0.0001, k++); 
	return(k);
}

/*
This function works exactly as Tet does; but it needs an input to tell it how deep to do the recursion, which is what k is for.
This function will allow you to grab Taylor series, where Tet doesn't.
*/

Tet_taylor(z,n,k) = {
	my(val = beta(z+k,n));
	for(i=1,k,val = log(val));
	return(val);
}

/*
This function will produce 50 terms of the Taylor series about a point A. The value n is the depth of iteration inherited from beta.
*/

TAYLOR_SERIES(A,n) = {
	my(ser = vector(50,i,0));
	my(k=Tet_GRAB_k(A,n));
	for(i=1,50, ser[i] = polcoeff(Tet_taylor(A+v,n,k),i-1,v));
	return(ser);
}

/*
This sums the first 50 terms of the Taylor series about A. 
The variable C is an array of Taylor coefficients; these can be grabbed with TAYLOR_SERIES.
*/

SUM_TAYLOR(z,A,C) = {
	sum(j=1,50, C[j]*(z-A)^(j-1));
}

/*
This function will evaluate the Taylor series if real(z) <= 1.5; and other wise will iterate the exponential. You can play with this if you want.
This produces Tetration. It requires an array of coefficients, and the point they're centered about.
WARNING: Don't forget about the radius of convergence when using this.
*/

TET_TAYLOR_ITER(z,A,C) = {
	if(real(z) <= 1.5, SUM_TAYLOR(z,A,C), exp(TET_TAYLOR_ITER(z-1,A,C)));
}
