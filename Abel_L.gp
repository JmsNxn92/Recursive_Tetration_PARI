
/*Set series precision to 100 terms*/
\ps 100

/*Set the digit precision to 100; this number is set so high particularly for grabbing Taylor series; as high precision is needed.*/
 
\p 100

/*Set stack size to 750000000; the iterated logarithms in this code sure love to eat memory*/

default(parisize,750000000);

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
****It catches whether exp(out) will over flow most of the time.
*/

beta_function(z,l,n,{v=0}) =
{
	my(out = 0);
	if(v==0,
		for(i=0,n-1,
			if(abs(out)<= 1E4, 
				out = exp(out)/(exp(l*(n-i-z)) +1),
				if(abs(out) <= 1E8,
					return(exp(out)/(exp(l*(n-i-z)) +1)),
					return(out);
				);
			);
		);
		out,
		Ser(out,v);
		for(i=0,n-1,
			if(abs(polcoef(out,0,v))<= 1E4, 
				out = exp(out)/(exp(l*(n-i-z)) +1),
				if(abs(polcoef(out,0,v)) <= 1E8,
					return(exp(out)/(exp(l*(n-i-z)) +1)),
					return(out);
				);
			);
		);
		out;
	);
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
The requirement real(z) <= 10 can be lowered for faster convergence; if increased, there's a higher chance of an overflow. You don't really need higher though, as this will produce 100 digit precision.

*********Beware, the normalization constant may move when altering these values!
*/


tau(z,l,n, {v=0})={
	if(v==0,
		if((real(beta_function(z,l,n)) <= 1E6)&&(real(z) <=10), 
			log(1 + tau(z+1,l,n)/beta_function(z+1,l,n)) - log(1+exp(-l*z)),
			-log(1+exp(-l*z))
		),
		if((real(polcoef(beta_function(z,l,n,v),0,v)) <= 1E6)&&(real(polcoef(z,0,v))<=10),
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

The functional equation this satisfies is exp(Abl(z,l,n)) = Abl(z+1,l,n); and this function approaches that solution for n to infinity
*/

Abl(z,l,n,{v=0}) = {
	beta_function(z,l,n,v) + tau(z,l,n,v);
}

/* This is the old code, and isn't really needed; it was used in an attempt to work around divergent values; this was solved another way though
Instead use the legal Abl which is just a sum.
Abl(z,l,n,{v=0}) = {
	if(v==0,
		if(real(beta_function(z,l,n)) <= 1E8, 
			beta_function(z,l,n) + tau(z,l,n),
			exp(Abl(z-1,l,n))
		),
		if(real(polcoef(beta_function(z,l,n,v),0,v)) <= 1E8,
			beta_function(z,l,n,v) + tau(z,l,n,v),
			exp(Abl(z-1,l,n,v))
		)
	);
}
*/



/*
This is the pasted together beta_function. This function combines all the multipliers to find a new asymptotic solution.
n should be set to about 100 or higher; this will diverge fairly fast for large values of z, as it grows like tetration; and is more volatile than the previous asymptotic solutions.

v flags, again, whether you want to grab the taylor series or not..
*/

beta(z,n,{v=0}) = {
	beta_function(z,1/sqrt(1+z),n,v);
}

/*
This function acts similarly to tau; but it is constructed with respect to the new asymptotic solution.

I have used the function log_safe rather than log; to catch errors of when we dip too close to zero in the iteration.
This equates to when tau_B = - beta + small_error  which only happens really in the iteration, less so at the final result. I'm unsure why it does this.

The parameters in the test <=1E100 and <=8 are to make sure we don't overflow when we call beta. 
beta is only necessarily accurate for these paramaters; otherwise we open ourselves up to a flat line error.
*/

log_safe(z) = {
	if(abs(z) <= 1E-100000000,
		log(1E-100000000),
		log(z)
	);
}


/* the good log_safe, I haven't made this perfect yet.
log_safe(z) = {
	if(abs(z) <= 1E-100000000,
		my(q=abs(z));
		my(count =0);
		while(q <= 1E-100000000,
			count++;
			q = 10*q;
		);
		log(1E-100000000)*10^(count),
		log(z);
	);

}
*/

tau_B(z,n,{v=0}) = {
	if(v==0,
		if((real(beta(z,n)) <= 1E100)&& (real(z) <=8), 
			log_safe(1 + tau_B(z+1,n)/beta(z+1,n)) +beta_function(z,1/sqrt(2+z),n)- beta(z,n)-log(1+exp(-z/sqrt(2+z))),
			0
		),
		if((real(polcoef(beta(z,n,v),0,v)) <= 1E100) &&(real(polcoef(z,0,v)) <=15),
			log_safe(1 + tau_B(z+1,n,v)/beta(z+1,n,v)) +beta_function(z,1/sqrt(2+z),n,v)- beta(z,n,v)-log(1+exp(-z/sqrt(2+z))),
			0
		)
	);
}

/*
This is the tetration function not normalized.
*/

Tet_B(z,n,{v=0}) =
{
	beta(z,n,v) + tau_B(z,n,v);
}

/*
This is the normalized tetration function. 
The normalization constant is found by 2 - polrootsreal(Pol(Tet_B(1+v,100,v),v))
*/

sexp(z,n,{v=0}) =
{
	Tet_B(z+1.975055755684131009317363653179973725420301320301153615046531831010730376234679136573077233845701568,n,v);
}

/*
This produces a similar function to tau_B, but it's defined off the Abel function. 
This function is guessing the error of the same tetration; but executes different code.
*/

tau_Abl(z,n,{v=0}) = {
	if(v==0,
		if((real(Abl(z,1/sqrt(1+z),n)) <= 1E8)&& (real(z) <=10), 
			log(1 + tau_Abl(z+1,n)/Abl(z+1,1/sqrt(2+z),n)) +Abl(z,1/sqrt(2+z),n)- Abl(z,1/sqrt(1+z),n),
			0
		),
		if((real(polcoef(Abl(z,1/sqrt(1+z),n,v),0,v)) <= 1E8) &&(real(polcoef(z,0,v)) <=10),
			log(1 + tau_Abl(z+1,n,v)/Abl(z+1,1/sqrt(2+z),n,v)) +Abl(z,1/sqrt(2+z),n,v)- Abl(z,1/sqrt(1+z),n,v),
			0
		)
	);
}

/*
This is the un-normalized tetration.
*/

Tet_Abl(z,n,{v=0}) = {
	Abl(z,1/sqrt(1+z),n,v) + tau_Abl(z,n,v);
}

/* This is the normalization constant; if you've fiddled with the code it may change.

To generate a new normalization constant, run the command,

Norm_Constant = 2 - polrootsreal(Pol(Tet_Abl(1+v,100,z),v)) */

Norm_Constant = 1.975055755684131009317363653179973725420301320301153615046531831010730376234679136573077233845701568;


/*
This is the normalized tetration off the Abel method.
*/
sexp_Abl(z,n,{v=0}) = {
	Tet_Abl(z+Norm_Constant,n,v);
}


/*************************************************************************************/

/*
This function will produce 100 terms of the Taylor series about a point A. The value n is the depth of iteration inherited from beta.
*/

TAYLOR_SERIES(A,n) = {
	my(ser = vector(100,i,0));
	my(out = sexp(A+v,n,v));
	for(i=1,100, ser[i] = polcoef(out,i-1,v));
	return(ser);
}

/*
This sums the first 100 terms of the Taylor series about A. 
The variable C is an array of Taylor coefficients; these can be grabbed with TAYLOR_SERIES.
*/

SUM_TAYLOR(z,A,C) = {
	sum(j=1,100, C[j]*(z-A)^(j-1));
}

/*
This function will evaluate the Taylor series if real(z) <= 0.5; and other wise will iterate the exponential. You can play with this if you want.
This produces Tetration. It requires an array of coefficients, and the point they're centered about.
WARNING: Don't forget about the radius of convergence when using this.
*/

sexp_T(z,A,C) = {
	if(real(z) <= -0.5, 
		SUM_TAYLOR(z,A,C), 
		exp(sexp_T(z-1,A,C))
	);
}

/*
This is the logarithmic solution to Tetration. This code will work fine away from the real line and on the real-line precisely.
This code will diverge in a neighborhood of the real-line, and produce incorrect results. To handle the neighborhood of the real-line, use the taylor series approach.

This means, to get the value of say Tet(1+0.8*I,100); it is better to get Tet(1+z,100,z) and expand the Taylor Series output.

This function is out of date; but still works if you want to use it; it has trouble in the complex plane; as the iterated logarithms tend to diverge between +/- Pi.
*/

/*
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
*/

/*
This is the normalized tetration function; which we call the super-exponential.
Again, this function is out of date.
*/

/*
Sexp(z,n,{v=0}) = {
	if(v==0,
		if(real(z) <= 0.5,
			if(real(z) < -0.5, 
				log(Sexp(z+1,n,v)),
				Tet(z+1.969637739698065306544624079350257708852542229771084623924562193889980567396859073585886113019833431,n,v)
			),
			exp(Sexp(z-1,n,v))
		),
		if(real(polcoef(z,0,v)) <= 0.5,
			if(real(polcoef(z,0,v)) < -0.5, 
				log(Sexp(z+1,n,v)),
				Tet(z+1.969637739698065306544624079350257708852542229771084623924562193889980567396859073585886113019833431,n,v)
			),
			exp(Sexp(z-1,n,v))
		)
	);	
}
*/

/*
This was a more involved version of sexp; I was troubleshooting whether it would feasibly run faster/more accurately; it doesn't

Sexp_Abl(z,n,{v=0}) = {
	if(v==0,
		if(real(z) <= 0.5,
			if(real(z) < -0.5, 
				log(Sexp_Abl(z+1,n,v)),
				Tet_Abl(z+1.969637740420574527022680792269958316936039199716124642340271732596898734287065726696866509279661244,n,v)
			),
			exp(Sexp_Abl(z-1,n,v))
		),
		if(real(polcoef(z,0,v)) <= 0.5,
			if(real(polcoef(z,0,v)) < -0.5, 
				log(Sexp_Abl(z+1,n,v)),
				Tet_Abl(z+1.969637740420574527022680792269958316936039199716124642340271732596898734287065726696866509279661244,n,v)
			),
			exp(Sexp_Abl(z-1,n,v))
		)
	);	
}
*/
