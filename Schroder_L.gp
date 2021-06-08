/*
This function constructs the approximate Schroder function
The variable w is the main variable we care about
The variable l is the multiplier of the approximate Schroder function
the variable n is the depth of iteration required

The variable v is a marker as to whether you want to grab taylor series or not.
For instance, writing Phi(Pi+z,1,100,z), will grab the taylor series with respect to z about the point Pi.
If left blank, then the function will run pointwise as t_COMPLEX, and grabbing taylor series will produce an error.

n can be set to 100, but produces enough accuracy for about 15 
The functional equation this satisfies is exp(Phi(w,l,n))/(1+w) =Phi(exp(-l)*w,l,n); and this program approaches the solution for n to infinity
*/

Phi(w,l,n,{v=0}) =
{
	my(out = 0);
	if(v==0,
		for(i=0,n-1,
			if(abs(out)<=1E6,
				out = exp(out)/(exp(l*(n-i))*w +1),
				return(out)
			)
		);
		out,
		Ser(out,v);
		for(i=0,n-1,
			if(abs(polcoef(out,0,v))<= 1E6, 
				out = exp(out+O(v^100))/(exp(l*(n-i))*w +1),
				return(out);
			)
		);
		out;
	);
}

/*
This function is the error term between the approximate Schroder function and the actual Schroder function
The variable w is the main variable we care about
The variable l is the multiplier
The variable n is the depth of iteration inherited from Phi
n can be set about 100, still; but 15 or 20 is more optimal. 

The variable v is a marker as to whether you want to grab taylor series or not.
For instance, writing err(Pi+z,1,100,z), will grab the taylor series with respect to z about the point Pi.
If left blank, then the function will run pointwise as t_COMPLEX, and grabbing taylor series will produce an error.

Setting the variable k above 10 will produce overflow errors. 
Precision of about 10 digits is acquired at k = 8 or 9
*/

err(w,l,n,{v=0})={
	if(v==0,
    	if(real(Phi(w,l,n)) <=1E6, 
    		log(1 + err(exp(-l)*w,l,n)/Phi(exp(-l)*w,l,n)) - log(1+w),
    		-log(1+w)
    	),
    	if(real(polcoef(Phi(w,l,n,v),0,v)) <= 1E6,
    		log(1 + err(exp(-l)*w,l,n,v)/Phi(exp(-l)*w,l,n,v)) - log(1+w),
    		-log(1+w)
    	)
    );
}

/*
This is the actual Schroder function
The variable w is the main variable we care about
The variable l is the multiplier
The variable n is the depth of iteration inherited from Phi

The variable v is a marker as to whether you want to grab taylor series or not.
For instance, writing Sch(Pi+z,1,100,z), will grab the taylor series with respect to z about the point Pi.
If left blank, then the function will run pointwise as t_COMPLEX, and grabbing taylor series will produce an error.

The functional equation this satisfies is exp(Sch(w,l,n)) = Sch(exp(-l)*w,l,n); and this function approaches that solution for n to infinity
*/

Sch(w,l,n,{v=0}) ={
	if(v==0,
		if(real(Phi(w,l,n)) <= 1E6, 
			Phi(w,l,n) + err(z,l,n),
			exp(Sch(exp(l)*w,l,n))
		),
		if(real(polcoef(Phi(w,l,n,v),0,v)) <= 1E6,
			Phi(w,l,n,v) + err(z,l,n,v),
			exp(Sch(exp(l)*w,l,n,v))
		)
	);
}
