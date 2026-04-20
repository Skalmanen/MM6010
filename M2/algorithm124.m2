load "algorithm123.m2";

buchberger = ps -> (
    g := ps;
    gHat := {0_(ring  ps#0)};

    while g != gHat do (
	gHat = g;

	-- generating all the possible pairs by looping through i and then from j -> end
        polynomialPairs := flatten for i from 0 to #gHat-1 list for j from i+1 to #gHat-1 list (gHat#i, gHat#j);
	  	  
	i := 0;
	while i < #polynomialPairs do (
	    p := (polynomialPairs#i)#0;
	    q := (polynomialPairs#i)#1;
	    
	    xDelta := lcm(leadMonomial(p), leadMonomial(q));

	    s := xDelta // leadTerm(p) * p - xDelta // leadTerm(q) * q;

	    (qs, r) = multivariatePolynomialDivision(s, gHat);

	    if r != 0_(ring ps#0) then (
		g = append(gHat,r);
		);
	    
	    i = i+1;
	    );
	);
    print g;
    return g;
    )

testAlgorithm = () -> (
    R = QQ[x,y, MonomialOrder => GLex];

    p1 := x^2 + y;
    p2 := x^4 + 2*x^2*y + y^2 + 3;

    myBasis := gb ideal (buchberger({p1, p2}));
    correctBasis := gb ideal(p1, p2);
    areTheyEqual := myBasis == correctBasis;
    
    print "My basis:";
    print myBasis;

    print "Built-in basis:";
    print gens correctBasis;

    print "Are they equal?";
    print areTheyEqual;
    )
    
--again, this is not needed and M2 implements LCM.
findXDelta = (p,q) -> (
    alpha := flatten exponents leadMonomial(p);
    gamma := flatten exponents leadMonomial(q);

    delta := for i from 0 to (#alpha-1) list (max(alpha#i,gamma#i)); --getting the max of each variables exponent
    
    ringVariables := toList gens ring p;
    xDelta := product(#ringVariables, i -> ringVariables#i ^ (delta#i)); -- reconstructing the variable

    return xDelta
    )

testGamma = () -> (
    R := RR[x,y, MonomialOrder => GLex];
    p := x^3*y^2-x^2*y^3+x;
    q := 3*x^4*y +y^2;

    xDelta := findXDelta(p,q);
    xDeltaShouldBe := x^4*y^2;

    print "-----------gamma results------------";
    print xDelta;
    print xDeltaShouldBe;
    print("the result is " | toString(xDelta == xDeltaShouldBe))
    )
