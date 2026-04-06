multivariatePolynomialDivision = (p, gs) -> (
    r := 0_(ring p);
    f := p;
    qs := for i from 0 to #gs-1 list 0_(ring p); --generating dummy variable for result

    print qs;
    while f != 0_(ring p) do (
	i := 0;
	divisionOccured := false;
 	--print "Outer  loop";
	
	while (i < #gs and divisionOccured == false) do (
	    --print "Inner loop";
	    lmf := leadMonomial f;
	    lmg := leadMonomial (gs#i);
	    
	    if (lmf % lmg == 0) then (
		monomialQuotient := lmf // lmg;

		coefficientQuotient := leadCoefficient f / leadCoefficient gs#i;

		t := coefficientQuotient * monomialQuotient;
    
		qs = replace(i, qs#i + t, qs);
		f = f - t*gs#i;
		
		divisionOccured = true;
		--print "we were divisible";
	    )
	    else (
	        i = i+1;

		--print "not divisible, function further along";
	    )
	);
    	
        if divisionOccured == false then (
	    r = r + leadTerm(f);
	    f = f- leadTerm(f);
	   --print "division did not occur, subtracting leading term and continuing";
	);
    );
    return (qs, r)
)

-- This below code is useless, it was an attempt to get around the fact that i couldnt devide with the leadTerm.
-- However, I needed to seperate out the divisions as well as use "//" instead of "/"
leadingTermDevision = (p,g) -> (
    lmP := leadMonomial(p);
    lmG := leadMonomial(g);
    lcP := leadCoefficient(p);
    lcG := leadCoefficient(g);
    exponentsP := flatten exponents lmP;
    exponentsG := flatten exponents lmG;
    
    divisible := all(#exponentsP, i -> exponentsP#i >= exponentsG#i); --checking that all exponents are greater or equal
    result := 0_(ring p);

    if divisible then(
	coefficientResult := lcP / lcG;
	
	exponentsQ := for i from 0 to (#exponentsP-1) list (exponentsP#i - exponentsG#i); --doing the division

	ringVariables := toList gens ring p; 
	
	Q := product(#ringVariables, i -> ringVariables#i ^ (exponentsQ#i)); --generating the the monomial
	
	result = coefficientResult * Q;
     );

    return (divisible, result);
 )

 runMultimultivariateTest = () -> (
     R = QQ[x,y,MonomialOrder=>Lex];
     p := -x^2*y-x^2+x*y+y^3;
     g1 := x*y-x+y-1;
     g2 := x+y+1;

     (resultG1, resultRemainder1) := multivariatePolynomialDivision(p,{g1,g2});
     (resultG2, resultRemainder2) := multivariatePolynomialDivision(p,{g2,g1});

     print "-----------{g1,g2}------------";
     print resultG1;
     print resultRemainder1;
     print (resultG1#0 == -x+4);
     print (resultG1#1 == -2*x+5);
     print (resultRemainder1 == y^3-9*y-1);
     
     print "-----------{g2,g1}------------";
     print resultG2;
     print resultRemainder2;
     print (resultG2#1 == 0);
     print (resultG2#0 == -x*y-x+y^2+3*y+1);
     print (resultRemainder2 == -4*y^2-4*y-1);
)
     
