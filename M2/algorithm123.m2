load "Algorithm_1.2.1.m2"

multivariatePolynomialDivision = (p, gs...) -> (
    r := 0;
    f := p;
    numberOfPolynomials = #gs;
    
    while f != 0_R do (
	i := 0;
	divisionOccured := false;
	
	while (i < numberOfPolynomials and divisionOccured == false) do  (
	    (_, rem) = univariatePolynomialDivision(leadTerm(f), leadTerm(g));
	    if rem == 0 then (
		t := leadTerm(f)/leadTerm(g);
		oldG := gs[i]; 
		gs[i] = oldG + t;
		f = f - t*oldG
		divisionOccured = true;
		);
	    else i = i+1;
	    );
	if !divisionOccured then (
	    r = r + leadTerm(f);
	    f = f- leadTerm(f);
	    );
	);
      (gs, r)
     )
