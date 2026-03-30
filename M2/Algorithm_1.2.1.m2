univariatePolynomialDivision = (p,g) -> (
    q := 0;
    r := p;
    while (r != 0_R  and degree(g) < degree(r)) do (
	t := leadTerm(r) / leadTerm(g);
	q = q+t;
	r = r-t*g;
	);
    (q,r)
    )
