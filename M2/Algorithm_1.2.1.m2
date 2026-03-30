univariatePolynomialDivision = (p,g) -> (
    q := 0;
    r := p;
    while (r != 0_R and degree(g) < degree(r)) do (
	t := leadTerm(r) / leadTerm(g);
	q = q+t;
	r = r-t*g;
	);
    (q,r)
    )

runUnivariateTest = () -> (
    R = QQ[x];
    p := x^6 - 1;
    g := x^4 - 1;
    
    (resultQ, resultR) := univariatePolynomialDivision(p,g);
    
    print(resultQ == x^2);
    print(resultR == x^2 - 1);
)
