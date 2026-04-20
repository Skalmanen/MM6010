runExample = () ->(
    R = QQ[z];

    g := z^3-3*z^2+2*z;

    I = ideal g;

    q := R / I;

    basisList := flatten entries basis q;
    polyBasis := for b in basisList list lift(b, R); -- ensuring that the actualy polynomials exists in R and not only in q

    rows := for b in  basisList list (
	currentPoly = z *lift(b,R);
    
	r := currentPoly % I;

	--taking #1 here since coefficients returns a bundle of want the coefficients, no the monomials
	coeffsMatrix := (coefficients(r, Monomials => polyBasis))#1;
	flatten entries coeffsMatrix
	);

    print rows;
    zmatrix := matrix rows;
    print zmatrix;

    finalZMatrix := zmatrix^2+zmatrix+id_(target zmatrix);
    print finalZMatrix;
)


 
