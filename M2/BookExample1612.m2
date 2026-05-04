runExample = () -> (
    R = RR[z];

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

    zmatrix := matrix rows;

    finalZMatrix := zmatrix^2+zmatrix+id_(target zmatrix);
    return finalZMatrix;
)


 runExampleWithFunction = () -> (
    R = RR[z];

    g := z^3-3*z^2+2*z;

    I = ideal g;

    load "algorithm161.m2";

    zmatrix := solveQMatrix(R, I, z);

    finalZMatrix := zmatrix^2+zmatrix+id_(target zmatrix);
    return finalZMatrix;
)

testFunction = () -> (
    preBuiltMatrix := runExample();
    matrixFromFunction := runExampleWithFunction();

    print "Matrix from the example";
    print preBuiltMatrix;

    print "Matrix using refactored function";
    print matrixFromFunction;
 )
