solveSystemOfPolynomialEquations = (theRing, theIdeal) -> (
    varsToLoopOver := gens theRing;
    
   -- get the list of eigenvalues
    eigenList := for var in varsToLoopOver list (
	qMatrix := solveQMatrix(theRing, theIdeal, var);
	
	 --ensure that the eigenvector can be solved and save it to the list
         toList eigenvalues(sub(qMatrix,CC))
	);

    -- use fold to go through all elements in the eigenList  outwards
    candidatePoints := fold((L1, L2) -> (
	    --use table to create a table of all the pairs and flatten down
	    flatten table(L1, L2, (a, b) -> (
		    --if we are past first instance, we append instead of create a list
		    if instance(a, List) then a | {b} else {a, b}
		    )	  
		)
	    ),
	eigenList);

    epsilon := 1e-6;
    functionList := flatten entries gens theIdeal;
    
    validSolutions := select(candidatePoints, pt -> (
	    -- ensure that x_n maps to the the nth point
	    subMap = apply(#varsToLoopOver, i -> varsToLoopOver#i => pt#i);
	    -- sub in the values
	    all(functionList, p -> abs(sub(p, subMap)) < epsilon)
	    ));

    -- something is wrong with my solveQMatrix which results in cluttered returns.
    -- so for now just to check if it works, we remove the duplicates
    uniqueSolutions := {};
    scan(validSolutions, v -> (
	    if not any(uniqueSolutions, u -> (
		    dist := sqrt(sum(0..#v-1, i -> abs(v#i - u#i)^2));
		    dist < epsilon
		    )) then uniqueSolutions = append(uniqueSolutions, v);
	    ));

    uniqueSolutions
    )


testAlgorithm2D = () -> (
    R := CC[x, y];
    I := ideal(x^2 + y^2 - 5, x*y - 2);

    myResults := solveSystemOfPolynomialEquations(R, I);
    solution := {{-2, -1}, {-1, -2}, {2, 1}, {1, 2}};
    
    print ("The 2D test output: " | toString  myResults);
    
    print ("Verified solution: " | toString solution);
)

testAlgorithm4D = () -> (
    R4 := CC[x1, x2, x3, x4];
    
    p1 := x1^2 + x2^2 - 2;
    p2 := x1 - x2;
    p3 := x3^2 + x4^2 - 13;
    p4 := x3*x4 - 6;
    I4 := ideal(p1, p2, p3, p4);

    results4D := solveSystemOfPolynomialEquations(R4, I4);
    
    if #results4D == 8 then (
        print "SUCCESS: Found all 8 points in 4D space.";
    ) else (
        print ("FAILURE: Expected 8 points, but found " |  #results4D);
    );
)

solveQMatrix = (theRing, theIdeal, thePolynomial) ->(
    q := theRing / theIdeal;

    basisList := flatten entries basis q;
    polyBasis := for b in basisList list lift(b,theRing); --basis to R

    rows := for b in polyBasis list (
	currentPoly := thePolynomial * b;
	r := currentPoly % theIdeal;

	coeffsMatrix := last coefficients(r, Monomials => polyBasis);
	flatten entries coeffsMatrix
	);

    return transpose matrix rows;
    )

print "-------2D CASE-------";
testAlgorithm2D();
print "------4D CASE---------";
testAlgorithm4D();
print "-------DONE-------";



    
