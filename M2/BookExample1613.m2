runExample = () -> (
    load "algorithm161.m2";
    
    R = RR[x,y];
    f := (y^2+x-3, x^2-4*x*y-y^2-10*y+15);
    I := ideal f;

    regMatrix := solveQMatrix(R,I,1);
    xMatrix := solveQMatrix(R,I,x);
    xyMatrix := solveQMatrix(R,I,x*y);
    yMatrix := solveQMatrix(R,I,y);

    --only care about the x,y individual ones
    xEigens := (eigenvectors sub(xMatrix,CC))#0;
    yEigens := (eigenvectors sub(yMatrix,CC))#0;
    
    varietyList := for i from 0 to  #xEigens-1 list {xEigens#i, yEigens#i};
    print varietyList;
    )

    

    

    
