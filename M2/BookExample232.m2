R=QQ[c_1,c_2,c_3,c_4, MonomialOrder=>Lex]
M = matrix{{1,-2,0,-1},{4,0,7,3}, {7,-6,7,0}}
C = matrix{{c_1},{c_2},{c_3},{c_4}}
I = ideal(M*C)
GI = value toString gens gb I
(Ce,Me) = value toString coefficients(GI, Monomials=>{c_1,c_2,c_3,c_4})

Mr = value toString transpose Me

inverse(substitute(submatrix(Mr, {0,1}),QQ))*Mr
