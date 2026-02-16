R = QQ[x_1..x_4, MonomialOrder => Lex]

p1 = x_1-x_2**2+x_3

p2 = x_1**2+x_2*x_3**2+x_4

I = ideal(p1,p2)

G = gens gb I

polyToDivide = x_1**2+x_2**2+x_3**2+x_4**2-1

rem = polyToDivide%I

print("the rem is: ")
print(rem)

