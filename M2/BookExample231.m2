R = QQ[z]

p = z^4+3*z^3-3*z^2+3*z-4

pc = -2*z^4+3*z^3-6*z^2+5*z-4

g = z^2+1

(q,r)= quotientRemainder(p,g)

(qc, rc) = quotientRemainder(pc,g)
