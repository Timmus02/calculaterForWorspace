from sympy import symbols, Matrix, cos, sin, pprint, simplify
import numpy as np
from JacobiBibliothek import *

l2 = "l2"

l1 = 100
l3 = 600
l4 = 600
l5 = 600
l6 = 300

q1, q2, q3, q4, l2 = symbols("q1 q2 q3 q4 l2")

_0T1 = Cdh(np.array([0, 0]), l2, l1+l3, 0)
_1T2 = Cdh(np.array(["-q1", 0]), 0, l4, 0)
_2T3 = Cdh(np.array(["-q2", 0]), 0, l5, 0)
_3T4 = Cdh(np.array(["-q3", 90]), 0, 0, 90)
_4T5 = Cdh(np.array(["q4", 0]), l6, 0, 0)
_0T2 = _0T1.getTrans() * _1T2.getTrans()
_0T3 = _0T2 * _2T3.getTrans()
_0T4 = _0T3 * _3T4.getTrans()
_0T5 = _0T1.getTrans() * _1T2.getTrans() * _2T3.getTrans() * _3T4.getTrans() * _4T5.getTrans()

print("######0T1#########")
pprint(_0T1.getTrans())
print("######1T2######")
pprint(_1T2.getTrans())

print("######2T3######")
pprint(_2T3.getTrans())
print("######3T4######")
pprint(_3T4.getTrans())
print("######4T5######")
pprint(_4T5.getTrans())

print("######0T5######")
pprint(simplify(_0T5))

#########Achse1##########################
trans1, Or1 = calcJacobiTrans("", 1)
###########Achse2########################
trans2, Or2 = calcJacobiRot(_0T5, _0T1.getTrans(), Matrix([0,0,0,1]), 0)
#########Achse3########
trans3, Or3 = calcJacobiRot(_0T5, _0T2,  Matrix([0,0,0,1]), 0)
###########Achse4#########
trans4, Or4 = calcJacobiRot(_0T5, _0T3,  Matrix([0,0,0,1]), 0)
##########Achse5###########
trans5, Or5 = calcJacobiRot(_0T5, _0T4,  Matrix([0,0,0,1]), 0)

print("######Jacobi######")
Jacobi = Matrix([ [trans1, trans2, trans3, trans4, trans5],
                  [Or1, Or2, Or3, Or4, Or5]])
pprint(simplify(Jacobi))
Jacobi = Jacobi.subs([(l2,0), (q1,0), (q2,0), (q3,0),(q4,0)])
print("######Jacobi mit Nulllage######")
pprint(Jacobi)