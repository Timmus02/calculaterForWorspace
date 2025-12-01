from sympy import symbols, Matrix, cos, sin, pprint
import numpy as np
from sympyCalcTransformationMatrix import Cdh

def calcJacobiRot(_0t5, _0ta,): #a=i-1
    TCPVec = Matrix([0, 0, 0, 1])
    _0rE = _0t5 * TCPVec
    _0rE.row_del(3)

    _0ri = _0ta.col(3) #4 Spalte
    _0ri.row_del(3)
    pprint(_0ri)
    _0ei = _0ta.col(2) #3 Spalte
    _0ei.row_del(3)

    pprint(_0ei)
    trans = _0ei.cross(_0rE-_0ri) #geht nur mit 3 Zeilen 
    Or = _0ei
    pprint(trans)
    return trans, Or

l2 = "l2"

l1 = 100
l3 = 600
l4 = 600
l5 = 600
l6 = 300

_0T1 = Cdh(0, l2, l1+l3, 0)
_1T2 = Cdh("-q1", 0, l4, 0)
_2T3 = Cdh("-q2", 0, l5, 0)
_3T4 = Cdh("-q3+90", 0, 0, 90)
_4T5 = Cdh("q4", l6, 0, 0)
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
pprint(_0T5)

##calc 0rE
TCPVec = Matrix([0, 0, 0, 1])
_0rE = _0T5 * TCPVec
_0rE.row_del(3)
#########Achse1##########################
#_0e1 = _0T1.getTrans().col(2) #3 Spalte
#_0e1.row_del(3)
_0e1 = Matrix([0, 0, 0]) #kann das Stimmen? Heißt das die Schubachse keine Einfluss hat?
trans1 = _0e1
Or1 = Matrix([0, 0, 0])
###########Achse2########################
trans2, Or2 = calcJacobiRot(_0T5, _0T1.getTrans())
#########Achse3########
trans3, Or3 = calcJacobiRot(_0T5, _0T2)
###########Achse4#########
trans4, Or4 = calcJacobiRot(_0T5, _0T3)
##########Achse5###########
trans5, Or5 = calcJacobiRot(_0T5, _0T4)

Jacobi = Matrix([ [trans1, trans2, trans3, trans4, trans5],
                  [Or1, Or2, Or3, Or4, Or5]])
Jacobi.simplify
print("Jacobi")
pprint(Jacobi)
print(Jacobi.shape)
