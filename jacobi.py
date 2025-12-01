import classesDH as dh
import numpy as np

def calcJacobiRot(_0t5, _0ta,): #a=i-1
    TCPVec = np.array([0, 0, 0, 1])
    _0rE = _0t5 @ TCPVec
    _0rE=np.delete(_0rE,3, 0)

    _0ra = _0ta[:,3] #4 Spalte
    _0ra = np.delete(_0ra, 3, 0)
    print(_0ra)

    _0ei = _0ta[:,2] #3 Spalte
    _0ei = np.delete(_0ei, 3, 0)
    print(_0ei)

    trans = np.cross(_0ei, (_0rE-_0ra)) #geht nur mit 3 Zeilen 
    Or = _0ei
    print(trans)
    return trans, Or

l1 = 100
l3 = 600
l4 = 600
l5 = 600
l6 = 300

_0t1 = dh.Cdh_trans(0,    l1+l3,    0,      200,        2*1434,   0)
#                  _d,      _a,     _alpha,      _stepSize, _max, _min
_1t2 = dh.Cdh_rot(0,        l4,      0,              10,     -5,     0)
_2t3 = dh.Cdh_rot(0,        l5,      0,              10,     168,    0)
_3t4 = dh.Cdh_rot(0,         0,      90,             10,     125+90, 90) 
_4t5 = dh.Cdh_rot(l6,        0,      0,              0,       0,      0) #Endeffektor

_0t1.setZero()
_1t2.setZero()
_2t3.setZero()
_3t4.setZero()
_4t5.setZero()

_0t2 = _0t1.getTrans() @ _1t2.getTrans()
_0t3 = _0t2 @ _2t3.getTrans()
_0t4 = _0t3 @ _3t4.getTrans()
_0t5 = _0t4 @ _4t5.getTrans()

TCPVec = np.array([0, 0, 0, 1])
_0rE = _0t5 * TCPVec
_0rE=np.delete(_0rE,3, 0)
#########Achse1##########################
#_0e1 = _0T1.getTrans().col(2) #3 Spalte
#_0e1.row_del(3)
_0e1 = np.array([0, 0, 0]) #kann das Stimmen? Heißt das die Schubachse keine Einfluss hat?
trans1 = _0e1
Or1 = np.array([0, 0, 0])
###########Achse2########################
trans2, Or2 = calcJacobiRot(_0t5, _0t1.getTrans())
#########Achse3########
trans3, Or3 = calcJacobiRot(_0t5, _0t2)
###########Achse4#########
trans4, Or4 = calcJacobiRot(_0t5, _0t3)
##########Achse5###########
trans5, Or5 = calcJacobiRot(_0t5, _0t4)

JacobiTrans = np.column_stack([trans1, trans2, trans3, trans4, trans5])
JacobiRot = np.column_stack([Or1, Or2, Or3, Or4, Or5])
print(JacobiTrans)

Jacobi = np.vstack([ JacobiTrans, JacobiRot])
print(Jacobi)
np.set_printoptions(precision=4, suppress=True)
print("Jacobi-Matrix:\n", Jacobi)

eigvals, eigvecs = np.linalg.eig(Jacobi@Jacobi.T)
# Diagonalmatrix D mit Eigenwerten
D = np.diag(eigvals)
print("Eigenwerte D (Diag-Matrix):")
print(D)