from sympy import symbols, Matrix, cos, sin, pprint, trigsimp, pi
import numpy as np

class Cdh:
    phi = np.empty(2)
    alpha = 0
    a = 0
    d = 0
    phi, d, a, alpha = symbols('phi d a alpha')

    def __init__(self, _phi:np.ndarray, _d, _a, _alpha): #array of phi and alph should be ["value or symbol", Offset[int]]
        self.alpha = _alpha
        self.a = _a
        self.d = _d
        self.phi = _phi
        self.calcTrans()

    def calcTrans(self):
        phi_offset = symbols("phi_offset")
        phi, phi_offset, alpha, a, d = self.phi[0], self.phi[1], self.alpha, self.a, self.d

        phi_offset = self._convertToRad(phi_offset)
        alpha = self._convertToRad(alpha)

        if not type(phi) == np.str_:
            phi = self._convertToRad(phi)

        cPhi = self._AdditionstheoremCOS(phi, phi_offset)
        sPhi = self._AdditionstheoremSIN(phi, phi_offset)
        cAlpha = cos(alpha)
        sAlpha = sin(alpha)
        self.tran = Matrix([
            [cPhi, -sPhi*cAlpha,  sPhi*sAlpha,  a*cPhi],
            [sPhi,  cPhi*cAlpha, -cPhi*sAlpha,  a*sPhi],
            [0,         sAlpha,           cAlpha,           d],
            [0,         0,                    0,                    1]
        ])
 
    def getTrans(self):
        return self.tran
    def _AdditionstheoremCOS(self, x1, x2):
        return cos(x1) * cos(x2) + -1 * sin(x1) * sin(x2)
    def _AdditionstheoremSIN(self, x1, x2):
        return sin(x1) * cos(x2) + cos(x1) * sin(x2)
    def _convertToRad(self, deg):
        #print(type(deg))
        if type(deg) == np.float16 or type(deg) == np.int64 or type(deg) == int:
            return (pi/180)*deg
        if type(deg) == np.str_:
            if deg.isdigit():
               return (pi/180)*float(deg) 
        print("Wert nicht convertierbar in RAD")
        exit()

def calcJacobiRot(_0t, _0ta, TCPVec, flac): #0t is whole Transformation; a=i-1; TCPVec= Matrix([0, 0, 0, 1]); if flac 1 --> first Element ez0 = [0,0,1] 0r0=[0,0,0]
    _0rE = _0t * TCPVec
    _0rE.row_del(3)
    if flac:
        _0r0 = Matrix([0,0,0])
        _0e0 = Matrix([0,0,1])
        trans = _0e0.cross(_0rE-_0r0) #geht nur mit 3 Zeilen 
        Or = _0e0
        return trans, Or

    _0ri = _0ta.col(3) #4 Spalte
    _0ri.row_del(3)
    #pprint(_0ri)
    _0ei = _0ta.col(2) #3 Spalte
    _0ei.row_del(3)

    #pprint(_0ei)
    trans = _0ei.cross(_0rE-_0ri) #geht nur mit 3 Zeilen 
    Or = _0ei
    #pprint(trans)
    return trans, Or
def calcJacobiTrans(_0ta, flac): #flac 1 --> first Element ez0 = [0,0,1]
    if flac:
        return Matrix([0,0,1]), Matrix([0,0,0])
    
    _0ei = _0ta.col(2) #3 Spalte
    _0ei.row_del(3)
    return _0ei, Matrix([0, 0, 0])
