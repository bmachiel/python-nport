from __future__ import division

import numpy as np
import nport

from nport import Z, Y, S, T, ABCD

class TwoNPortMatrix(nport.NPortMatrixBase):
    """Class representing a 2n-port matrix (Z, Y, S, T, G, H or ABCD)"""
    # See also "block matrices" in the Matrix Cookbook
    def __new__(cls, matrix, type, z0=None):
        matrix = np.asarray(matrix, dtype=complex)
        obj = nport.NPortMatrixBase.__new__(cls, matrix, type, z0)
        if len(obj.shape) != 4:
            raise ValueError("the matrix should be four-dimensional")
        if obj.shape[0] != 2 or obj.shape[0] != obj.shape[1]:
            raise ValueError("the matrices should be 2x2 matrices")
        if obj.shape[2] != obj.shape[3]:
            raise ValueError("the submatrices should be square")
        return obj

    @property
    def ports(self):
        return 2 * self.shape[2]

    def nportmatrix(self):
        matrix = np.vstack((np.hstack((self[0,0], self[0,1])),
            np.hstack((self[1,0], self[1,1]))))
        return nport.NPortMatrix(matrix, self.type, self.z0)

    def renormalize(self, z0):
        """Renormalize the 2n-port parameters to `z0`"""
        assert self.type in (S, T)
        if z0 == self.z0:
            result = self
        elif self.type == S:
            result = self.nportmatrix().renormalize(z0).twonportmatrix()
        elif self.type == T:
            result = self.convert(S).renormalize(z0).convert(T)
        return result

    def convert(self, type, z0=None):
        """Convert from one 2n-port matrix representation to another"""
        z0 = self.convert_z0test(type, z0)
        idty = np.identity(self.shape[2], dtype=complex)
        invert = np.linalg.inv
        
        # TODO: check for singularities
        if type == self.type:
            if type in (S, T):
                return self.renormalize(z0)
            else:
                return self                
        elif self.type == S and type == ABCD:
            s11 = np.asmatrix(self[0,0])
            s12 = np.asmatrix(self[0,1])
            s21 = np.asmatrix(self[1,0])
            s22 = np.asmatrix(self[1,1])
            s21i = np.linalg.inv(s21)
            i_plus_s11__times_s21i = (idty + s11) * s21i
            i_minus_s11__times_s21i = (idty - s11) * s21i
            i_plus_s22 = (idty + s22)
            i_minus_s22 = (idty - s22)
            a = (i_plus_s11__times_s21i * i_minus_s22 + s12)
            b = (i_plus_s11__times_s21i * i_plus_s22 - s12) * self.z0
            c = (i_minus_s11__times_s21i * i_minus_s22 - s12) / self.z0
            d = (i_minus_s11__times_s21i * i_plus_s22 + s12)
            result = 0.5 * self.__class__([[a, b], [c, d]], ABCD, z0)
        elif self.type == ABCD and type == S:
            # via scattering transfer parameters
            # TODO
            t = self.convert(T, z0)
            result = t.convert(S)
        elif self.type == S and type == T:
            # first renormalize
            if z0 != self.z0:
                result = self.renormalize(z0).convert(T, z0)
            else:
                s11 = np.asmatrix(self[0,0])
                s12 = np.asmatrix(self[0,1])
                s21 = np.asmatrix(self[1,0])
                s22 = np.asmatrix(self[1,1])
                s21i = np.linalg.inv(s21)
                t11 = s21i
                t12 = - s21i * s22
                t21 = s11 * s21i
                t22 = s12 - s11 * s21i * s22
                t = self.__class__([[t11, t12], [t21, t22]], T, z0)
                result = t.renormalize(z0)
        elif self.type == T and type == S:
            t11 = np.asmatrix(self[0,0])
            t12 = np.asmatrix(self[0,1])
            t21 = np.asmatrix(self[1,0])
            t22 = np.asmatrix(self[1,1])
            t11i = np.linalg.inv(t11)
            s11 = t21 * t11i
            s12 = t22 - t21 * t11i * t12
            s21 = t11i
            s22 = - t11i * t12
            result = self.__class__([[s11, s12], [s21, s22]], S, self.z0)
            if z0 != self.z0:
                result = result.renormalize(z0)
        elif self.type == T and type == ABCD:
            t11 = np.asmatrix(self[0,0])
            t12 = np.asmatrix(self[0,1])
            t21 = np.asmatrix(self[1,0])
            t22 = np.asmatrix(self[1,1])
            a = t11 + t12 + t21 + t22
            b = (t11 - t12 + t21 - t22) * self.z0
            c = (t11 + t12 - t21 -t22) / self.z0
            d = t11 - t12 -t21 + t22
            result = 0.5 * self.__class__([[a, b], [c, d]], ABCD, z0)
        elif self.type == ABCD and type == T:
            a = np.asmatrix(self[0,0])
            b0 = np.asmatrix(self[0,1]) / z0
            c0 = np.asmatrix(self[1,0]) * z0
            d = np.asmatrix(self[1,1])
            t11 = a + b0 + c0 + d
            t12 = a - b0 + c0 - d
            t21 = a + b0 - c0 - d
            t22 = a - b0 - c0 + d
            result = 0.5 * self.__class__([[t11, t12], [t21, t22]], T, z0)
        elif self.type == ABCD and type == Z:
            a = np.asmatrix(self[0,0])
            b = np.asmatrix(self[0,1])
            c = np.asmatrix(self[1,0])
            d = np.asmatrix(self[1,1])
            ci = np.linalg.inv(c)
            z11 = a * ci
            z12 = a * ci * d - b
            z21 = ci
            z22 = ci * d
            result = self.__class__([[z11, z12], [z21, z22]], Z, z0)
        elif self.type == Z and type == ABCD:
            z11 = np.asmatrix(self[0,0])
            z12 = np.asmatrix(self[0,1])
            z21 = np.asmatrix(self[1,0])
            z22 = np.asmatrix(self[1,1])
            z21i = np.linalg.inv(z21)
            a = z11 * z21i
            b = z11 * z21i * z22 - z12
            c = z21i
            d = z21i * z22
            result = self.__class__([[a, b], [c, d]], ABCD, z0)
        elif self.type == ABCD and type == Y:
            a = np.asmatrix(self[0,0])
            b = np.asmatrix(self[0,1])
            c = np.asmatrix(self[1,0])
            d = np.asmatrix(self[1,1])
            bi = np.linalg.inv(b)
            y11 = d * bi
            y12 = c - d * bi * a
            y21 = - bi
            y22 = bi * a
            result = self.__class__([[y11, y12], [y21, y22]], Y, z0)
        elif self.type == Y and type == ABCD:
            y11 = np.asmatrix(self[0,0])
            y12 = np.asmatrix(self[0,1])
            y21 = np.asmatrix(self[1,0])
            y22 = np.asmatrix(self[1,1])
            y21i = np.linalg.inv(y21)
            a = - y21i * y22
            b = - y21i
            c = y12 - y11 * y21i * y22
            d = - y11 * y21i
            result = self.__class__([[a, b], [c, d]], ABCD, z0)
        elif self.type == Z and type == T:
            z11 = np.asmatrix(self[0,0]) / z0
            z12 = np.asmatrix(self[0,1]) / z0
            z21 = np.asmatrix(self[1,0]) / z0
            z22 = np.asmatrix(self[1,1]) / z0
            z21i = np.linalg.inv(z21)
            i_plus_z11__times_z21i = (idty + z11) * z21i
            z11_minus_i__times_z21i = (z11 - idty) * z21i
            i_plus_z22 = (idty + z22)
            i_minus_z22 = (idty - z22)
            t11 = i_plus_z11__times_z21i * i_plus_z22 - z12
            t12 = i_plus_z11__times_z21i * i_minus_z22 + z12
            t21 = z11_minus_i__times_z21i * i_plus_z22 - z12
            t22 = z11_minus_i__times_z21i * i_minus_z22 + z12
            result = 0.5 * self.__class__([[t11, t12], [t21, t22]], T, z0)
        elif self.type == T and type == Z:
            # via scattering parameters
            s = self.convert(S).nportmatrix()
            z = s.convert(Z)
            return z.twonportmatrix()            
        elif self.type == Y and type == T:
            y11 = np.asmatrix(self[0,0]) * z0
            y12 = np.asmatrix(self[0,1]) * z0
            y21 = np.asmatrix(self[1,0]) * z0
            y22 = np.asmatrix(self[1,1]) * z0
            y21i = np.linalg.inv(y21)
            i_plus_y11__times_y21i = (idty + y11) * y21i
            i_minus_y11__times_y21i = (idty - y11) * y21i
            i_plus_y22 = (idty + y22)
            i_minus_y22 = (idty - y22)
            t11 = - i_plus_y11__times_y21i * i_plus_y22 + y12
            t12 = i_plus_y11__times_y21i * i_minus_y22 + y12
            t21 = - i_minus_y11__times_y21i * i_plus_y22 - y12
            t22 = i_minus_y11__times_y21i * i_minus_y22 - y12
            result = 0.5 * self.__class__([[t11, t12], [t21, t22]], T, z0)
        elif self.type == T and type == Y:
            # via scattering parameters
            s = self.convert(S).nportmatrix()
            y = s.convert(Y)
            return y.twonportmatrix()
        else:
            nportmatrix = self.nportmatrix()
            # catch infinite recursion for TwoPortMatrix
            nportmatrix.__class__ = nport.NPortMatrix
            result = nportmatrix.convert(type, z0)
            return result.twonportmatrix()

        return result


class TwoNPort(nport.NPortBase):
    """Class representing a 2n-port across a list of frequencies"""
    def __new__(cls, freqs, matrices, type, z0=None):
        matrices = np.asarray(matrices, dtype=complex)
        obj = nport.NPortBase.__new__(cls, freqs, matrices, type, z0)
        if len(obj[0].shape) != 4:
            raise ValueError("the matrices should be four-dimensional")
        if obj[0].shape[0] != 2 or obj[0].shape[0] != obj[0].shape[1]:
            raise ValueError("the matrices should be a 2x2 matrix")
        if obj[0].shape[2] != obj[0].shape[3]:
            raise ValueError("the submatrices should be square")
        return obj

    def __init__(self, freqs, matrices, type, z0=50):
        """
        Initialize an instance, specifying the frequency samples, all
        corresponsing matrices, the matrix type and the reference impedance

        """
        pass

    def __getitem__(self, arg):
        if type(arg) == int:
            return TwoNPortMatrix(np.asarray(self).__getitem__(arg), self.type, self.z0)
        else:
            return np.asarray(self).__getitem__(arg)

    @property
    def ports(self):
        """Return the number of ports of this TwoNPort"""
        return 2 * self[0].shape[2]

    def add(self, freq, matrix):
        """
        Return a TwoNPort with the specified frequency sample added

        If matrix is a TwoNPortMatrix, its elements will be converted to this
        TwoNPort's type and characteristic impedance.
        If matrix is a complex array, it is assumed the elements are in the
        format of this TwoNPort

        """
        if type(matrix) == TwoNPortMatrix:
            if matrix.type != self.type or matrix.z0 != self.z0:
                matrix = matrix.convert(self.type, self.z0)
        index = self.freqs.searchsorted(freq)
        freqs = np.insert(self.freqs, index, freq)
        matrices = np.insert(self.matrices, index, matrix, 0)
        return TwoNPort(freqs, matrices, self.type, self.z0)

    def nport(self):
        """Convert this TwoNPort to an NPort"""
        nportmatrices = []
        for matrix in self:
            twonportmatrix = TwoNPortMatrix(matrix, self.type, self.z0)
            nportmatrices.append(twonportmatrix.nportmatrix())
        return nport.NPort(self.freqs, nportmatrices, self.type, self.z0)

    def renormalize(self, z0):
        """Renormalize the 2n-port parameters to `z0`"""
        assert self.type in (S, T)
        if z0 == self.z0:
            result = self
        else:
            renormalized = []
            for matrix in self:
                twonportmatrix = TwoNPortMatrix(matrix, self.type, self.z0)
                renormalized.append(twonportmatrix.renormalize(z0))
            result = self.__class__(self.freqs, renormalized, self.type, z0)
        return result

    def convert(self, type, z0=None):
        """Convert to another 2n-port matrix representation"""
        converted = []
        for matrix in self:
            twonportmatrix = TwoNPortMatrix(matrix, self.type, self.z0)
            converted.append(twonportmatrix.convert(type, z0))
        return TwoNPort(self.freqs, converted, type, z0)

    #~ def invert(self):
        #~ """Return a TwoNPort described by the inverse of this TwoNPort's
        #~ matrices
        
        #~ """
        #~ inverted = []
        #~ for matrix in self:
            #~ nportmatrix = NPortMatrix(matrix, self.type, self.z0)
            #~ inverted.append(np.linalg.inv(nportmatrix))
        #~ return NPort(self.freqs, inverted, self.type, self.z0)
