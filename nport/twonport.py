from __future__ import division

import numpy as np
import nport


class TwoNPortMatrix(nport.NPortMatrixBase):
    """Class representing a 2n-port matrix (Z, Y, S, T, G, H or ABCD)"""
    # See also "block matrices" in the Matrix Cookbook
    def __new__(cls, matrix, mtype, z0=None):
        obj = nport.NPortMatrixBase.__new__(cls, matrix, mtype, z0)
        if len(obj.shape) != 4:
            raise ValueError("the matrix should be four-dimensional")
        if obj.shape[0] != 2 or obj.shape[0] != obj.shape[1]:
            raise ValueError("the matrix should be a 2x2 matrix")
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

    def convert(self, mtype, z0=None):
        """Convert from one 2n-port matrix representation to another"""
        z0 = nport.NPortMatrixBase.__z0test__(mtype, z0)
        identity = np.identity(self.shape[2], dtype=complex)

        # TODO: check for singularities
        if self.type == nport.SCATTERING and mtype == nport.TRANSMISSION:
                s11 = np.asmatrix(self[0,0])
                s12 = np.asmatrix(self[0,1])
                s21 = np.asmatrix(self[1,0])
                s22 = np.asmatrix(self[1,1])
                s21i = np.linalg.inv(s21)
                a = 0.5 * ((identity + s11) * s21i * (identity - s22) + s12)
                b = self.z0 / 2 * ((identity + s11) * s21i * (identity + s22) - s12)
                c = 0.5 / self.z0 * ((identity - s11) * s21i * (identity - s22) - s12)
                d = 0.5 * ((identity - s11) * s21i * (identity + s22) + s12)
                result = [[a, b], [c, d]]
        elif self.type == nport.SCATTERING and mtype == nport.SCATTERING_TRANSFER:
            raise NotImplementedError
        elif self.type == nport.TRANSMISSION and mtype == nport.SCATTERING:
                # see PhDwiki
                a = np.asmatrix(self[0,0])
                b0 = np.asmatrix(self[0,1]) / z0
                c0 = np.asmatrix(self[1,0]) * z0
                d = np.asmatrix(self[1,1])
                x = np.linalg.inv(a + b0 + c0 + d)
                s11 = (a + b0 - c0 -d) * x
                s12 = 0.5 * ((a - b0 - c0 + d) -
                    (a + b0 - c0 - d) * x * (a - b0 + c0 -d))
                s21 = 2 * x
                s22 = x * (-a + b0 - c0 + d)
                result = [[s11, s12], [s21, s22]]
                # alternative formulas [KHA06] p.21
        elif self.type == nport.TRANSMISSION and mtype == nport.SCATTERING_TRANSFER:
            raise NotImplementedError
        else:
            result = self.nportmatrix().convert(mtype, z0)
            return result.twonportmatrix()
            # ABCD -> Z and Y conversions in [PAU08] section 7.5.5

        return TwoNPortMatrix(result, mtype, z0)


class TwoNPort(nport.NPortBase):
    """Class representing a 2n-port across a list of frequencies"""
    def __new__(cls, freqs, matrices, mtype, z0=None):
        obj = nport.NPortBase.__new__(cls, freqs, matrices, mtype, z0)
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
