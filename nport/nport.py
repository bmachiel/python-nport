
from __future__ import division

import numpy as np

from .base import Z, Y, S, T, H, G, ABCD
from .base import IMPEDANCE, ADMITTANCE, SCATTERING
from .base import NPortMatrixBase, NPortBase
from .parameter import rad


class NPortMatrix(NPortMatrixBase):
    """Class representing an n-port matrix (Z, Y or S; for 2-ports also T, G,
    H or ABCD)

    :param matrix: matrix elements
    :type matrix: *n* by *n* array
    :param type: matrix type
    :type type: :data:`Z`, :data:`Y`, :data:`S`, :data:`T`, :data:`G`, :data:`H`
                or :data:`ABCD`
    :param z0: normalizing impedance (only :class:`S` and :class:`T`)
    :type z0: :class:`float`
    
    """
    def __new__(cls, matrix, type, z0=None):
        matrix = np.asarray(matrix, dtype=complex)
        obj = NPortMatrixBase.__new__(cls, matrix, type, z0)
        if len(obj.shape) != 2:
            raise ValueError("the matrix should be two-dimensional")
        if obj.shape[0] != obj.shape[1]:
            raise ValueError("the matrix should be square")
        if obj.shape[0] == 2:
            obj.__class__ = TwoPortMatrix
        return obj

    @property
    def ports(self):
        """The number of ports of this :class:`NPortMatrix`
        
        :rtype: :class:`int`
        
        """
        return self.shape[0]

    def power(self, n):
        """Return this :class:`NPortMatrix` raised to the `n`\ th power
        
        :rtype: :class:`NPortMatrix`
        
        """
        return np.linalg.matrix_power(self, n)

    def twonportmatrix(self, inports=None, outports=None):
        """Return the 2n-port matrix represented by this n-port matrix
        
        :param inports: the list of ports that make up the inputs of the 2n-port
        :type inports: class:`tuple` or :class:`list`
        :param outports: the list of ports that make up the outputs
        :type outports: :class:`tuple` or :class:`list`
        :rtype: :class:`TwoNPortMatrix`
        
        """
        if self.ports % 2 != 0:
            raise TypeError("this NPortMatrix' number of ports is not a "
                            "multiple of 2")
        n = int(self.ports / 2)
        if inports is None and outports is None:
            inports = range(n)
            outports = range(n, 2*n)
            matrix = self
        else:
            # check whether the given sets of ports are valid
            assert inports is not None and outports is not None
            assert len(inports) == n
            assert len(outports) == n
            allports = set(inports).union(set(outports))
            assert len(allports) == 2*n
            assert min(allports) == 1 and max(allports) == 2*n
            inports = [port - 1 for port in inports]
            outports = [port - 1 for port in outports]

            ports = inports + outports

            # shuffle rows and columns to obtain 2n-port format
            matrix = []
            for row in ports:
                matrix.append(self[row])
            matrix = np.asarray(matrix).T
            matrix2 = []
            for column in ports:
                matrix2.append(matrix[column])
            matrix = np.asarray(matrix2).T

        matrix = np.asmatrix(matrix)
        x11 = matrix[0:n  , 0:n  ]
        x12 = matrix[0:n  , n:2*n]
        x21 = matrix[n:2*n, 0:n  ]
        x22 = matrix[n:2*n, n:2*n]
        matrix = np.array([[x11, x12], [x21, x22]])
        return TwoNPortMatrix(matrix, self.type, self.z0)

    def reverse(self):
        """Return an :class:`NPortMatrix` with the ports reversed. This simply
        flips the matrix horizontally and vertically.
        
        :rtype: :class:`NPortMatrix`
        
        """
        return self[::-1, ::-1]

    def renormalize(self, z0):
        """Renormalize the n-port parameters to `z0`
        
        :param z0: new normalizing impedance
        :type z0: :class:`float`
        :rtype: :class:`NPortMatrix`
        
        """
        # http://qucs.sourceforge.net/tech/node98.html
        # "Renormalization of S-parameters to different port impedances"
        if self.type not in (S, T):
            raise TypeError("Only S and T matrices can be renormalized")

        if z0 == self.z0:
            result = self
        elif self.type == S:
            idty = np.identity(len(self), dtype=complex)

            r = (z0 - self.z0) / (z0 + self.z0)
            result = self.__class__(np.dot(self - idty * r,
                                           np.linalg.inv(idty - r * self)),
                                    self.type, z0)
        elif self.type == T:
            result = self.convert(S, z0).convert(T)

        return result

    def convert(self, type, z0=None):
        """Convert to another n-port matrix representation
        
        :param type: new matrix type
        :type type: :class:`Z`, :class:`Y` or :class:`S`
        :param z0: normalizing impedance (only :data:`S` and :data:`T`)
        :type z0: :class:`float`
        :rtype: :class:`NPortMatrix`

        """
        # references:
        #  * http://qucs.sourceforge.net/tech/node98.html
        #           "Transformations of n-Port matrices"
        #  * Multiport S-Parameter Measurements of Linear Circuits With Open
        #           Ports by Reimann et al. (only S <-> Z)
        #  * Agilent AN 154 - S-Parameter Design (S <-> T)
        #  * MATLAB S-parameter toolbox (Z, Y, H, G, ABCD, T)
        #           http://www.mathworks.com/matlabcentral/fileexchange/6080
        z0 = self.convert_z0test(type, z0)
        idty = np.identity(len(self), dtype=complex)
        invert = np.linalg.inv

        if type in (ABCD, T):
            raise TypeError("Cannot convert an NPort to %s-parameter "
                            "representation. Convert to a TwoNPort first" %
                            type)
        elif type in (H, G):
            if self.ports != 2:
                raise TypeError("Can only convert 2x2 matrices to %s-parameter "
                                "representation" % type)
        elif type not in (Z, Y, S):
            raise TypeError("Unknown n-port parameter type")

        # TODO: check for singularities
        if self.type == SCATTERING:
            if type == SCATTERING:
                return self.renormalize(z0)
            elif type == IMPEDANCE:
                result = (2 * invert(idty - self) - idty) * self.z0
            elif type == ADMITTANCE:
                result = (2 * invert(idty + self) - idty) / self.z0
        elif self.type == IMPEDANCE:
            if type == SCATTERING:
                result = idty - 2 * invert(idty + self / z0)
            elif type == ADMITTANCE:
                result = invert(self)
            elif type == IMPEDANCE:
                result = self
        elif self.type == ADMITTANCE:
            if type == SCATTERING:
                result = 2 * invert(idty + self * z0) - idty
            elif type == IMPEDANCE:
                result = invert(self)
            elif type == ADMITTANCE:
                result = self

        return NPortMatrix(result, type, z0)

    def submatrix(self, ports):
        """Keep only the parameters corresponding to the given ports,
        discarding the others. For a Z-matrix and a Y-matrix this corresponds to
        terminating the discarded ports in an open or short circuit
        respectively.

        :param ports: list of ports to keep
        :type ports: iterable
        :rtype: :class:`NPortMatrix`
        
        """
        indices = [port - 1 for port in ports]
        submatrix = self[indices, :][:, indices]
        return NPortMatrix(submatrix, self.type, self.z0)

    def recombine(self, portsets):
        """Recombine ports, reducing the number of ports of this NPortMatrix.
        
        :param portsets: an :class:`int` specifies the number of a port that
                         needs to be kept as-is. If the :class:`int` is
                         negative, the port's polarity will be reversed. A
                         :class:`tuple` specifies a pair of ports that are to be
                         recombined into a single port. The second element of
                         this :class:`tuple` acts as the ground reference to the
                         first element.
        :type portsets: iterable of :class:`int`\s and :class:`tuple`\s
        :rtype: :class:`NPortMatrix`
        
        >>> recombine([(1,3), (2,4), 5, -6]
        
        will generate a four-port where:
        
        * port 1 is original port 1 referenced to port 3
        * port 2 is original port 2 referenced to port 4
        * port 3 is original port 5
        * port 4 is original port 6, but with reversed polarity

        """
        if self.type == IMPEDANCE:
            m = []
            for i, ports in enumerate(portsets):
                row = [0 for port in range(self.ports)]
                try:
                    if isinstance(ports, tuple):
                        assert len(ports) == 2
                        row[ports[0] - 1] = 1
                        row[ports[1] - 1] = -1
                    else:
                        assert isinstance(ports, int)
                        assert ports != 0
                        if ports > 0:
                            row[ports - 1] = 1
                        else:
                            row[-ports - 1] = -1
                except IndexError:
                    raise IndexError("specified port number is higher than "
                                     "number of ports")
                m.append(row)
            m = np.matrix(m, dtype=float)
            result = m * np.asarray(self) * m.T
            return self.__class__(result, self.type, self.z0)
        else:
            z_recombined = self.convert(IMPEDANCE).recombine(portsets)
            return z_recombined.convert(self.type)

    def shunt(self, portsets):
        """Connect ports together, reducing the number of ports of this
        :class:`NPortMatrix`.
        
        :param portsets: an :class:`int` specifies the number of a port that
                         needs to be kept as-is. A :class:`tuple` specifies a
                         set of ports that are to be connected together.
        :type portsets: iterable of :class:`int`\s and :class:`tuple`\s
        :rtype: :class:`NPortMatrix`

        >>> shunt([1, (2, 3), (4, 5, 6)]
        
        will generate a three-port where:
        
        * port 1 is original port 1
        * port 2 is original port 2 and 3 connected together
        * port 3 is original ports 4, 5 and 6 connected together

        """
        if self.type == ADMITTANCE:
            # columns
            tmp = np.zeros((self.ports, len(portsets)), dtype=complex)
            for i, ports in enumerate(portsets):
                try:
                    if isinstance(ports, tuple):
                        ports = [port - 1 for port in ports]
                        tmp[:, i] = np.sum(self[:, ports], 1)
                    else:
                        assert isinstance(ports, int)
                        assert ports > 0
                        tmp[:, i] = self[:, ports - 1]
                except IndexError:
                    raise IndexError("specified port number is higher than "
                                     "number of ports")
            # rows
            result = np.zeros((len(portsets), len(portsets)), dtype=complex)
            for i, ports in enumerate(portsets):
                try:
                    if isinstance(ports, tuple):
                        ports = [port - 1 for port in ports]
                        result[i, :] = np.sum(tmp[ports, :], 0)
                    else:
                        assert isinstance(ports, int)
                        assert ports > 0
                        result[i, :] = tmp[ports - 1, :]
                except IndexError:
                    raise IndexError("specified port number is higher than "
                                     "number of ports")
            return self.__class__(result, self.type, self.z0)
        else:
            y_shunted = self.convert(ADMITTANCE).shunt(portsets)
            return y_shunted.convert(self.type)

    def is_passive(self):
        """Check whether this n-port matrix is passive
        
        :rtype: :class:`bool`
        
        """
        if self.type != SCATTERING:
            return self.convert(SCATTERING).ispassive()
        else:
            return np.max(np.sum(np.abs(np.asarray(self))**2, 1)) <= 1

    def is_reciprocal(self):
        """Check whether this n-port matrix is reciprocal
        
        :rtype: :class:`bool`
        
        """
        raise NotImplementedError

    def is_symmetrical(self):
        """Check whether this n-port matrix is symmetrical
        
        :rtype: :class:`bool`
        
        """
        raise NotImplementedError


class NPort(NPortBase):
    """Class representing an n-port across a list of frequencies
    
    :param freqs: list of frequency samples
    :type freqs: :class:`list`
    :param matrix: list of matrix elements
    :type matrix: list of *n* by *n* arrays
    :param type: matrix type
    :type type: :data:`Z`, :data:`Y`, :data:`S`, :data:`T`, :data:`G`, :data:`H`
                or :data:`ABCD`
    :param z0: normalizing impedance (only :class:`S` and :class:`T`)
    :type z0: :class:`float`

    """
    matrix_cls = NPortMatrix
    
    def __new__(cls, freqs, matrices, type, z0=None):
        matrices = np.asarray(matrices, dtype=complex)
        obj = NPortBase.__new__(cls, freqs, matrices, type, z0)
        if len(obj[0].shape) != 2:
            raise ValueError("the matrices should be two-dimensional")
        if obj[0].shape[0] == 2:
            obj.__class__ = TwoPort
        return obj

    @property
    def ports(self):
        """The number of ports of this :class:`NPort`
        
        :rtype: :class:`int`
        
        """
        return self[0].ports

    def power(self, n):
        """Return this :class:`NPort`\'s matrices raised to the `n`\ th power
        
        :rtype: :class:`NPort`
        
        """
        matrices = [np.linalg.matrix_power(m, n) for m in self]
        return self.__class__(self.freqs, matrices, self.type, self.z0)

    def add(self, freq, matrix):
        """Return an NPort with the specified frequency sample added.

        :param freq: frequency at which to insert `matrix`
        :type freq: :class:`float`
        :param matrix: matrix to insert at `freq`
        :type matrix: :class:`NPortMatrix` or a complex array
        :rtype: :class:`NPort`

        If matrix is an :class:`NPortMatrix`, its elements will be converted to
        this :class:`NPort`'s type and characteristic impedance. If `matrix` is
        a complex array, it is assumed the elements are in the format of this
        :class:`NPort`.

        """
        if type(matrix) == NPortMatrix:
            if matrix.type != self.type or matrix.z0 != self.z0:
                matrix = matrix.convert(self.type, self.z0)
        index = self.freqs.searchsorted(freq)
        freqs = np.insert(self.freqs, index, freq)
        matrices = np.insert(self.matrices, index, matrix, 0)
        return self.__class__(freqs, matrices, self.type, self.z0)

    def twonport(self, inports=None, outports=None):
        """Convert this :class:`NPort` to a :class:`TwoNPort` using `inports`
        as the input ports and `outports` as the output ports.
        
        :param inports: the list of ports that make up the inputs of the 2n-port
        :type inports: :class:`tuple` or :class:`list`
        :param outports: the list of ports that make up the outputs
        :type outports: :class:`tuple` or :class:`list`
        :rtype: :class:`TwoNPort`

        """
        twonportmatrices = [m.twonportmatrix(inports, outports) for m in self]
        return TwoNPort(self.freqs, twonportmatrices, self.type, self.z0)

    def renormalize(self, z0):
        """Renormalize the n-port parameters to `z0`

        :param z0: new normalizing impedance
        :type z0: :class:`float`
        :rtype: :class:`NPort`

        """
        if self.type not in (S, T):
            raise TypeError("Only S and T matrices can be renormalized")

        if z0 == self.z0:
            result = self
        else:
            renormalized = [matrix.renormalize(z0) for matrix in self]
            result = self.__class__(self.freqs, renormalized, self.type, z0)
        return result

    def reverse(self):
        """Return an :class:`NPort` with the ports reversed. This simply flips
        the matrices horizontally and vertically.
        
        :rtype: :class:`NPort`
        
        """
        reversed = self[:, ::-1, ::-1]
        return NPort(self.freqs, reversed, self.type, self.z0)

    def convert(self, type, z0=None):
        """Convert to another n-port matrix representation
        
        :param type: n-port representation type to convert to
        :type type: :data:`Z`, :data:`Y`, :data:`S`, :data:`T`, :data:`G`,
                    :data:`H` or :data:`ABCD`
        :param z0: normalizing impedance (only :data:`S` and :data:`T`)
        :type z0: :class:`float`

        """
        converted = [matrix.convert(type, z0) for matrix in self]
        return NPort(self.freqs, converted, type, z0)

    def submatrix(self, ports):
        """Keep only the parameters corresponding to the given ports, discarding
        the others. For a Z-matrix and a Y-matrix this corresponds to
        terminating the discarded ports in an open or short circuit
        respectively.
        
        :param ports: list of ports to keep
        :type ports: iterable
        
        """
        indices = [port - 1 for port in ports]
        submatrices = self[:, indices, :][:, :, indices]
        return NPort(self.freqs, submatrices, self.type, self.z0)

    def invert(self):
        """Return an :class:`NPort` described by the inverse of this
        :class:`NPort`'s matrices
        
        :rtype: :class:`NPort`

        """
        # TODO: determine type of output
        inverted = [np.linalg.inv(matrix) for matrix in self]
        return NPort(self.freqs, inverted, self.type, self.z0)

    def recombine(self, portsets):
        """Recombine ports, reducing the number of ports of this :class:`NPort`.
        
        :param portsets: an :class:`int` specifies the number of a port that
                         needs to be kept as-is. If the :class:`int` is
                         negative, the port's polarity will be reversed. A
                         :class:`tuple` specifies a pair of ports that are to be
                         recombined into a single port. The second element of
                         this :class:`tuple` acts as the ground reference to the
                         first element.
        :type portsets: iterable of :class:`int`\s and :class:`tuple`\s
        :rtype: :class:`NPort`

        >>> recombine([(1,3), (2,4), 5, -6]
        
        will generate a four-port where:
        
        * port 1 is original port 1 referenced to port 3
        * port 2 is original port 2 referenced to port 4
        * port 3 is original port 5
        * port 4 is original port 6, but with reversed polarity

        """
        recombined = [matrix.recombine(portsets) for matrix in self]
        return self.__class__(self.freqs, recombined, self.type, self.z0)

    def shunt(self, portsets):
        """Connect ports together, reducing the number of ports of this
        :class:`NPort`.
        
        :param portsets: an :class:`int` specifies the number of a port that
                         needs to be kept as-is. A :class:`tuple` specifies a
                         set of ports that are to be connected together.
        :type portsets: iterable of :class:`int`\s and :class:`tuple`\s
        :rtype: :class:`NPort`

        >>> shunt([1, (2, 3), (4, 5, 6)]
        
        will generate a three-port where:
        
        * port 1 is original port 1
        * port 2 is original port 2 and 3 connected together
        * port 3 is original ports 4, 5 and 6 connected together

        """
        shunted = [matrix.shunt(portsets) for  matrix in self]
        return self.__class__(self.freqs, shunted, self.type, self.z0)

    def is_passive(self):
        """Check whether this :class:`NPort` is passive
        
        :rtype: :class:`bool`
        
        """
        for matrix in self:
            if not matrix.is_passive():
                return False
        return True

    def group_delay(self, port1, port2):
        """Return the group delay of the parameter as specified by the indices
        `port1` and `port2`
        
        :param port1: index of input port
        :type port1: :class:`int`
        :param port2: index of output port
        :type port2: :class:`int`
        :returns: group delay of parameter at indices `port1` and `port2`
        :rtype: :class:`ndarray`
        
        """
        if self.type not in (S, ):
            raise TypeError("Group delay only makes sense for S-parameters (?)")
        phase = np.unwrap(rad(self.get_parameter(port1, port2)))
        dphase = np.gradient(phase)
        dfreq = np.gradient(self.freqs)
        return - dphase / (2 * np.pi * dfreq)


def dot(arg1, arg2):
    """Matrix multiplication for :class:`NPort`s and :class:`TwoNPort`s
    
    :type left: :class:`NPort`, :class:`TwoNPort`, :class:`NPortMatrix`,
                :class:`TwoNPortMatrix` or :class:`ndarray`
    :type right: :class:`NPort`, :class:`TwoNPort`, :class:`NPortMatrix`,
                 :class:`TwoNPortMatrix` or :class:`ndarray`
    
    """
    def merge_freqs(freqs1, freqs2):
        minf = max(freqs1[0], freqs2[0])
        maxf = min(freqs1[-1], freqs2[-1])
        result = list(set(freqs1).union(set(freqs2)))
        result.sort()
        return np.array(result[result.index(minf):result.index(maxf)])
    
    if isinstance(arg1, NPort):
        if isinstance(arg2, NPort):
            result_freqs = merge_freqs(arg1.freqs, arg2.freqs)
            arg1_matrices = arg1.at(result_freqs)
            arg2_matrices = arg2.at(result_freqs)
            result_matrices = np.asarray([np.dot(a, b) for (a, b) in
                zip(arg1_matrices, arg2_matrices)])
        else:
            result_freqs = arg1.freqs
            result_matrices = np.array([np.dot(matrix, other)
                for matrix in arg1.matrices])
        return NPort(result_freqs, result_matrices, arg1.type, arg1.z0)
    elif isinstance(arg1, TwoNPort):
        if isinstance(arg2, TwoNPort):
            def twonport_dot(l, r):
                a = np.dot(l[0, 0], r[0, 0]) + np.dot(l[0, 1], r[1, 0])
                b = np.dot(l[0, 0], r[0, 1]) + np.dot(l[0, 1], r[1, 1])
                c = np.dot(l[1, 0], r[0, 0]) + np.dot(l[1, 1], r[1, 0])
                d = np.dot(l[1, 0], r[0, 1]) + np.dot(l[1, 1], r[1, 1])
                return np.asarray([[a, b], [c, d]])
            
            result_freqs = merge_freqs(arg1.freqs, arg2.freqs)
            arg1_matrices = arg1.at(result_freqs)
            arg2_matrices = arg2.at(result_freqs)
            result_matrices = np.asarray([twonport_dot(a, b) for (a, b) in
                zip(arg1_matrices, arg2_matrices)])
        else:
            raise NotImplementedError
        return TwoNPort(result_freqs, result_matrices, arg1.type,
                                 arg1.z0)
    elif type(arg2) == NPort:
        raise NotImplementedError
    else:
        return np.dot(arg1, arg2)


from .twoport import TwoPortMatrix, TwoPort
from .twonport import TwoNPortMatrix, TwoNPort
