
from __future__ import division

import numpy as np
from scipy import interpolate
import parameter
import twonport


# parameter matrix types
Z = IMPEDANCE = TYPE_Z = 'Z'
Y = ADMITTANCE = TYPE_Y = 'Y'
S = SCATTERING = TYPE_S = 'S'
T = SCATTERING_TRANSFER = TYPE_T = 'T'
H = HYBRID = TYPE_H = 'H'
G = INVERSE_HYBRID = TYPE_G = 'G'
ABCD = TRANSMISSION = TYPE_ABCD = 'ABCD'


#~ class MatrixTypeError(Exception):
    #~ """Inappropriate parameter matrix type"""
    #~ def __init__(self, given, expected):
        #~ self.given = given
        #~ self.expected = expected
    #~ def __str__(self):
        #~ return "given type %s does not match expected type %s" % (self.given, self.expected)


class NPortMatrixBase(np.ndarray):
    """Base class representing an n-port matrix (Z, Y, S, T, G, H or ABCD)
    
    :type matrix: *n* by *n* array
    :param type: matrix type
    :type type: :data:`Z`, :data:`Y`, :data:`S`, :data:`T`, :data:`G`,
                :data:`H` or :data:`ABCD`
    :param z0: normalizing impedance (only S)
    :type z0: :class:`float`

    """
    # TODO: implement type checking for operations (as in NPortBase)
    def __new__(cls, matrix, type, z0=None):
        # http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
        obj = matrix.view(cls)
        if type not in (Z, Y, S, T, H, G, ABCD):
            raise ValueError("illegal matrix type specified")
        obj.type = type
        obj.z0 = obj._z0test(type, z0)
        return obj

    def _z0test(self, type, z0):
        """Verify whether a normalizing impedance may be specified for the given
        n-port parameter type
        
        :param type: n-port parameter type
        :type type: :data:`Z`, :data:`Y`, :data:`S`, :data:`T`, :data:`G`,
                    :data:`H` or :data:`ABCD`
        :param z0: normalizing impedance
        :type z0: :class:`float`
        
        """
        if type in (S, T):
            if z0 is None:
                z0 = 50.0
        elif z0 is not None:
            raise ValueError("the specified n-port representation (%s) does "
                             "not require a reference impedance" % type)
        return z0

    def __repr__(self):
        array_repr = repr(self.view(np.ndarray))[6:-1].replace('\n     ', '\n')
        attrs = 'type=%s' % self.type
        if self.type in (S, T):
            attrs += ', z0=%g' % self.z0
        return '%s(%s)\n %s' % (self.__class__.__name__,
                                attrs,
                                array_repr)
    
    def __array_finalize__(self, obj):
        if obj is None: return
        self.type = getattr(obj, 'type', None)
        self.z0 = getattr(obj, 'z0', None)

    @property
    def ports(self):
        """The number of ports of this :class:`NPortMatrixBase`
        
        :rtype: :class:`int`
        
        """
        raise NotImplementedError
    
    def convert_z0test(self, type, z0):
        """Check supplied `z0` for conversion and set default value depending on
        the type
        
        :param type: target type for conversion
        :type type: :data:`Z`, :data:`Y`, :data:`S`, :data:`T`, :data:`G`,
                    :data:`H` or :data:`ABCD`
        :param z0: target normalization impedance (only :data:`S` and :data:`T`)
        :type z0: :class:`float`
        :returns: `z0` or default value when `z0 == None`
        :rtype: :class:`float`
        
        """
        if type in (S, T):
            if z0 is None:
                if self.type in (S, T):
                    z0 = self.z0
                else:
                    z0 = 50.0
        elif z0 is not None:
            raise ValueError("the specified n-port representation (%s) does"
                             " not require a reference impedance" % type)
        return z0


class NPortBase(NPortMatrixBase):
    """Base class representing an n-port across a list of frequencies

    :param freqs: list of frequency samples
    :type freqs: :class:`list`
    :param matrices: list of matrix elements
    :type matrix: list of 2 by 2 arrays in which each element is an 
                  *n* by *n* array
    :param type: matrix type
    :type type: :data:`Z`, :data:`Y`, :data:`S`, :data:`T`, :data:`G`, :data:`H`
                or :data:`ABCD`
    :param z0: normalizing impedance (only :data:`S` and :data:`T`)
    :type z0: :class:`float`

    """
    def __metaclass__(cname, cbases, cdict):
        # http://thread.gmane.org/gmane.comp.python.general/373788/focus=373795
        P = type(cname, cbases, cdict)
        def make_op(name, func):
            def op(this, other):
                if issubclass(type(other), NPortBase):
                    if other.type == this.type:
                        if other.z0 != this.z0:
                            raise ValueError("operands have different Z0")
                        result_type = this.type
                        result_z0 = this.z0
                    else:
                        raise ValueError("operands have different types")
                    min_freq = max(this.freqs[0], other.freqs[0])
                    max_freq = min(this.freqs[-1], other.freqs[-1])
                    result_freqs = list(set(this.freqs).union(set(other.freqs)))
                    result_freqs.sort()
                    result_freqs = \
                        np.array(result_freqs[result_freqs.index(min_freq):
                                              result_freqs.index(max_freq)+1])
                    this_matrices = this.at(result_freqs)
                    other_matrices = other.at(result_freqs)
                    result_matrices = func(this_matrices, other_matrices)
                else:
                    result_freqs = this.freqs
                    result_matrices = func(this, other)
                    result_type = this.type
                    result_z0 = this.z0
                subclass = type(this)
                return subclass(result_freqs, result_matrices, result_type,
                                result_z0)
            op.__name__ = name
            return op
        for name, func in np.ndarray.__dict__.items():
            if callable(func) and name in ['__add__', '__sub__', '__mul__',
                '__div__', '__radd__', '__rsub__', '__rmul__', '__rdiv__']:
                setattr(P, name, make_op(name, func))
        return P

    def __new__(cls, freqs, matrices, type, z0=None):
        if len(freqs) != len(matrices):
            raise ValueError("the list of frequencies and the list of "
                             "matrices should have equal lenghts")
        if matrices[0].shape[0] != matrices[0].shape[1]:
            raise ValueError("the matrices should be square")
        obj = NPortMatrixBase.__new__(cls, matrices, type, z0)
        obj.freqs = np.asarray(freqs)
        return obj

    def __array_finalize__(self, obj):
        NPortMatrixBase.__array_finalize__(self, obj)
        self.freqs = getattr(obj, 'freqs', None)

    def __getitem__(self, arg):
        if type(arg) == int:
            return NPortMatrix(np.asarray(self).__getitem__(arg), self.type,
                               self.z0)
        else:
            return np.asarray(self).__getitem__(arg)

    #~ def __repr__(self):
        #~ return "freqs: " + repr(self.freqs) + "\n" + NPortMatrixBase.__repr__(self)

    #~ def __str__(self):
        #~ return self.__repr__()

    def get_parameter(self, port1, port2):
        """Return the parameter as specified by the indices `port1` and `port2`
        as an ndarray
        
        :type port1: :class:`int`
        :type port2: :class:`int`
        :returns: parameter at indices `port1` and `port2`
        :rtype: :class:`ndarray`
        
        """
        return np.asarray(self[:, port1 - 1, port2 - 1])

    def get_element(self, port1, port2):
        """Return the submatrices made up of the element as specified by the
        indices `port1` and `port2`
        
        """
        subclass = type(self)
        return subclass(self.freqs,
            np.asarray(self[:, (port1 - 1, ), :][:, :, (port2 - 1, )]),
            self.type, self.z0)

    def at(self, freqs):
        """Return the interpolated n-port data at the given
        * list of frequencies (`freqs` is iterable), or
        * at a single frequency (`freqs` is a value)
        
        :param freqs: frequency point or list of frequencies at which to return
                      the n-port matrices
        :type freqs: :class:`float` or iterable of :class:`float`s
        :rtype: `type(self)`
        
        """
        # check whether freqs is iterable
        try:
            it = iter(freqs)
            single = False
        except TypeError:
            single = True
            freqs = [freqs]
        subclass = type(self)
        func = interpolate.interp1d(self.freqs, self, axis=0)
        interpolated = func(freqs)
        interpolated_nport = subclass(freqs, interpolated, self.type, self.z0)
        if single:
            return interpolated_nport[0]
        else:
            return interpolated_nport
            
    def average(self, n):
        """Take a moving average over `n` frequency samples
        
        :param n: number of samples to average over
        :type n: :class:`int`
        :rtype: `type(self)`
        
        """
        averaged = np.zeros_like(np.asarray(self))
        for i in range(len(self)):
            for j in range(- int(np.floor(n/2)), int(np.ceil(n/2))):
                index = i + j
                if index < 0:
                    index = 0
                elif index >= len(self):
                    index = len(self) - 1
                averaged[i] += self[index] / n        
        # numpy convolve
        #~ ones = np.ones(n) / n
        #~ averaged = np.zeros_like(np.asarray(self))
        #~ for i in range(self.ports):
            #~ for j in range(self.ports):
                #~ averaged[:, i, j] = np.convolve(ones, self[:, i, j])
        # or scipy convolve
        #~ averaged = scipy.signal.convolve(np.ones((n, ) + self.shape[1:])/n, np.asarray(self), 'same')
        subclass = type(self)
        return subclass(self.freqs, averaged, self.type, self.z0)


class NPortMatrix(NPortMatrixBase):
    """Class representing an n-port matrix (Z, Y or S; for 2-ports also T, G,
    H or ABCD)

    :param matrix: matrix elements
    :type matrix: *n* by *n* array
    :param type: matrix type
    :type type: :data:`Z`, :data:`Y`, :data:`S`, :data:`T`, :data:`G`, :data:`H`
                or :data:`ABCD`
    :param z0: normalizing impedance (only :class:`S`)
    :type z0: :class:`float`
    
    """
    def __new__(cls, matrix, type, z0=None):
        matrix = np.asarray(matrix, dtype=complex)
        obj = NPortMatrixBase.__new__(cls, matrix, type, z0)
        if len(obj.shape) != 2:
            raise ValueError("the matrix should be two-dimensional")
        if obj.shape[0] != obj.shape[1]:
            raise ValueError("the matrix should be square")
        if matrix.shape[0] == 2:
            obj.__class__ = TwoPortMatrix
        return obj

    @property
    def ports(self):
        """The number of ports in this NPortMatrix
        
        :rtype: :class:`int`
        
        """
        return self.shape[0]

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
        return twonport.TwoNPortMatrix(matrix, self.type, self.z0)

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
        discarding the others

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

    def ispassive(self):
        """Check whether this n-port matrix is passive
        
        :rtype: :class:`bool`
        
        """
        if self.type != SCATTERING:
            return self.convert(SCATTERING).ispassive()
        else:
            if np.max(np.sum(np.abs(np.asarray(self))**2, 1)) > 1:
                return False
            else:
                return True

    def isreciprocal(self):
        """Check whether this n-port matrix is reciprocal
        
        :rtype: :class:`bool`
        
        """
        raise NotImplementedError

    def issymmetrical(self):
        """Check whether this n-port matrix is symmetrical
        
        :rtype: :class:`bool`
        
        """
        raise NotImplementedError


class TwoPortMatrix(NPortMatrix):
    """Class representing a 2-port matrix (Z, Y, S, T, G, H or ABCD)
    
    :param matrix: matrix elements
    :type matrix: 2 by 2 array
    :param type: matrix type
    :type type: :data:`Z`, :data:`Y`, :data:`S`, :data:`T`, :data:`G`, :data:`H`
                or :data:`ABCD`
    :param z0: normalizing impedance (only :data:`S` and :data:`T`)
    :type z0: :class:`float`

    """
    def convert(self, type, z0=None):
        """Convert to another 2-port matrix representation
        
        :param type: new matrix type
        :type type: :data:`Z`, :data:`Y`, :data:`S`, :data:`T`, :data:`G`,
                    :data:`H` or :data:`ABCD`
        :param z0: normalizing impedance (only :data:`S` and :data:`T`)
        :type z0: :class:`float`
        :rtype: :class:`NPortMatrix`

        """
        # references:
        #  * http://qucs.sourceforge.net/tech/node98.html
        #           "Transformations of n-Port matrices"
        z0 = self.convert_z0test(type, z0)
        idty = np.identity(len(self), dtype=complex)
        invert = np.linalg.inv

        if type not in (Z, Y, S, T, H, G, ABCD):
            raise TypeError("Unknown 2-port parameter type")

        # TODO: check for singularities (condition number)
        result = None
        if self.type == SCATTERING:
            s11 = self[0, 0]
            s12 = self[0, 1]
            s21 = self[1, 0]
            s22 = self[1, 1]
            if type == HYBRID:
                h11 = ((1 + s11) * (1 + s22) - s12 * s21) * self.z0
                h12 = 2 * s12
                h21 = -2 * s21
                h22 = ((1 - s11) * (1 - s22) - s12 * s21) / self.z0
                d = (1 - s11) * (1 + s22) + s12 * s21
                result = np.asarray([[h11, h12], [h21, h22]]) / d
            elif type == INVERSE_HYBRID:
                g11 = ((1 - s11) * (1 - s22) - s12 * s21) / self.z0
                g12 = -2 * s12
                g21 = 2 * s21
                g22 = ((1 + s11) * (1 + s22) - s12 * s21) * self.z0
                d = (1 + s11) * (1 - s22) + s12 * s21
                result = np.asarray([[g11, g12], [g21, g22]]) / d
        elif self.type == SCATTERING_TRANSFER:
            t11 = self[0, 0]
            t12 = self[0, 1]
            t21 = self[1, 0]
            t22 = self[1, 1]
            if type == HYBRID:
                h11 = (- t11 + t12 - t21 + t22) * self.z0
                h12 = - 2 * (t11 * t22 - t12 * t21)
                h21 = - 2
                h22 = (t11 + t12 - t21 - t22) / self.z0
                d = - t11 + t12 + t21 - t22
                result = np.asarray([[h11, h12], [h21, h22]]) / d
            elif type == INVERSE_HYBRID:
                g11 = (t11 + t12 - t21 - t22) / self.z0
                g12 = 2 * (t11 * t22 - t12 * t21)
                g21 = 2
                g22 = (- t11 + t12 - t21 + t22) * self.z0
                d = t11 + t12 + t21 + t22
                result = np.asarray([[g11, g12], [g21, g22]]) / d
        elif self.type == IMPEDANCE:
            z11 = self[0, 0]
            z12 = self[0, 1]
            z21 = self[1, 0]
            z22 = self[1, 1]
            dz = z11 * z22 - z12 * z21
            if type == HYBRID:
                result = np.asarray([[dz, z12], [-z21, 1]]) / z22
            elif type == INVERSE_HYBRID:
                result = np.asarray([[1, -z12], [z21, dz]]) / z11
        elif self.type == ADMITTANCE:
            y11 = self[0, 0]
            y12 = self[0, 1]
            y21 = self[1, 0]
            y22 = self[1, 1]
            dy = y11 * y22 - y12 * y21
            if type == HYBRID:
                result = np.asarray([[1, -y12], [y21, dy]]) / y11
            elif type == INVERSE_HYBRID:
                result = np.asarray([[dy, y12], [-y21, 1]]) / y22
        if self.type == TRANSMISSION:
            a11 = self[0, 0]
            a12 = self[0, 1]
            a21 = self[1, 0]
            a22 = self[1, 1]
            da = a11 * a22 - a12 * a21
            if type == HYBRID:
                result = np.asarray([[a12, da], [-1, a21]]) / a22
            elif type == INVERSE_HYBRID:
                result = np.asarray([[a21, -da], [1, a12]]) / a11
        elif self.type == HYBRID:
            h11 = self[0, 0]
            h12 = self[0, 1]
            h21 = self[1, 0]
            h22 = self[1, 1]
            dh = h11 * h22 - h12 * h21
            if type == SCATTERING:
                h11_ = h11 / z0
                h22_ = h22 * z0
                s11 = (h11_ - 1) * (1 + h22_) - h12 * h21
                s12 = 2 * h12
                s21 = -2 * h21
                s22 = (h11_ + 1) * (1 - h22_) + h12 * h21
                d = (h11_ + 1) * (1 + h22_) - h12 * h21
                result = np.asarray([[s11, s12], [s21, s22]]) / d
            elif type == SCATTERING_TRANSFER:
                h11_ = h11 / z0
                h22_ = h22 * z0
                t11 = + (h11_ + 1) * (1 - h22_) + h12 * h21
                t12 = - (h11_ + 1) * (1 + h22_) + h12 * h21
                t21 = + (h11_ - 1) * (1 - h22_) + h12 * h21
                t22 = - (h11_ - 1) * (1 + h22_) + h12 * h21
                result = np.asarray([[t11, t12], [t21, t22]]) / (2 * h21)
            elif type == IMPEDANCE:
                result = np.asarray([[dh, h12], [-h21, 1]]) / h22
            elif type == ADMITTANCE:
                result = np.asarray([[1, -h12], [h21, dh]]) / h11
            elif type == TRANSMISSION:
                result = np.asarray([[-dh, -h11], [-h22, -1]]) / h21
            elif type == HYBRID:
                result = self
            elif type == INVERSE_HYBRID:
                result = np.asarray([[h22, -h12], [-h21, h11]]) / dh
        elif self.type == INVERSE_HYBRID:
            g11 = self[0, 0]
            g12 = self[0, 1]
            g21 = self[1, 0]
            g22 = self[1, 1]
            dg = g11 * g22 - g12 * g21
            if type == SCATTERING:
                g11_ = g11 * z0
                g22_ = g22 / z0
                s11 = (1 - g11_) * (g22_ + 1) + g12 * g21
                s12 = -2 * g12
                s21 = 2 * g21
                s22 = (1 + g11_) * (g22_ - 1) - g12 * g21
                d = (1 + g11_) * (g22_ + 1) - g12 * g21
                result = np.asarray([[s11, s12], [s21, s22]]) / d
            elif type == SCATTERING_TRANSFER:
                g11_ = g11 * z0
                g22_ = g22 / z0
                t11 = + (g11_ + 1) * (1 - g22_) + g12 * g21
                t12 = + (g11_ + 1) * (1 + g22_) - g12 * g21
                t21 = - (g11_ - 1) * (1 - g22_) - g12 * g21
                t22 = - (g11_ - 1) * (1 + g22_) + g12 * g21
                result = np.asarray([[t11, t12], [t21, t22]]) / (2 * g21)
            elif type == IMPEDANCE:
                result = np.asarray([[1, -g12], [g21, dg]]) / g11
            elif type == ADMITTANCE:
                result = np.asarray([[dg, g12], [-g21, 1]]) / g22
            elif type == TRANSMISSION:
                result = np.asarray([[1, g22], [g11, dg]]) / g21
            elif type == HYBRID:
                result = np.asarray([[g22, -g12], [-g21, g11]]) / dg
            elif type == INVERSE_HYBRID:
                result = self
                
        if result is not None:
            return NPortMatrix(result, type, z0)
        else:
            nportmatrix = super(self.__class__, self)
            return nportmatrix.twonportmatrix().convert(type, z0).nportmatrix()


class NPort(NPortBase):
    """Class representing an n-port across a list of frequencies
    
    :param freqs: list of frequency samples
    :type freqs: :class:`list`
    :param matrix: list of matrix elements
    :type matrix: list of *n* by *n* arrays
    :param type: matrix type
    :type type: :data:`Z`, :data:`Y`, :data:`S`, :data:`T`, :data:`G`, :data:`H`
                or :data:`ABCD`
    :param z0: normalizing impedance (only :data:`S`)
    :type z0: :class:`float`

    """
    def __new__(cls, freqs, matrices, type, z0=None):
        matrices = np.asarray(matrices, dtype=complex)
        obj = NPortBase.__new__(cls, freqs, matrices, type, z0)
        if len(obj[0].shape) != 2:
            raise ValueError("the matrices should be two-dimensional")
        return obj

    @property
    def ports(self):
        """The number of ports of this :class:`NPort`
        
        :rtype: :class:`int`
        
        """
        return self[0].shape[0]

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
        """Convert this NPort to a TwoNPort using inports as the input ports
        and outports as the output ports.
        
        :param inports: the list of ports that make up the inputs of the 2n-port
        :type inports: :class:`tuple` or :class:`list`
        :param outports: the list of ports that make up the outputs
        :type outports: :class:`tuple` or :class:`list`
        :rtype: :class:`TwoNPort`

        """
        twonportmatrices = []
        for matrix in self:
            nportmatrix = NPortMatrix(matrix, self.type, self.z0)
            twonportmatrices.append(nportmatrix.twonportmatrix(inports,
                                                               outports))
        return twonport.TwoNPort(self.freqs, twonportmatrices, self.type,
                                 self.z0)

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
            renormalized = []
            for matrix in self:
                nportmatrix = NPortMatrix(matrix, self.type, self.z0)
                renormalized.append(nportmatrix.renormalize(z0))
            result = self.__class__(self.freqs, renormalized, self.type, z0)
        return result

    def convert(self, type, z0=None):
        """Convert to another n-port matrix representation
        
        :param type: n-port representation type to convert to
        :type type: :data:`Z`, :data:`Y`, :data:`S`, :data:`T`, :data:`G`,
                    :data:`H` or :data:`ABCD`
        :param z0: normalizing impedance (only :data:`S` and :data:`T`)
        :type z0: :class:`float`

        """
        converted = []
        for matrix in self:
            nportmatrix = NPortMatrix(matrix, self.type, self.z0)
            converted.append(nportmatrix.convert(type, z0))
        return NPort(self.freqs, converted, type, z0)

    def submatrix(self, ports):
        """Keep only the parameters corresponding to the given ports, discarding
        the others
        
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
        inverted = []
        for matrix in self:
            nportmatrix = NPortMatrix(matrix, self.type, self.z0)
            inverted.append(np.linalg.inv(nportmatrix))
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
        recombined = []
        for i, matrix in enumerate(self):
            nportmatrix = NPortMatrix(matrix, self.type, self.z0)
            recomb = nportmatrix.recombine(portsets)
            recombined.append(recomb)
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
        shunted = []
        for i, matrix in enumerate(self):
            nportmatrix = NPortMatrix(matrix, self.type, self.z0)
            shunt = nportmatrix.shunt(portsets)
            shunted.append(shunt)
        return self.__class__(self.freqs, shunted, self.type, self.z0)

    def ispassive(self):
        """Check whether this :class:`NPort` is passive
        
        :rtype: :class:`bool`
        
        """
        for nportmatrix in self:
            if not NPortMatrix(nportmatrix, self.type, self.z0).ispassive():
                return False
        return True


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
    
    if type(arg1) == NPort:
        if type(arg2) == NPort:
            result_freqs = merge_freqs(arg1.freqs, arg2.freqs)
            arg1_matrices = arg1.at(result_freqs)
            arg2_matrices = arg2.at(result_freqs)
            result_matrices = np.asarray([np.dot(a, b) for (a,b) in
                zip(arg1_matrices, arg2_matrices)])
        else:
            result_freqs = arg1.freqs
            result_matrices = np.array([np.dot(matrix, other)
                for matrix in arg1.matrices])
        return NPort(result_freqs, result_matrices, arg1.type, arg1.z0)
    elif type(arg1) == twonport.TwoNPort:
        if type(arg2) == twonport.TwoNPort:
            def twonport_dot(l, r):
                a = np.dot(l[0,0], r[0,0]) + np.dot(l[0,1], r[1,0])
                b = np.dot(l[0,0], r[0,1]) + np.dot(l[0,1], r[1,1])
                c = np.dot(l[1,0], r[0,0]) + np.dot(l[1,1], r[1,0])
                d = np.dot(l[1,0], r[0,1]) + np.dot(l[1,1], r[1,1])
                return np.asarray([[a, b], [c, d]])
            
            result_freqs = merge_freqs(arg1.freqs, arg2.freqs)
            arg1_matrices = arg1.at(result_freqs)
            arg2_matrices = arg2.at(result_freqs)
            result_matrices = np.asarray([twonport_dot(a, b) for (a,b) in
                zip(arg1_matrices, arg2_matrices)])
        else:
            raise NotImplementedError
        return twonport.TwoNPort(result_freqs, result_matrices, arg1.type,
                                 arg1.z0)
    elif type(arg2) == NPort:
        raise NotImplementedError
    else:
        np.dot(arg1, arg2)
