
from __future__ import division

import numpy as np
from scipy import interpolate


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
    :param z0: normalizing impedance (only :class:`S` and :class:`T`)
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

    def __reduce__(self):
        object_state = list(np.ndarray.__reduce__(self))
        subclass_state = (self.type, self.z0, )
        object_state[2] = (object_state[2], subclass_state)
        return tuple(object_state)
    
    def __setstate__(self,state):
        nd_state, own_state = state
        np.ndarray.__setstate__(self, nd_state)
        self.type, self.z0 = own_state

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
    matrix_cls = NPortMatrixBase
    
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
        obj = super(NPortBase, cls).__new__(cls, matrices, type, z0)
        obj.freqs = np.asarray(freqs)
        return obj

    def __array_finalize__(self, obj):
        super(NPortBase, self).__array_finalize__(obj)
        self.freqs = getattr(obj, 'freqs', None)

    def __getitem__(self, arg):
        if type(arg) == int:
            return self.matrix_cls(np.asarray(self).__getitem__(arg), self.type,
                                   self.z0)
        else:
            return np.asarray(self).__getitem__(arg)

    def __reduce__(self):
        object_state = list(super(NPortBase, self).__reduce__())
        subclass_state = (self.freqs, )
        object_state[2] = (object_state[2], subclass_state)
        return tuple(object_state)
    
    def __setstate__(self,state):
        super_state, own_state = state
        super(NPortBase, self).__setstate__(super_state)
        self.freqs, = own_state

    #~ def __repr__(self):
        #~ return "freqs: " + repr(self.freqs) + "\n" + NPortMatrixBase.__repr__(self)

    #~ def __str__(self):
        #~ return self.__repr__()

    def get_parameter(self, port1, port2):
        """Return the parameter as specified by the indices `port1` and `port2`
        as an ndarray
        
        :param port1: first index
        :type port1: :class:`int`
        :param port2: second index
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
