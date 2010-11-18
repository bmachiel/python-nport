from __future__ import division

import numpy as np
import nport
from twonport import TwoNPort


class TransmissionLine(object):
    """Class representing a two-conductor transmission line. It takes the
    2-port parameters of a transmission line and its length, and it calculates
    the transmission lines':

    * propagation constant :attr:`gamma`, 
    * characteristic impedance :attr:`z0`, and
    * per-unit-length parameters (:attr:`rpm`, :attr:`lpm`, :attr:`gpm`, :attr:`cpm`)

    """
    def __init__(self, twonport, length):
        """
        :param twonport: 2-port that represents a transmission line
        :type twonport: TwoNPort
        :param length: physical length of the transmission line in meters
        :type length: float
        
        """
        assert twonport.ports == 2
        self.twonport = twonport.convert(nport.TRANSMISSION)
        self.freqs = self.twonport.freqs
        self.length = length
        
        # retrieve ABCD parameters
        self.a = self.twonport.get_parameter(1, 1)[:,0,0]
        self.b = self.twonport.get_parameter(1, 2)[:,0,0]
        self.c = self.twonport.get_parameter(2, 1)[:,0,0]
        self.d = self.twonport.get_parameter(2, 2)[:,0,0]

        # the arccosh can be implemented in a number of different ways
        numpy_acosh = np.arccosh

        def mathworld_acosh(z):
            # http://mathworld.wolfram.com/InverseHyperbolicCosine.html
            return np.log(z + np.sqrt(z + 1) * np.sqrt(z - 1))

        def matlab_acosh(z):
            #MATLAB
            return np.log(z + np.sqrt(z*z - 1))
            
        acosh = matlab_acosh

        tmp = (self.a + self.d) / 2.0
        self._gamma = acosh(tmp) / self.length
        self._z0 = np.sqrt(self.b / self.c)

        # TODO: unwrap gamma (*before* division by length)
        
        # extract RLGC parameters [EIS92]
        self._rpm = (self.gamma * self.z0).real
        self._lpm = (self.gamma * self.z0).imag / (2 * np.pi * self.freqs)
        self._gpm = (self.gamma / self.z0).real
        self._cpm = (self.gamma / self.z0).imag / (2 * np.pi * self.freqs)

    @property
    def gamma(self):
        """Propagation constant"""
        return self._gamma

    @property
    def z0(self):
        """Characteristic impedance"""
        return self._z0

    @property
    def rpm(self):
        """Resistance per meter"""
        return self._rpm

    @property
    def lpm(self):
        """Inductance per meter"""
        return self._lpm

    @property
    def gpm(self):
        """Conductance per meter"""
        return self._gpm

    @property
    def cpm(self):
        """Capacitance per meter"""
        return self._cpm
