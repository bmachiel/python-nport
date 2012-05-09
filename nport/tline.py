from __future__ import division

import numpy as np
from .base import TRANSMISSION
from twonport import TwoNPort


class TransmissionLine(object):
    """Class representing a two-conductor transmission line. It takes the
    2-port parameters of a transmission line and its length, and it calculates
    the transmission lines':

    * propagation constant :attr:`gamma`, 
    * characteristic impedance :attr:`z0`, and
    * per-unit-length parameters (:attr:`rpm`, :attr:`lpm`, :attr:`gpm`,
                                  :attr:`cpm`)

    """
    def __init__(self, twonport, length, reciprocal=True):
        """
        :param twonport: 2-port that represents a transmission line
        :type twonport: TwoNPort
        :param length: physical length of the transmission line in meters
        :type length: float
        
        """
        if not isinstance(twonport, TwoNPort) or twonport.ports != 2:
            raise TypeError
        self.twonport = twonport.convert(TRANSMISSION)
        self.freqs = self.twonport.freqs
        self.length = length
        
        # retrieve ABCD parameters
        self.a = self.twonport.get_parameter(1, 1)[:,0,0]
        self.b = self.twonport.get_parameter(1, 2)[:,0,0]
        self.c = self.twonport.get_parameter(2, 1)[:,0,0]
        self.d = self.twonport.get_parameter(2, 2)[:,0,0]

        sum = (self.a + self.d) / 2.0
        diff = (self.d - self.a) / 2.0
        ad_bc = self.a * self.d - self.b * self.c
       
        delta = unwrap_sqrt(sum**2 - ad_bc)
        exp_gl_forward = (sum + delta) / ad_bc
        exp_gl_backward = (sum - delta) / ad_bc

        self._gamma_forward = unwrap_log(exp_gl_forward) / self.length
        self._gamma_backward = - unwrap_log(exp_gl_backward) / self.length

        self._z0_forward = - self.b / (diff - delta)
        self._z0_backward = self.b / (diff + delta)

        # extract RLGC parameters [EIS92]
        two_pi_f = 2 * np.pi * self.freqs
        gamma_mul_z0_forward = self._gamma_forward * self._z0_forward
        gamma_div_z0_forward = self._gamma_forward / self._z0_forward
        gamma_mul_z0_backward = self._gamma_backward * self._z0_backward
        gamma_div_z0_backward = self._gamma_backward / self._z0_backward

        self._rpm_forward = gamma_mul_z0_forward.real
        self._lpm_forward = gamma_mul_z0_forward.imag / two_pi_f
        self._gpm_forward = gamma_div_z0_forward.real
        self._cpm_forward = gamma_div_z0_forward.imag / two_pi_f

        self._rpm_backward = gamma_mul_z0_backward.real
        self._lpm_backward = gamma_mul_z0_backward.imag / two_pi_f
        self._gpm_backward = gamma_div_z0_backward.real
        self._cpm_backward = gamma_div_z0_backward.imag / two_pi_f

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


def unwrap_sqrt(arg):
    """square root of complex numbers (first unwrap)"""
    mag, ang = np.abs(arg), np.unwrap(np.angle(arg))
    return np.sqrt(mag) * np.exp(1j * ang / 2)


def unwrap_log(arg):
    """natural logarithm of complex numbers (first unwrap)"""
    mag, ang = np.abs(arg), np.unwrap(np.angle(arg))
    return np.log(mag) + 1j * ang
