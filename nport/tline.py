from __future__ import division

import numpy as np
from .base import TRANSMISSION
from twonport import TwoNPort



class cached_property(property):
    def __init__(self, function, *args, **kwargs):
        super(cached_property, self).__init__(function, *args, **kwargs)
        self._function_name = function.__name__
        
    def __get__(self, obj, *args):
        cache_variable = '_' + self._function_name
        try:
            return getattr(obj, cache_variable)
        except AttributeError:
            cache_value = super(cached_property, self).__get__(obj, *args)
            setattr(obj, cache_variable, cache_value)
            return cache_value


class property_if_reciprocal(property):
    def __init__(self, function, *args, **kwargs):
        super(property_if_reciprocal, self).__init__(function, *args, **kwargs)
        
    def __get__(self, obj, *args):
        if obj.reciprocal:
            return super(property_if_reciprocal, self).__get__(obj, *args)
        else:
            raise AttributeError('Non-reciprocal transmission line. Use the '
                                 'forward and backward properties instead')


class TransmissionLine(object):
    """
    """
    def __init__(self, freqs, reciprocal=False):
        self.freqs = np.asarray(freqs)
        self.reciprocal = reciprocal

    @cached_property
    def two_pi_f(self):
        return 2 * np.pi * self.freqs

    @property_if_reciprocal
    def gamma(self):
        return self.gamma_forward

    @property_if_reciprocal
    def z0(self):
        return self.z0_forward

    @property_if_reciprocal
    def z(self):
        return self.z_forward

    @property_if_reciprocal
    def y(self):
        return self.y_forward

    @property_if_reciprocal
    def r(self):
        return self.r_forward

    @property_if_reciprocal
    def l(self):
        return self.l_forward

    @property_if_reciprocal
    def g(self):
        return self.g_forward

    @property_if_reciprocal
    def c(self):
        return self.c_forward
        
    def twoport(self, length):
        from .twoport import TwoPort
        if not self.reciprocal:
            egl_fw = np.exp(self.gamma_forward * length)
            emgl_bw = np.exp(- self.gamma_backward * length)
            z0_fw, z0_bw = self.z0_forward, self.z0_backward
            denom = z0_fw + z0_bw
            a = (egl_fw * z0_fw + emgl_bw * z0_bw) / denom
            c = (egl_fw - emgl_bw) / denom
            b = c * z0_fw * z0_bw
            d = (egl_fw * z0_bw + emgl_bw * z0_fw) / denom
            abcd_matrices = np.hstack((np.array([a]).T, np.array([b]).T,
                                       np.array([c]).T, np.array([d]).T))
        else:
            gamma_length = self.gamma * length
            # FIXME: A = D only for symmetric transmission lines!
            a = d = np.cosh(gamma_length)
            b = self.z0 * np.sinh(gamma_length)
            c = 1 / self.z0 * np.sinh(gamma_length)
            abcd_matrices = np.hstack((np.array([a]).T, np.array([b]).T,
                                       np.array([c]).T, np.array([d]).T))
        return TwoPort(self.freqs, abcd_matrices.reshape(-1,2,2), TRANSMISSION)


class GammaZ0TransmissionLine(TransmissionLine):
    def __init__(self, freqs, gamma_forward, z0_forward,
                 gamma_backward=None, z0_backward=None):
        reciprocal = all(map(lambda x: x is None,
                             (gamma_backward, z0_backward)))
        super(GammaZ0TransmissionLine, self).__init__(freqs, reciprocal)
        self.gamma_forward = np.asarray(gamma_forward)
        self.z0_forward = np.asarray(z0_forward)
        if not reciprocal:
            if None in (gamma_backward, z0_backward):
                raise ValueError('Both gamma_backward and z0_backward need to '
                                 'be specified')
            self.gamma_backward = np.asarray(gamma_backward)
            self.z0_backward = np.asarray(z0_backward)
        else:
            self.gamma_backward = gamma_forward
            self.z0_backward = z0_forward

    @cached_property
    def z_forward(self):
        return self.gamma_forward * self.z0_forward
        
    @cached_property
    def y_forward(self):
        return self.gamma_forward / self.z0_forward
        
    @cached_property
    def z_backward(self):
        return self.gamma_backward * self.z0_backward
        
    @cached_property
    def y_backward(self):
        return self.gamma_backward / self.z0_backward

    @cached_property
    def r_forward(self):
        return self.z_forward.real

    @cached_property
    def l_forward(self):
        return self.z_forward.imag / self.two_pi_f

    @cached_property
    def g_forward(self):
        return self.y_forward.real

    @cached_property
    def c_forward(self):
        return self.y_forward.imag / self.two_pi_f

    @cached_property
    def r_backward(self):
        return self.z_backward.real

    @cached_property
    def l_backward(self):
        return self.z_backward.imag / self.two_pi_f

    @cached_property
    def g_backward(self):
        return self.y_backward.real

    @cached_property
    def c_backward(self):
        return self.y_backward.imag / self.two_pi_f


class RLGCTransmissionLine(TransmissionLine):
    def __init__(self, freqs, r_forward, l_forward, g_forward, c_forward,
                 r_backward=None, l_backward=None, g_backward=None,
                 c_backward=None):
        rlgc_backward = (r_backward, l_backward, g_backward, c_backward)
        reciprocal = all(map(lambda x: x is None, rlgc_backward))
        super(RLGCTransmissionLine, self).__init__(freqs, reciprocal)
        self.r_forward = np.asarray(r_forward)
        self.l_forward = np.asarray(l_forward)
        self.g_forward = np.asarray(g_forward)
        self.c_forward = np.asarray(c_forward)
        if not reciprocal:
            if None in rlgc_backward:
                raise ValueError('One or more backward parameters are missing')
            self.r_backward = np.asarray(r_backward)
            self.l_backward = np.asarray(l_backward)
            self.g_backward = np.asarray(g_backward)
            self.c_backward = np.asarray(c_backward)
        else:
            self.r_backward = self.r_forward
            self.l_backward = self.l_forward
            self.g_backward = self.g_forward
            self.c_backward = self.c_forward

    @cached_property
    def z_forward(self):
        return self.r_forward + 1j * self.two_pi_f * self.l_forward

    @cached_property
    def y_forward(self):
        return self.g_forward + 1j * self.two_pi_f * self.c_forward

    @cached_property
    def gamma_forward(self):
        return unwrap_sqrt(self.z_forward * self.y_forward)
        
    @cached_property
    def z0_forward(self):
        return unwrap_sqrt(self.z_forward / self.y_forward)

    @cached_property
    def z_backward(self):
        return self.r_backward + 1j * self.two_pi_f * self.l_backward

    @cached_property
    def y_backward(self):
        return self.g_backward + 1j * self.two_pi_f * self.c_backward

    @cached_property
    def gamma_backward(self):
        return unwrap_sqrt(self.z_backward * self.y_backward)
        
    @cached_property
    def z0_backward(self):
        return unwrap_sqrt(self.z_backward / self.y_backward)


def unwrap_sqrt(arg):
    """square root of complex numbers (first unwrap)"""
    mag, ang = np.abs(arg), np.unwrap(np.angle(arg))
    return np.sqrt(mag) * np.exp(1j * ang / 2)


def unwrap_log(arg):
    """natural logarithm of complex numbers (first unwrap)"""
    mag, ang = np.abs(arg), np.unwrap(np.angle(arg))
    return np.log(mag) + 1j * ang
