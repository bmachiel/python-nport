from __future__ import division

import numpy as np
from .base import TRANSMISSION
from twonport import TwoNPort
from eigenshuffle import eigenshuffle



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



class MulticonductorTransmissionLine(object):
    """Class representing a multiconductor transmission line. It calculates:
    
    * per-unit-length matrices (R, L, G, C), 
    * modal propagation constant matrices and characteristic impedance matrices,
      and
    * natural propagation constant matrices and characteristic impedance
      matrices

    For non-uniform lines, the characteristic impedance matrices differ for
    forward and backward propagation:
    
    * :attr:`modal_z0_forward` and :attr:`natural_z0_forward`
    * :attr:`modal_z0_backward` and :attr:`natural_z0_backward`

    For reciprocal lines, the propagation constant matrices are the same for
    forward and backward propagation, while for non-reciprocal lines they
    differ.

    """
    def __init__(self, twonport, length, reciprocal=True):
        """
        :param twonport: 2n-port that represents a multiconductor transmission
                         line
        :type twonport: TwoNPort
        :param length: physical length of the transmission line in meters
        :type length: float
        :param reciprocal: True if `twonport` represents a reciprocal line
        :type reciprocal: bool
        
        """
        self.twonport = twonport.convert(TRANSMISSION)
        self.freqs = self.twonport.freqs
        self.length = length
        
        # modal analysis of an MTL
        # reference:
        # [FAR04] "A new generalized modal analysis theory for nonuniform
        #   multiconductor transmission lines" by J.A. Brandao Faria
        # other references:
        # [FAR93] "Multiconductor Transmission-Line Structures - Modal Analysis
        #   Techniques" by J.A. Brandao Faria
        # [PAU08] "Analysis of Multiconductor Transmission Lines", 2nd edition
        #   by Clayton R. Paul
        
        # retrieve ABCD parameters
        self.a = self.twonport.get_parameter(1, 1)
        self.b = self.twonport.get_parameter(1, 2)
        self.c = self.twonport.get_parameter(2, 1)
        self.d = self.twonport.get_parameter(2, 2)
        
        # calculate eigenvalues and eigenvectors
        if reciprocal:
            self.b_dot_ct = np.array([np.dot(b, c.T)
                                      for b, c in zip(self.b, self.c)])
            self.bt_dot_c = np.array([np.dot(b.T, c)
                                      for b, c in zip(self.b, self.c)])

            self.e0, self.t0 = eigenshuffle(self.b_dot_ct)
            self.el, self.tl = eigenshuffle(self.bt_dot_c)
        else:
            # [FAR04] assumes reciprocal lines, so we need to use more complex equations
            inv = True
            if not inv:
                # these result in the inverse eigenvalues of the method for
                # reciprocal lines and ai_dot_b and d_dot_bi below
                # TODO: understand why and explain
                a_dot_ci = np.array([np.dot(a, np.linalg.inv(c))
                                     for a, c in zip(self.a, self.c)])
                d_dot_bi = np.array([np.dot(d, np.linalg.inv(b))
                                     for d, b in zip(self.d, self.b)])
                t0_prod = np.array([np.dot(aci, dbi)
                                    for aci, dbi in zip(a_dot_ci, d_dot_bi)])
                w0_prod = np.array([np.dot(dbi, aci)
                                    for dbi, aci in zip(d_dot_bi, a_dot_ci)])
                # hence we need to invert the eigenvalues ...
                t0_prod = np.array([np.linalg.inv(t0p) for t0p in t0_prod])
                w0_prod = np.array([np.linalg.inv(w0p) for w0p in w0_prod])
            else:
                # ... or calculate the eigenvalues from the inverse matrix
                c_dot_ai = np.array([np.dot(c, np.linalg.inv(a))
                                     for a, c in zip(self.a, self.c)])
                b_dot_di = np.array([np.dot(b, np.linalg.inv(d))
                                     for b, d in zip(self.b, self.d)])
                t0_prod = np.array([np.dot(bdi, cai)
                                    for cai, bdi in zip(c_dot_ai, b_dot_di)])
                w0_prod = np.array([np.dot(cai, bdi)
                                    for cai, bdi in zip(c_dot_ai, b_dot_di)])

            ai_dot_b = np.array([np.dot(np.linalg.inv(a), b)
                                 for a, b in zip(self.a, self.b)])
            di_dot_c = np.array([np.dot(np.linalg.inv(d), c)
                                 for d, c in zip(self.d, self.c)])
            tl_prod = np.array([np.dot(aib, dic)
                                for aib, dic in zip(ai_dot_b, di_dot_c)])
            wl_prod = np.array([np.dot(dic, aib)
                                for dic, aib in zip(di_dot_c, ai_dot_b)])

            # shift eigenvalues for numeric stability
            t0_prod_shifted, t0_traces = shift_eigenvalues(t0_prod)
            tl_prod_shifted, tl_traces = shift_eigenvalues(tl_prod)
            w0_prod_shifted, w0_traces = shift_eigenvalues(w0_prod)
            wl_prod_shifted, wl_traces = shift_eigenvalues(wl_prod)
    
            e0_shifted, self.t0 = eigenshuffle(t0_prod_shifted)
            el_shifted, self.tl = eigenshuffle(tl_prod_shifted)
            ew0_shifted, self.w0 = eigenshuffle(w0_prod_shifted)
            ewl_shifted, self.wl = eigenshuffle(wl_prod_shifted)

            self.e0 = np.array([eig + trace
                                for eig, trace in zip(e0_shifted, t0_traces)])
            self.el = np.array([eig + trace
                                for eig, trace in zip(el_shifted, tl_traces)])
            self.ew0 = np.array([eig + trace
                                 for eig, trace in zip(ew0_shifted, w0_traces)])
            self.ewl = np.array([eig + trace
                                 for eig, trace in zip(ewl_shifted, wl_traces)])

        self.t0inv = np.array([np.linalg.inv(t0) for t0 in self.t0])
        self.tlinv = np.array([np.linalg.inv(tl) for tl in self.tl])

        if reciprocal:
            self.w0 = np.array([np.transpose(t0inv) for t0inv in self.t0inv])
            self.wl = np.array([np.transpose(tlinv) for tlinv in self.tlinv])

        self.w0inv = np.array([np.linalg.inv(w0) for w0 in self.w0])
        #self.wlinv = np.array([np.linalg.inv(wl) for wl in self.wl])

        self.am = np.array([np.dot(np.dot(t0inv, a), tl)
                            for t0inv, a, tl
                            in zip(self.t0inv, self.a, self.tl)])
        self.bm = np.array([np.dot(np.dot(t0inv, b), wl)
                            for t0inv, b, wl
                            in zip(self.t0inv, self.b, self.wl)])
        self.cm = np.array([np.dot(np.dot(w0inv, c), tl)
                            for w0inv, c, tl
                            in zip(self.w0inv, self.c, self.tl)])
        self.dm = np.array([np.dot(np.dot(w0inv, d), wl)
                            for w0inv, d, wl
                            in zip(self.w0inv, self.d, self.wl)])
        
        # calculate modal propagation factors and characteristic impedances
        diag_am = np.asarray([np.diag(am) for am in self.am])
        diag_bm = np.asarray([np.diag(bm) for bm in self.bm])
        diag_cm = np.asarray([np.diag(cm) for cm in self.cm])
        diag_dm = np.asarray([np.diag(dm) for dm in self.dm])

        sum = (diag_am + diag_dm) / 2.0
        if reciprocal:
            # TODO detect uniform MTLs
            #   based on T0, Tl, W0, Wl ?
            #   np.max(self.modal_z0_forward - self.modal_z0_backward)
            delta = unwrap_sqrt(sum**2 - 1)
        else:
            ad_bc = diag_am * diag_dm - diag_bm * diag_cm
            delta = unwrap_sqrt(sum**2 - ad_bc)
           
        exp_gl_forward = sum + delta
        exp_gl_backward = sum - delta
        
        self.modal_gamma_forward = unwrap_log(exp_gl_forward) / self.length
        self.modal_gamma_backward = - unwrap_log(exp_gl_backward) / self.length

        self.modal_z0_forward = (exp_gl_forward - diag_dm) / diag_cm
        self.modal_z0_backward = (diag_dm - exp_gl_backward) / diag_cm
        #~ if reciprocal:
            #~ self.modal_z0_forward = diag_bm / (exp_gl_forward - diag_am)
            #~ self.modal_z0_backward = diag_bm / (diag_am - exp_gl_backward)

        self.modal_y0_forward = 1.0 / self.modal_z0_forward
        self.modal_y0_backward = 1.0 / self.modal_z0_backward
        
        # TODO: unwrap gamma

        # calculate (natural) longitudinal RLGC matrices
        # NOTE: verify!
        #   z0_forward and z0_backward?
        #   t0, tl, w0, wl?
        # 1) compute gamma and Z0 matrices in natural coordinates
        #     [FAR93] eq 1.59, [PAU08] eq 7.77 (?)
        self.natural_gamma_forward = \
            np.asarray([np.dot(np.dot(t0, np.diag(gam)), tlinv)
                        for t0, gam, tlinv
                        in zip(self.t0, self.modal_gamma_forward, self.tlinv)])
#        self.natural_gamma_backward = \  # this is a guess
#            np.asarray([np.dot(np.dot(tl, np.diag(gam)), t0inv)
#                        for tl, gam, t0inv
#                        in zip(self.tl, self.modal_gamma, self.t0inv)])

        self.modal_y0m_forward = np.asarray([np.diag(y0)
                                             for y0
                                             in self.modal_y0_forward])
        self.modal_y0m_backward = np.asarray([np.diag(y0)
                                              for y0
                                              in self.modal_y0_backward])
            # the modal wave admittance matrix (diagonal)

        self.natural_y0_forward = np.asarray([np.dot(np.dot(w0, y0), tlinv)
                                              for w0, y0, tlinv
                                              in zip(self.w0,
                                                     self.modal_y0m_forward,
                                                     self.tlinv)])
#        self.natural_y0_backward = np.asarray([np.dot(np.dot(wl, y0), t0inv) # this is a guess
#            for wl, y0, t0inv in zip(self.wl, self.modal_y0m_backward, self.t0inv)])
            # conversion to natural wave admittance matrices [FAR93] eq 1.57
        self.natural_z0_forward = np.asarray([np.linalg.inv(y0) 
                                              for y0
                                              in self.natural_y0_forward])
            # natural wave impedance matrices

        # 2) calculate per-unit-length Z, Y and RLGC matrices
        #     [FAR93] eq 1.58 & 1.60
        self.zpm_forward = np.asarray([np.dot(gam, zw)
                                       for gam, zw
                                       in zip(self.natural_gamma_forward,
                                              self.natural_z0_forward)])
        self.ypm_forward = np.asarray([np.dot(yw, gam)
                                       for yw, gam
                                       in zip(self.natural_y0_forward,
                                              self.natural_gamma_forward)])
        
        two_pi_freqs = 2 * np.pi * self.freqs
        self._rpm_forward = self.zpm_forward.real
        self._lpm_forward = (self.zpm_forward.imag.T / two_pi_freqs).T
        self._gpm_forward = self.ypm_forward.real
        self._cpm_forward = (self.ypm_forward.imag.T / two_pi_freqs).T

        self.modal_gamma = self.modal_gamma_forward
        
    @property
    def gamma_forward(self):
        """Forward propagation constant"""
        return self._gamma_forward

    @property
    def z0_forward(self):
        """Forward characteristic impedance"""
        return self._z0_forward

    @property
    def rpm_forward(self):
        """Resistance per meter for forward traveling waves"""
        return self._rpm_forward

    @property
    def lpm_forward(self):
        """Inductance per meter for forward traveling waves"""
        return self._lpm_forward

    @property
    def gpm_forward(self):
        """Conductance per meter for forward traveling waves"""
        return self._gpm_forward

    @property
    def cpm_forward(self):
        """Capacitance per meter for forward traveling waves"""
        return self._cpm_forward


def unwrap_sqrt(arg):
    """square root of complex numbers (first unwrap)"""
    mag, ang = np.abs(arg), np.unwrap(np.angle(arg))
    return np.sqrt(mag) * np.exp(1j * ang / 2)


def unwrap_log(arg):
    """natural logarithm of complex numbers (first unwrap)"""
    mag, ang = np.abs(arg), np.unwrap(np.angle(arg))
    return np.log(mag) + 1j * ang


def shift_eigenvalues(matrices):
    """Shift eigenvalues of a matrix to make eigenvalue calculation numerically
    more stable
    
    """
    # [FAR04] 1.3.2
    n = matrices.shape[-1]
    traces = np.array([np.trace(m) / n for m in matrices])
    shifted = np.array([m - np.identity(n) * t
                        for m, t in zip(matrices, traces)])
    return shifted, traces
