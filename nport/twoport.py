from __future__ import division

import numpy as np

from .base import Z, Y, S, T, H, G, ABCD
from .base import IMPEDANCE, ADMITTANCE, SCATTERING, SCATTERING_TRANSFER
from .base import HYBRID, INVERSE_HYBRID, TRANSMISSION 
from .nport import NPortMatrix, NPort


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
            nportmatrix = super(TwoPortMatrix, self)
            return nportmatrix.twonportmatrix().convert(type, z0).nportmatrix()

    @property
    def delta(self):
        return np.linalg.det(self)

    def stability_circle_source(self):
        if self.type != SCATTERING:
            return self.convert(SCATTERING).stability_circle_source()
        else:
            s11 = self[0, 0]
            s12 = self[0, 1]
            s21 = self[1, 0]
            s22 = self[1, 1]
            delta = self.delta
            denom = np.abs(s11)**2 - np.abs(delta)**2
            center = np.conj(s11 - delta * np.conj(s22)) / denom
            radius = np.abs(s12 * s21 / denom)
            return center, radius
        
    def stability_circle_load(self):
        if self.type != SCATTERING:
            return self.convert(SCATTERING).stability_circle_source()
        else:
            s11 = self[0, 0]
            s12 = self[0, 1]
            s21 = self[1, 0]
            s22 = self[1, 1]
            delta = self.delta
            denom = np.abs(s22)**2 - np.abs(delta)**2
            center = np.conj(s22 - delta * np.conj(s11)) / denom
            radius = np.abs(s12 * s21 / denom)
            return center, radius

    def stability_k(self):
        """Returns the Rollet stability factor K and Delta
        
        :rtype: tuple
        
        """
        if self.type != SCATTERING:
            return self.convert(SCATTERING).stability_k()
        else:
            s11 = self[0, 0]
            s12 = self[0, 1]
            s21 = self[1, 0]
            s22 = self[1, 1]
            delta = self.delta
            k = ((1 - np.abs(s11)**2 - np.abs(s22)**2 + np.abs(delta)**2) /
                 (2 * np.abs(s12 * s21)))
            return k, delta

    def stability_mu(self):
        """Returns the mu stability factor
        
        :rtype: :class:`float`
        
        """
        if self.type != SCATTERING:
            return self.convert(SCATTERING).stability_k()
        else:
            s11 = self[0, 0]
            s12 = self[0, 1]
            s21 = self[1, 0]
            s22 = self[1, 1]
            delta = self.delta
            return (1 - np.abs(s11)**2) / (np.abs(s22 - delta * np.conj(s11)) +
                                           np.abs(s12 * s21))

    def isstable_k(self):
        k, delta = self.stability_k()
        return k > 1 and np.abs(delta) < 1            
    
    def isstable_mu(self):
        return self.stability_mu() > 1


class TwoPort(NPort):
    matrix_cls = TwoPortMatrix

    def stability_circle_source(self):
        center = np.zeros(len(self), complex)
        radius = np.zeros(len(self), complex)
        for i, twoportmatrix in enumerate(self):
            center[i], radius[i] = twoportmatrix.stability_circle_source()
        return center, radius

    def stability_circle_load(self):
        center = np.zeros(len(self), complex)
        radius = np.zeros(len(self), complex)
        for i, twoportmatrix in enumerate(self):
            center[i], radius[i] = twoportmatrix.stability_circle_load()
        return center, radius

    def stability_k(self):
        k = np.zeros(len(self), complex)
        delta = np.zeros(len(self), complex)
        for i, twoportmatrix in enumerate(self):
            k[i], delta[i] = twoportmatrix.stability_k()
        return k, delta

    def stability_mu(self):
        mu = np.zeros_like(self.freqs)
        for i, twoportmatrix in enumerate(self):
            mu[i] = twoportmatrix.stability_mu()
        return mu

    def isstable_k(self):
        for twoportmatrix in self:
            if not TwoPort(twoportmatrix, self.type, self.z0).isstable_k():
                return False
        return True
    
    def isstable_mu(self):
        for twoportmatrix in self:
            if not TwoPort(twoportmatrix, self.type, self.z0).isstable_mu():
                return False
        return True
