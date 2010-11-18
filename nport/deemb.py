import nport
import numpy as np


class Deembedder(object):
    def __init__(self):
        """
        Set up the de-embedder by passing the nport parameters of the dummy
        structures.
        
        """
        raise NotImplementedError
    
    def deembed(self, measurement):
        """
        De-embed the given measurement parameters

        :param measurement: full measurement structure
        :type measurement: :class:`nport.NPort`
        :returns: de-embedded parameters
        :rtype: :class:`nport.NPort` (S)
        
        """
        raise NotImplementedError


class TwoStep(Deembedder):
    """A simple two-step (open-short) de-embedder"""
    def __init__(self, open, short):
        """
        :param open: open structure
        :type open: :class:`nport.NPort`
        :param short: short structure
        :type short: :class:`nport.NPort`

        """
        y_short = short.convert(nport.Y)
        self.y_open = open.convert(nport.Y)
        self.z_so = (y_short - self.y_open).convert(nport.Z)

    def deembed(self, measurement):
        y_full = measurements.convert(nport.Y)
        z_ = (y_full - self.y_open).convert(nport.Z)
        z_dut = z_ - self.z_so
        
        return z_dut.convert(nport.S)

    deembed.__doc__ = Deembedder.deembed.__doc__


class Vandamme01(Deembedder):
    """
    "Improved Three-Step De-Embedding Method to Accurately Account for the
    Influence of Pad Parasitics in Silicon On-Wafer RF Test-Structures"
    by Ewout P. Vandamme, Dominique M. M.-P. Schreurs, and Cees van Dinther in
    *IEEE Transactions on Electron Devices*, vol. 48, no. 4, pp. 737-742, 2001
    
    """
    def __init__(self, open, short1, short2, through):
        """
        :param open: open structure
        :type open: :class:`nport.NPort`
        :param short1: port 1 short structure
        :type short1: :class:`nport.NPort`
        :param short2: port 2 short structure
        :type short2: :class:`nport.NPort`
        :param through: through structure
        :type through: :class:`nport.NPort`

        """
        # convert S to Y parameters
        y_open = open.convert(nport.Y)
        y_short1 = short1.convert(nport.Y)
        y_short2 = short2.convert(nport.Y)
        y_through = through.convert(nport.Y)

        # extract parameters from dummy structures
        y11_open = y_open.get_element(1, 1)
        y12_open = y_open.get_element(1, 2)
        y22_open = y_open.get_element(2, 2)
        y11_short1 = y_short1.get_element(1, 1)
        y22_short2 = y_short2.get_element(2, 2)
        y12_through = y_through.get_element(1, 2)

        self.g1 = y11_open + y12_open
        self.g2 = y22_open + y12_open
        self.g3 = (- y12_open.convert(nport.Z) +
                     y12_through.convert(nport.Z)).convert(nport.Y)

        z_x = y12_through.convert(nport.Z)
        z_y = (y11_short1 - self.g1).convert(nport.Z)
        z_z = (y22_short2 - self.g2).convert(nport.Z)

        self.z1 = 0.5 * (- z_x + z_y - z_z)
        self.z2 = 0.5 * (- z_x - z_y + z_z)
        self.z3 = 0.5 * (+ z_x + z_y + z_z)

    def deembed(self, measurement):
        y_meas = measurement.convert(nport.Y)
        y_a = y_meas - (np.asarray([[1, 0], [0, 0]]) * self.g1 +
                        np.asarray([[0, 0], [0, 1]]) * self.g2)
        z_a = y_a.convert(nport.Z)
        z_b = z_a - (np.asarray([[1., 0], [0, 0]]) * self.z1 +
                     np.asarray([[0, 0], [0, 1.]]) * self.z2 +
                     np.ones((2,2)) * self.z3)
        y_b = z_b.convert(nport.Y)
        y_dut = y_b - (np.asarray([[1., -1.], [-1., 1.]]) * self.g3)
        
        return y_dut.convert(nport.S)

    deembed.__doc__ = Deembedder.deembed.__doc__


kol00_asym = True

class Kolding00(Deembedder):
    """
    "A Four-Step Method for De-Embedding Gigahertz On-Wafer CMOS
    Measurements" by Troels Emil Kolding in *IEEE Transactions on Electron 
    Devices*, vol. 47, no. 4, pp. 734-740, 2000
    
    """
    def __init__(self, simple_open, simple_short, open, short1, short2, alpha=0.0):
        """        
        :param simpleopen: simple open structure
        :type simpleopen: :class:`nport.NPort`
        :param simpleshort: simple short structure
        :type simpleshort: :class:`nport.NPort`
        :param open: open structure
        :type open: :class:`nport.NPort`
        :param short1: port 1 short structure
        :type short1: :class:`nport.NPort`
        :param short2: port 2 short structure
        :type short2: :class:`nport.NPort`
        :param through: through structure
        :type through: :class:`nport.NPort`
        :param alpha: compensation parameter to account for the imperfections
                      of the short standard when gaps become large
        :type alpha: float

        """
        # convert S to Z parameters
        z_simpleopen = simple_open.convert(nport.Z)
        z_simpleshort = simple_short.convert(nport.Z)
        z_open = open.convert(nport.Z)
        z_short1 = short1.convert(nport.Z)
        z_short2 = short2.convert(nport.Z)

        # extract parameters from dummy structures
        z11_ss = z_simpleshort.get_element(1, 1)
        z22_ss = z_simpleshort.get_element(2, 2)
        z11_so = z_simpleopen.get_element(1, 1)
        z22_so = z_simpleopen.get_element(2, 2)

        # pads
        if kol00_asym:
            self.zc = 2.0/3.0 * z_simpleshort
            self.zp = (z_simpleopen - z_simpleshort)
        else:
            self.zc = 2.0/3.0 * np.identity(2) * z_simpleshort
            self.zp = np.identity(2) * (z_simpleopen - z_simpleshort)
        
        # remove pads from short1 and open
        z_short1__ = self._remove_pads(z_short1)
        z_short2__ = self._remove_pads(z_short2)
        z_open__ = self._remove_pads(z_open)

        z11_s1__ = z_short1__.get_element(1, 1)
        z21_s1__ = z_short1__.get_element(2, 1)
        z12_s1__ = z_short1__.get_element(1, 2)
        z22_s2__ = z_short2__.get_element(2, 2)
        z12_s2__ = z_short2__.get_element(1, 2)
        z21_s2__ = z_short2__.get_element(2, 1)
        z11_o__ = z_open__.get_element(1, 1)
        z21_o__ = z_open__.get_element(2, 1)
        z12_o__ = z_open__.get_element(1, 2)
        z22_o__ = z_open__.get_element(2, 2)

        # calculate remaining parameters
        if kol00_asym:
            z21 = 0.5 * (z21_s1__ + z12_s1__)
            z22 = 0.5 * (z21_s2__ + z12_s2__)
            self.z2 = np.ones((2, 2)) * (z21 + z22) * 0.5
            zi1_plus_z11 = (z11_s1__ - z21) / (1.0 + alpha)
            zi2_plus_z12 = (z22_s2__ - z22) / (1.0 + alpha)
            self.zi_plus_z1 = np.asarray([[1, 0], [0, 0]]) * zi1_plus_z11 + \
                              np.asarray([[0, 0], [0, 1]]) * zi2_plus_z12
            z31 = z21_o__ + z11_o__ - 2.0 * z21 - zi1_plus_z11
            z32 = z12_o__ + z22_o__ - 2.0 * z22 - zi2_plus_z12
            self.y3 = np.asarray([[1, 0], [0, 0]]) * z31.convert(nport.Y) + \
                      np.asarray([[0, 0], [0, 1]]) * z32.convert(nport.Y)
            zf1 = z31 * (z31 / (z21_o__ - z21) - 2)
            zf2 = z32 * (z32 / (z12_o__ - z22) - 2)
            #self.yf = np.asarray([[0, 1], [0, 0]]) * zf1.convert(nport.Y) + \
            #          np.asarray([[0, 0], [1, 0]]) * zf2.convert(nport.Y)
            self.yf = 0.5 * np.asarray([[0, 1], [1, 0]]) * \
                      (zf1 + zf2).convert(nport.Y)
        else:
            self.z2 = 0.5 * (z21_s1__ + z12_s1__)
            self.zi_plus_z1 = (z11_s1__ - self.z2) / (1.0 + alpha)
            z3 = z21_o__ + z11_o__ - 2.0 * self.z2 - self.zi_plus_z1
            zf = z3 * (z3 / (z21o__ - self.z2) - 2.0)
            self.y3 = z3.convert(nport.Y)
            self.yf = Zf.convert(nport.Y)

    def deembed(self, measurement):
        z_full = measurement.convert(nport.Z)
        z__ = self._remove_pads(z_full)

        if kol00_asym:
            print("WARNING: de-embedding method adjusted for asymmetrical fixtures")
            print("WARNING: might NOT be CORRECT")
            z___ = z__ - (self.zi_plus_z1 + self.z2)
            y___ = z___.convert(nport.Y)
            y_dut = y___ - (self.y3 - self.yf)
        else:
            z___ = z__ - (self.zi_plus_z1 * np.identity(2) +
                          self.z2 * np.ones((2,2)))
            y___ = z___.convert(nport.Y)
            y_dut = y___ - (self.y3 * np.identity(2) +
                            - self.yf * np.asarray([[0,1],[1,0]]))
        
        return y_dut.convert(nport.S)

    deembed.__doc__ = Deembedder.deembed.__doc__

    def _remove_pads(self, z_structure):
        """De-embed the pads from z_structure
        
        :param z_structure: structure to remove pads from
        :type z_structure: :class:`nport.NPort` (Z)
        :returns: structure with pads removed
        :rtype: :class:`nport.NPort` (Z)
        
        """
        assert z_structure.type == nport.Z
        z_ = z_structure - (3.0/2.0 * self.zc)
        y_ = z_.convert(nport.Y)

        y_nopads = y_ - self.zp.convert(nport.Y)
        z_nopads = y_nopads.convert(nport.Z)
        return z_nopads

