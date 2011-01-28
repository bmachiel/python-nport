import os
import numpy as np
from nport import touchstone
from nport.deemb import TwoStep, Vandamme01, Kolding00
import unittest

data_dir = os.path.split(__file__)[0] + os.path.sep + 'data' + os.path.sep


class TestDeemb(unittest.TestCase):
    def setUp(self):
        # the device under test to be extracted using the deembedding techniques
        self.dut = touchstone.read(data_dir + 'deemb_dut.s2p')
        
        # 2-port parameters for the ADS Momentum tests
        self.mom_embedded = touchstone.read(data_dir + 'deemb_mom.s2p')
        self.mom_simple_open = touchstone.read(data_dir + 'deemb_mom_simple_open.s2p')
        self.mom_simple_short = touchstone.read(data_dir + 'deemb_mom_simple_short.s2p')
        self.mom_open = touchstone.read(data_dir + 'deemb_mom_open.s2p')
        self.mom_short = touchstone.read(data_dir + 'deemb_mom_short.s2p')
        self.mom_short1 = touchstone.read(data_dir + 'deemb_mom_short1.s2p')
        self.mom_short2 = touchstone.read(data_dir + 'deemb_mom_short2.s2p')
        self.mom_through = touchstone.read(data_dir + 'deemb_mom_through.s2p')

    # Lumped test structures that have a 1-1 correspondence to the parameters
    # in the deembedding algorithms
    def test_twostep(self):
        embedded = touchstone.read(data_dir + 'deemb_twostep.s2p')
        open = touchstone.read(data_dir + 'deemb_twostep_open.s2p')
        short = touchstone.read(data_dir + 'deemb_twostep_short.s2p')
        deembedder = TwoStep(open, short)
        deembedded = deembedder.deembed(embedded)
        touchstone.write(deembedded, data_dir + 'deemb_twostep_deembedded', 'MA')
        maxerror = self.error(deembedded)
        self.assertAlmostEqual(maxerror, 0, 4)

    def test_vandamme01(self):
        embedded = touchstone.read(data_dir + 'deemb_vandamme01.s2p')
        open = touchstone.read(data_dir + 'deemb_vandamme01_open.s2p')
        short1 = touchstone.read(data_dir + 'deemb_vandamme01_short1.s2p')
        short2 = touchstone.read(data_dir + 'deemb_vandamme01_short2.s2p')
        through = touchstone.read(data_dir + 'deemb_vandamme01_through.s2p')
        deembedder = Vandamme01(open, short1, short2, through)
        deembedded = deembedder.deembed(embedded)
        touchstone.write(deembedded, data_dir + 'deemb_vandamme01_deembedded', 'MA')
        maxerror = self.error(deembedded)
        self.assertAlmostEqual(maxerror, 0, 3)

    def test_kolding00(self):
        embedded = touchstone.read(data_dir + 'deemb_kolding00.s2p')
        simple_open = touchstone.read(data_dir + 'deemb_kolding00_simple_open.s2p')
        simple_short = touchstone.read(data_dir + 'deemb_kolding00_simple_short.s2p')
        open = touchstone.read(data_dir + 'deemb_kolding00_open.s2p')
        short1 = touchstone.read(data_dir + 'deemb_kolding00_short1.s2p')
        short2 = touchstone.read(data_dir + 'deemb_kolding00_short2.s2p')
        deembedder = Kolding00(simple_open, simple_short, open, short1, short2)
        deembedded = deembedder.deembed(embedded)
        touchstone.write(deembedded, data_dir + 'deemb_kolding00_deembedded', 'MA')
        maxerror = self.error(deembedded)
        self.assertAlmostEqual(maxerror, 0, 2)

    # Deembedding from a realistic test fixture simulated in ADS Momentum
    def test_mom_twostep(self):
        deembedder = TwoStep(self.mom_open, self.mom_short)
        deembedded = deembedder.deembed(self.mom_embedded)
        touchstone.write(deembedded, data_dir + 'deemb_mom_deembedded_twostep', 'MA')
        maxerror = self.error(deembedded)
        self.assertAlmostEqual(maxerror, 0, 4)

    def test_mom_vandamme01(self):
        deembedder = Vandamme01(self.mom_open, self.mom_short1,
                                self.mom_short2, self.mom_through)
        deembedded = deembedder.deembed(self.mom_embedded)
        touchstone.write(deembedded, data_dir + 'deemb_mom_deembedded_vandamme01', 'MA')
        maxerror = self.error(deembedded)
        self.assertAlmostEqual(maxerror, 0, 4)

    def test_mom_kolding00(self):
        deembedder = Kolding00(self.mom_simple_open, self.mom_simple_short,
                               self.mom_open, self.mom_short1, self.mom_short2)
        deembedded = deembedder.deembed(self.mom_embedded)
        touchstone.write(deembedded, data_dir + 'deemb_mom_deembedded_kolding00', 'MA')
        maxerror = self.error(deembedded)
        self.assertAlmostEqual(maxerror, 0, 4)

    # support methods
    def error(self, deembedded):
        error = np.abs(self.dut - deembedded)
        maxerror = np.max(error)
        return maxerror


if __name__ == '__main__':
    unittest.main()
