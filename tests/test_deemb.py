import os
import numpy as np
from nport import touchstone
from nport.deemb import TwoStep, Vandamme01, Kolding00
import unittest

data_dir = os.path.split(__file__)[0] + os.path.sep + 'data' + os.path.sep


class TestDeemb(unittest.TestCase):
    def setUp(self):
        self.dut = touchstone.read(data_dir + 'deemb_dut.s2p')

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

    # support methods
    def error(self, deembedded):
        error = np.abs(self.dut - deembedded)
        maxerror = np.max(error)
        return maxerror


if __name__ == '__main__':
    unittest.main()
