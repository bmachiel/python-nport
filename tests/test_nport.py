import numpy as np
import nport
import unittest

class TestNPort(unittest.TestCase):

    def setUp(self):
        matrices1 = np.array([[[ 1.2+5.2j,  2.4+6.9j,  5.3+8.1j,  7.1+2.4j],
                               [ 1.3+7.1j,  0.4+8.5j,  2.4+6.4j,  6.7+2.7j],
                               [ 8.4+3.1j,  6.4+6.7j,  8.1+2.7j,  6.4+8.4j],
                               [ 1.8+8.4j,  6.2+5.1j,  8.1+2.9j,  5.3+3.7j]],

                              [[ 5.2+0.5j,  6.9+0.2j,  8.1+6.9j,  2.4+2.4j],
                               [ 7.1+3.3j,  8.5+4.6j,  6.4+5.4j,  2.7+2.8j],
                               [ 3.1+6.7j,  6.7+8.1j,  2.7+1.6j,  8.4+6.9j],
                               [ 8.4+4.5j,  5.1+6.8j,  2.9+1.6j,  3.7+3.7j]],

                              [[ 0.5+1.2j,  0.2+2.4j,  6.9+5.3j,  2.4+7.1j],
                               [ 3.3+1.3j,  4.6+0.4j,  5.4+2.4j,  2.8+6.7j],
                               [ 6.7+8.4j,  8.1+6.4j,  1.6+8.1j,  6.9+6.4j],
                               [ 4.5+1.8j,  6.8+6.2j,  1.6+8.1j,  3.7+5.3j]]])
        freqs = [1, 2, 3]
        self.z1 = nport.NPort(freqs, matrices1, nport.Z)
        self.y1 = nport.NPort(freqs, matrices1, nport.Y)
        self.s1 = nport.NPort(freqs, matrices1, nport.S)

    def test_s_renormalize(self):
        s2 = self.s1.convert(nport.S, 60)
        s3 = s2.convert(nport.S, 50)
        maxerror = np.max(np.abs(self.s1 - s3))
        self.assertAlmostEqual(maxerror, 0, 12)

    def test_convert_z_to_z(self):
        z2 = self.z1.convert(nport.Z)
        maxerror = np.max(np.abs(self.z1 - z2))
        self.assertAlmostEqual(maxerror, 0, 12)

    def test_convert_y_to_y(self):
        y2 = self.y1.convert(nport.Y)
        maxerror = np.max(np.abs(self.y1 - y2))
        self.assertAlmostEqual(maxerror, 0, 12)

    def test_convert_s_to_s(self):
        s2 = self.s1.convert(nport.S)
        maxerror = np.max(np.abs(self.s1 - s2))
        self.assertAlmostEqual(maxerror, 0, 12)

    def test_convert_z_to_y_to_z(self):
        maxerror = self._convert_max_error(self.z1, nport.Y)
        self.assertAlmostEqual(maxerror, 0, 12)
       
    def test_convert_y_to_z_to_y(self):
        maxerror = self._convert_max_error(self.y1, nport.Z)
        self.assertAlmostEqual(maxerror, 0, 12)
       
    def test_convert_z_to_s_to_z(self):
        maxerror = self._convert_max_error(self.z1, nport.S)
        self.assertAlmostEqual(maxerror, 0, 12)
       
    def test_convert_s_to_z_to_s(self):
        maxerror = self._convert_max_error(self.s1, nport.Z)
        self.assertAlmostEqual(maxerror, 0, 12)
       
    def test_convert_y_to_s_to_y(self):
        maxerror = self._convert_max_error(self.y1, nport.S)
        self.assertAlmostEqual(maxerror, 0, 11)
       
    def test_convert_s_to_y_to_s(self):
        maxerror = self._convert_max_error(self.s1, nport.Y)
        self.assertAlmostEqual(maxerror, 0, 12)
       
    def test_s_to_z_to_s_diffz0_renormalize(self):
        z = self.s1.convert(nport.Z)
        s2 = z.convert(nport.S, 60)
        s3 = s2.renormalize(50)
        maxerror = np.max(np.abs(self.s1 - s3))
        self.assertAlmostEqual(maxerror, 0, 12)       
       
    def test_s_to_y_to_s_diffz0_renormalize(self):
        y = self.s1.convert(nport.Y)
        s2 = y.convert(nport.S, 60)
        s3 = s2.renormalize(50)
        maxerror = np.max(np.abs(self.s1 - s3))
        self.assertAlmostEqual(maxerror, 0, 12)       
       
    def _convert_max_error(self, input, type):
        x = input.convert(type)
        input2 = x.convert(input.type)
        error = np.abs(input - input2)
        maxerror = np.max(error)
        return maxerror
    

if __name__ == '__main__':
    unittest.main()
