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

    # Z-parameters -------------------------------------------------------------
    def test_convert_z_to_z(self):
        maxerror = convert_same(self.z1)
        self.assertAlmostEqual(maxerror, 0, 15)

    def test_convert_z_to_y_to_z(self):
        maxerror = convert(self.z1, nport.Y)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_z_to_s_to_z(self):
        maxerror = convert(self.z1, nport.S)
        self.assertAlmostEqual(maxerror, 0, 13)

    # Y-parameters -------------------------------------------------------------
    def test_convert_y_to_z_to_y(self):
        maxerror = convert(self.y1, nport.Z)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_y_to_y(self):
        maxerror = convert_same(self.y1)
        self.assertAlmostEqual(maxerror, 0, 15)

    def test_convert_y_to_s_to_y(self):
        maxerror = convert(self.y1, nport.S)
        self.assertAlmostEqual(maxerror, 0, 12)
       
    # S-parameters -------------------------------------------------------------
    def test_convert_s_to_z_to_s(self):
        maxerror = convert(self.s1, nport.Z)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_s_to_y_to_s(self):
        maxerror = convert(self.s1, nport.Y)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_s_to_s(self):
        maxerror = convert_same(self.s1)
        self.assertAlmostEqual(maxerror, 0, 15)

    # Z-parameter renormalization ----------------------------------------------
    def test_convert_z_to_s_renormalize_to_z(self):
        maxerror = convert_renormalize_convert(self.z1, nport.S, 60)
        self.assertAlmostEqual(maxerror, 0, 13)

    # Y-parameter renormalization ----------------------------------------------
    def test_convert_y_to_s_renormalize_to_y(self):
        maxerror = convert_renormalize_convert(self.y1, nport.S, 60)
        self.assertAlmostEqual(maxerror, 0, 11)

    # S-parameter renormalization ----------------------------------------------
    def test_convert_s_to_z_to_s_renormalize(self):
        maxerror = convert_convert_renormalize(self.s1, nport.Z, 60)
        self.assertAlmostEqual(maxerror, 0, 13)       
       
    def test_convert_s_to_y_to_s_renormalize(self):
        maxerror = convert_convert_renormalize(self.s1, nport.Y, 60)
        self.assertAlmostEqual(maxerror, 0, 13)       

    def test_convert_s_renormalize(self):
        maxerror = convert_convert_renormalize(self.s1, nport.S, 60)
        self.assertAlmostEqual(maxerror, 0, 13)


class TestTwoPort(unittest.TestCase):
    def setUp(self):
        matrices1 = np.array([[[2.4+6.9j, 5.3+8.1j],
                               [8.4+6.4j, 6.7+2.7j]],

                              [[5.2+0.5j, 6.9+0.2j],
                               [3.1+6.7j, 8.4+6.9j]],

                              [[3.3+1.3j, 5.4+2.4j],
                               [6.7+8.4j, 1.6+8.1j]]])
        freqs = [1, 2, 3]
        self.z1 = nport.NPort(freqs, matrices1, nport.Z)
        self.y1 = nport.NPort(freqs, matrices1, nport.Y)
        self.abcd1 = nport.NPort(freqs, matrices1, nport.ABCD)
        self.h1 = nport.NPort(freqs, matrices1, nport.H)
        self.g1 = nport.NPort(freqs, matrices1, nport.G)
        self.s1 = nport.NPort(freqs, matrices1, nport.S)
        self.t1 = nport.NPort(freqs, matrices1, nport.T)

    # Z-parameters -------------------------------------------------------------
    def test_convert_z_to_z(self):
        maxerror = convert_same(self.z1)
        self.assertAlmostEqual(maxerror, 0, 15)

    def test_convert_z_to_y_to_z(self):
        maxerror = convert(self.z1, nport.Y)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_z_to_abcd_to_z(self):
        maxerror = convert(self.z1, nport.ABCD)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_z_to_h_to_z(self):
        maxerror = convert(self.z1, nport.H)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_z_to_g_to_z(self):
        maxerror = convert(self.z1, nport.G)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_z_to_s_to_z(self):
        maxerror = convert(self.z1, nport.S)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_z_to_t_to_z(self):
        maxerror = convert(self.z1, nport.T)
        self.assertAlmostEqual(maxerror, 0, 13)

    # Y-paramerers -------------------------------------------------------------
    def test_convert_y_to_z_to_y(self):
        maxerror = convert(self.y1, nport.Z)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_y_to_y(self):
        maxerror = convert_same(self.y1)
        self.assertAlmostEqual(maxerror, 0, 15)

    def test_convert_y_to_abcd_to_y(self):
        maxerror = convert(self.y1, nport.ABCD)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_y_to_h_to_y(self):
        maxerror = convert(self.y1, nport.H)
        self.assertAlmostEqual(maxerror, 0, 15)

    def test_convert_y_to_g_to_y(self):
        maxerror = convert(self.y1, nport.G)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_y_to_s_to_y(self):
        maxerror = convert(self.y1, nport.S)
        self.assertAlmostEqual(maxerror, 0, 12)

    def test_convert_y_to_t_to_y(self):
        maxerror = convert(self.y1, nport.T)
        self.assertAlmostEqual(maxerror, 0, 9)

    # ABCD-parameters ----------------------------------------------------------
    def test_convert_abcd_to_z_to_abcd(self):
        maxerror = convert(self.abcd1, nport.Z)
        self.assertAlmostEqual(maxerror, 0, 14)
       
    def test_convert_abcd_to_y_to_abcd(self):
        maxerror = convert(self.abcd1, nport.Y)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_abcd_to_abcd(self):
        maxerror = convert_same(self.abcd1)
        self.assertAlmostEqual(maxerror, 0, 15)

    def test_convert_abcd_to_h_to_abcd(self):
        maxerror = convert(self.abcd1, nport.H)
        self.assertAlmostEqual(maxerror, 0, 14)
       
    def test_convert_abcd_to_g_to_abcd(self):
        maxerror = convert(self.abcd1, nport.G)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_abcd_to_s_to_abcd(self):
        maxerror = convert(self.abcd1, nport.S)
        self.assertAlmostEqual(maxerror, 0, 11)

    def test_convert_abcd_to_t_to_abcd(self):
        maxerror = convert(self.abcd1, nport.T)
        self.assertAlmostEqual(maxerror, 0, 11)
    
    # H-parameters -------------------------------------------------------------
    def test_convert_h_to_z_to_h(self):
        maxerror = convert(self.h1, nport.Z)
        self.assertAlmostEqual(maxerror, 0, 15)

    def test_convert_h_to_y_to_h(self):
        maxerror = convert(self.h1, nport.Y)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_h_to_abcd_to_h(self):
        maxerror = convert(self.h1, nport.ABCD)
        self.assertAlmostEqual(maxerror, 0, 15)

    def test_convert_h_to_h(self):
        maxerror = convert_same(self.h1)
        self.assertAlmostEqual(maxerror, 0, 15)

    def test_convert_h_to_g_to_h(self):
        maxerror = convert(self.h1, nport.G)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_h_to_s_to_h(self):
        maxerror = convert(self.h1, nport.S)
        self.assertAlmostEqual(maxerror, 0, 15)

    def test_convert_h_to_t_to_h(self):
        maxerror = convert(self.h1, nport.T)
        self.assertAlmostEqual(maxerror, 0, 9)

    # G-parameters -------------------------------------------------------------
    def test_convert_h_to_z_to_h(self):
        maxerror = convert(self.g1, nport.Z)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_g_to_y_to_g(self):
        maxerror = convert(self.g1, nport.Y)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_h_to_abcd_to_h(self):
        maxerror = convert(self.g1, nport.ABCD)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_y_to_h_to_y(self):
        maxerror = convert(self.g1, nport.H)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_g_to_g(self):
        maxerror = convert_same(self.g1)
        self.assertAlmostEqual(maxerror, 0, 15)

    def test_convert_h_to_s_to_h(self):
        maxerror = convert(self.g1, nport.S)
        self.assertAlmostEqual(maxerror, 0, 12)

    def test_convert_h_to_t_to_h(self):
        maxerror = convert(self.g1, nport.T)
        self.assertAlmostEqual(maxerror, 0, 11)

    # S-parameters -------------------------------------------------------------
    def test_convert_s_to_z_to_s(self):
        maxerror = convert(self.s1, nport.Z)
        self.assertAlmostEqual(maxerror, 0, 13)
       
    def test_convert_s_to_y_to_s(self):
        maxerror = convert(self.s1, nport.Y)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_s_to_abcd_to_s(self):
        maxerror = convert(self.s1, nport.ABCD)
        self.assertAlmostEqual(maxerror, 0, 15)

    def test_convert_s_to_h_to_s(self):
        maxerror = convert(self.s1, nport.H)
        self.assertAlmostEqual(maxerror, 0, 13)
       
    def test_convert_s_to_g_to_s(self):
        maxerror = convert(self.s1, nport.G)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_s_to_s(self):
        maxerror = convert_same(self.s1)
        self.assertAlmostEqual(maxerror, 0, 15)

    def test_convert_s_to_t_to_s(self):
        maxerror = convert(self.s1, nport.T)
        self.assertAlmostEqual(maxerror, 0, 14)

    # T-parameters -------------------------------------------------------------
    def test_convert_t_to_z_to_t(self):
        maxerror = convert(self.t1, nport.Z)
        self.assertAlmostEqual(maxerror, 0, 14)
       
    def test_convert_t_to_y_to_t(self):
        maxerror = convert(self.t1, nport.Y)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_s_to_abcd_to_s(self):
        maxerror = convert(self.t1, nport.ABCD)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_t_to_h_to_t(self):
        maxerror = convert(self.t1, nport.H)
        self.assertAlmostEqual(maxerror, 0, 14)
       
    def test_convert_t_to_g_to_t(self):
        maxerror = convert(self.t1, nport.G)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_t_to_s_to_t(self):
        maxerror = convert(self.t1, nport.S)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_t_to_t(self):
        maxerror = convert_same(self.t1)
        self.assertAlmostEqual(maxerror, 0, 15)

    # Z-parameter renormalization ----------------------------------------------
    def test_convert_z_to_s_renormalize_to_z(self):
        maxerror = convert_renormalize_convert(self.z1, nport.S, 60)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_z_to_t_renormalize_to_z(self):
        maxerror = convert_renormalize_convert(self.z1, nport.T, 60)
        self.assertAlmostEqual(maxerror, 0, 13)

    # Y-parameter renormalization ----------------------------------------------
    def test_convert_y_to_s_renormalize_to_y(self):
        maxerror = convert_renormalize_convert(self.y1, nport.S, 60)
        self.assertAlmostEqual(maxerror, 0, 11)

    def test_convert_y_to_t_renormalize_to_y(self):
        maxerror = convert_renormalize_convert(self.y1, nport.T, 60)
        self.assertAlmostEqual(maxerror, 0, 9)

    # ABCD-parameter renormalization -------------------------------------------
    def test_convert_abcd_to_s_renormalize_to_abcd(self):
        maxerror = convert_renormalize_convert(self.abcd1, nport.S, 60)
        self.assertAlmostEqual(maxerror, 0, 11)

    def test_convert_abcd_to_t_renormalize_to_abcd(self):
        maxerror = convert_renormalize_convert(self.abcd1, nport.T, 60)
        self.assertAlmostEqual(maxerror, 0, 11)

    # S-parameter renormalization ----------------------------------------------
    def test_convert_s_to_z_to_s_renormalize(self):
        maxerror = convert_convert_renormalize(self.s1, nport.Z, 60)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_s_to_y_to_s_renormalize(self):
        maxerror = convert_convert_renormalize(self.s1, nport.Y, 60)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_s_to_abcd_to_s_renormalize(self):
        maxerror = convert_convert_renormalize(self.s1, nport.ABCD, 60)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_s_to_h_to_s_renormalize(self):
        maxerror = convert_convert_renormalize(self.s1, nport.H, 60)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_s_to_g_to_s_renormalize(self):
        maxerror = convert_convert_renormalize(self.s1, nport.G, 60)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_s_renormalize(self):
        maxerror = renormalize(self.s1, 60)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_s_to_t_to_s_renormalize(self):
        maxerror = convert_convert_renormalize(self.s1, nport.T, 60)
        self.assertAlmostEqual(maxerror, 0, 14)

    # T-parameter renormalization ----------------------------------------------
    def test_convert_t_to_z_to_t_renormalize(self):
        maxerror = convert_convert_renormalize(self.t1, nport.Z, 60)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_t_to_y_to_t_renormalize(self):
        maxerror = convert_convert_renormalize(self.t1, nport.Y, 60)
        self.assertAlmostEqual(maxerror, 0, 13)

    def test_convert_t_to_abcd_to_t_renormalize(self):
        maxerror = convert_convert_renormalize(self.t1, nport.ABCD, 60)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_t_to_h_to_t_renormalize(self):
        maxerror = convert_convert_renormalize(self.t1, nport.H, 60)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_t_to_g_to_t_renormalize(self):
        maxerror = convert_convert_renormalize(self.t1, nport.G, 60)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_t_to_s_to_t_renormalize(self):
        maxerror = convert_convert_renormalize(self.t1, nport.S, 60)
        self.assertAlmostEqual(maxerror, 0, 14)

    def test_convert_t_renormalize(self):
        maxerror = renormalize(self.t1, 60)
        self.assertAlmostEqual(maxerror, 0, 14)

    # H-parameter renormalization ----------------------------------------------
    def test_convert_h_to_s_renormalize_to_h(self):
        maxerror = convert_renormalize_convert(self.h1, nport.S, 60)
        self.assertAlmostEqual(maxerror, 0, 12)

    def test_convert_h_to_t_renormalize_to_h(self):
        maxerror = convert_renormalize_convert(self.h1, nport.T, 60)
        self.assertAlmostEqual(maxerror, 0, 11)

    # G-parameter renormalization ----------------------------------------------
    def test_convert_g_to_s_renormalize_to_g(self):
        maxerror = convert_renormalize_convert(self.g1, nport.S, 60)
        self.assertAlmostEqual(maxerror, 0, 12)

    def test_convert_g_to_t_renormalize_to_g(self):
        maxerror = convert_renormalize_convert(self.g1, nport.T, 60)
        self.assertAlmostEqual(maxerror, 0, 11)


# support methods
def convert(input, type):
    x = input.convert(type)
    input2 = x.convert(input.type)
    error = np.abs(input - input2)
    maxerror = np.max(error)
    return maxerror

def convert_same(input):
    x = input.convert(input.type)
    maxerror = np.max(np.abs(input - x))
    return maxerror

def convert_convert_renormalize(input, type, z0):
    x = input.convert(type)
    y = x.convert(input.type, z0)
    input2 = y.renormalize(input.z0)
    error = np.abs(input - input2)
    maxerror = np.max(error)
    return maxerror

def convert_renormalize_convert(input, type, z0):
    x = input.convert(type)
    y = x.renormalize(z0)
    input2 = y.convert(input.type)
    error = np.abs(input - input2)
    maxerror = np.max(error)
    return maxerror

def renormalize(input, z0):
    x = input.renormalize(z0)
    input2 = x.renormalize(input.z0)
    error = np.abs(input - input2)
    maxerror = np.max(error)
    return maxerror


if __name__ == '__main__':
    unittest.main()
