import numpy as np


def parameter(real=None, imag=None, mag=None, db=None, deg=None, rad=None):
    """Initialize a parameter

    Specify:
    * real and (optionally) imag, or
    * mag/db and deg/rad

    """
    if real is not None:
        assert (mag is None) and (db is None) and (deg is None) and (rad is None)
        if imag is None:
            imag = 0
    else: # real is None
        assert imag is None
        if mag is not None:
            assert db is None
        if db is not None:
            assert mag is None
            mag = np.pow(10, db/20.0)
        if rad is not None:
            assert deg is None
        if deg is not None:
            assert rad is None
            rad = np.radians(deg)
        real = mag * np.cos(rad)
        imag = mag * np.sin(rad)
    return complex(real, imag)

def real(arg):
    """Return the real part of the elements of the array"""
    return np.real(arg)

def imag(arg):
    """Return the imaginary part of the elements of the array"""
    return np.imag(arg)

def mag(arg):
    """Return the magnitude of the elements of the array"""
    return np.abs(arg)

def db(arg):
    """Return the magnitude in decibels of the elements of the array"""
    return 20.0 * np.log10(mag(arg))

def rad(arg):
    """Return the phase in radians of the elements of the array"""
    return np.angle(arg)

def deg(arg):
    """Return the phase in degrees of the elements of the array"""
    return np.angle(arg, deg=True)
