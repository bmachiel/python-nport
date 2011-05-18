import numpy as np


def parameter(real=None, imag=None, mag=None, db10=None, db20=None, deg=None,
              rad=None):
    """Initialize a parameter

    Specify:
    * real and (optionally) imag, or
    * mag/db and deg/rad

    """
    if real is not None:
        if (mag or db10 or db20 or deg or rad) is not None:
            raise ValueError('Illegal combination of arguments.')
        if imag is None:
            imag = 0
    else: # real is None
        if ((imag is not None) and
            (mag and db10 and db20 is not None) and 
            (deg and rad is not None)):
            raise ValueError('Illegal combination of arguments.')

        if db10 is not None:
            mag = np.power(10, db10 / 10.0)
        elif db20 is not None:
            mag = np.power(10, db20 / 20.0)
        if deg is not None:
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


def db10(arg):
    """Return the magnitude in decibels (power) of the elements of the array"""
    return 10.0 * np.log10(mag(arg))


def db20(arg):
    """Return the magnitude in decibels (V or I) of the elements of the array"""
    return 20.0 * np.log10(mag(arg))


def rad(arg):
    """Return the phase in radians of the elements of the array"""
    return np.angle(arg)


def deg(arg):
    """Return the phase in degrees of the elements of the array"""
    return np.angle(arg, deg=True)
