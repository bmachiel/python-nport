
try:
    from version import __version__
except ImportError:
    __version__ = 'unknown (package not built using setuptools)'


from .parameter import parameter, real, imag, mag, db10, db20, rad, deg
from .base import Z, Y, S, T, H, G, ABCD
from .base import IMPEDANCE, ADMITTANCE, SCATTERING, SCATTERING_TRANSFER
from .base import HYBRID, INVERSE_HYBRID, TRANSMISSION 
from .nport import NPortMatrix, NPort, dot
from .twonport import TwoNPortMatrix, TwoNPort

from . import deemb
from . import tline
from . import touchstone
from . import citi
