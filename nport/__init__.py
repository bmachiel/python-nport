# Copyright (c) 2010- 2011 Brecht Machiels <brecht.machiels@esat.kuleuven.be>
#                          ESAT-MICAS, K.U.Leuven
#
# This file is part of python-nport (http://github.com/bmachiel/python-nport).
#
# python-nport is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python-nport is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python-nport.  If not, see <http://www.gnu.org/licenses/>.

"""python-nport

nport is a Python package for handling n-port data. It provides abstraction
for handling frequency dependent n-port parameters. It loads and saves
Touchstone and CITI files. Additionally, nport provides a set of functions
relating to RF theory, transmission lines and deembedding.
"""

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
