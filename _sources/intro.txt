.. _intro:

Introduction
============

In circuit theory, n-port network parameters are widely used in modeling both
lumped and distributed networks. Often, the networks are frequency-dependent and
are represented by a list of network parameters corresponding to a list of
sampled frequency points. This data is typically stored in Touchstone or CITI
files.

:mod:`nport` tries to provide an intuitive interface to n-port network
parameters. The frequency dependence of the parameters is handled transparantly
so that the user does not have to take care of this manually. This way,
performing calculations on n-port data is simplified significantly. In addition,
a number of functions are provided to manipulate the n-port data, such as
converting between different n-port representations.

Next to providing an abstraction for n-port parameter data, :mod:`nport`
includes extra modules building on top of this basic functionality:

* reading and writing Touchstone and CITI files
* extraction of transmission line parameters from 2-port data
* RF de-embedding algorithms

Refer to the :ref:`tutorial` for a more in-depth introduction to using
:mod:`nport`.
