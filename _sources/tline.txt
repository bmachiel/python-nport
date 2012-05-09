:mod:`nport.tline` --- Transmission lines
=========================================

.. module:: nport.tline
    :synopsis: Transmission line tools


A distinction is made between simple two-conductor transmission lines and 
multiconductor transmission lines. The latter can be analyzed in terms of 
propagating modes.


Two-conductor transmission lines
--------------------------------

.. autoclass:: TransmissionLine
    :members: gamma, z0, rpm, lpm, gpm, cpm
    :undoc-members:

The RLGC parameters are calculated from the propagation constant and
characteristic impedance as detailed in [EIS92]_.



Multi-conductor transmission lines
----------------------------------

.. todo:: complete class interface and documentation

.. autoclass:: nport.tline.MulticonductorTransmissionLine
    :members:
    :undoc-members:

The implementation of :class:`MulticonductorTransmissionLine` follows [FAR04]_. 
An example of the application of this technique can be found in [FAR04b]_. The
technique was generalized in order to support non-reciprocal transmission lines.

More information on the analysis of (uniform) multiconductor transmission lines 
can be found in [FAR93]_ and [PAU08]_.


.. [EIS92] "S-Parameter-Based IC Interconnect Transmission Line
           Characterization"
           by William R. Eisenstadt and Yungseon Eo
           in *IEEE Transactions on Components, Hybrids, and Manufacturing
           Technology*, 
           vol. 15, no. 4, pp.483-490, 1992

.. [FAR04] "A new generalized modal analysis theory for nonuniform
           multiconductor transmission lines"
           by J.A. Brandao Faria, 
           in *IEEE Transactions on Power Systems*, 
           vol. 19, no. 2, pp.926-933, 2004

.. [FAR04b] "A New Modal Analysis Theory for Multiconductor Nonuniform
            Transmission-Line Structures: Application to the Analysis of Line 
            Junctions"
            by J.A. Brandao Faria, 
            in *IEEE Transactions on Power Systems*, 
            vol. 19, no. 3, pp.1380-1386, 2004


.. [FAR93] "Multiconductor Transmission-Line Structures - Modal Analysis
           Techniques"
           by J.A. Brandao Faria, 
           Wiley-Interscience, 1993


.. [PAU08] "Analysis of Multiconductor Transmission Lines", 2nd edition
           by Clayton R. Paul
           Wiley-IEEE, 2008
