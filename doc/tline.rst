:mod:`nport.tline` --- Transmission lines
=========================================

.. module:: nport.tline
    :synopsis: Transmission line tools


For now, only two-conductor transmission lines are supported.

Two-conductor transmission lines
--------------------------------

.. autoclass:: nport.tline.TransmissionLine
    :members: gamma, z0, rpm, lpm, gpm, cpm
    :undoc-members:

The RLGC parameters are calculated from the propagation constant and
characteristic impedance as detailed in [EIS92]_.


References
----------

.. [EIS92] "S-Parameter-Based IC Interconnect Transmission Line
           Characterization"
           by William R. Eisenstadt and Yungseon Eo
           in *IEEE Transactions on Components, Hybrids, and Manufacturing
           Technology*, 
           vol. 15, no. 4, pp.483-490, 1992
