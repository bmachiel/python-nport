:mod:`nport.deemb` --- De-embedding tools
=========================================

.. module:: nport.deemb
   :synopsis: De-embedding tools

Introduction
------------

This module implements a number of de-embedding methods. Each de-embedding
method is implemented as a subclass of :class:`nport.deemb.Deembedder`.
Currently, only a few de-embedding methods are included with :mod:`nport` (see 
:ref:`de-embedding_methods` below).

Usage
-----

A :class:`Deembedder` object is instantiated by passing the n-port parameters of 
the method's dummy structures. The resulting :class:`Deembedder` object can then
de-embed any measurement:

.. automethod:: Deembedder.deembed

The following example de-embeds some transistor measurements using the two-step
de-embedding method::

    from nport import touchstone
    from nport.deemb import TwoStep

    # read in S-parameters of de-embedding structures
    open_ = touchstone.read("deemb_open.s2p")
    short = touchstone.read("deemb_short.s2p")

    # set up the de-embedder
    twostep = TwoStep(open_, short)

    # read in S-parameters of devices
    nmos16 = touchstone.read("nmos_W16_L80_Vg10_Vd12.s2p")
    nmos32 = touchstone.read("nmos_W32_L80_Vg10_Vd12.s2p")
    nmos64 = touchstone.read("nmos_W64_L80_Vg10_Vd12.s2p")

    # de-embed the measurements
    deemb16 = twostep.deembed(nmos16)
    deemb32 = twostep.deembed(nmos32)
    deemb64 = twostep.deembed(nmos64)

    # write the de-embedded S-parameters to touchstone files
    touchstone.write(deemb16, "nmos_W16_L80_Vg10_Vd12_deembedded")
    touchstone.write(deemb32, "nmos_W32_L80_Vg10_Vd12_deembedded")
    touchstone.write(deemb64, "nmos_W64_L80_Vg10_Vd12_deembedded")
    
Other de-embedding methods require more dummy structures. Aside from having to
pass more arguments when instantiating the :class:`Deembedder` object, the
de-embedding script looks the same.

.. _de-embedding_methods:

De-embedding methods
--------------------

Two-step de-embedding
^^^^^^^^^^^^^^^^^^^^^

This class implements basic open-short deembedding.

.. autoclass:: nport.deemb.TwoStep

Vandamme three-step de-embedding
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This class implements the three-step deembedding algorithm by Vandamme et al.
[VAN01]_.

.. autoclass:: nport.deemb.Vandamme01

Kolding four-step de-embedding
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This class implements the three-step deembedding algorithm by Kolding [KOL00]_.

.. autoclass:: nport.deemb.Kolding00


.. [VAN01] "Improved Three-Step De-Embedding Method to Accurately Account for
           the Influence of Pad Parasitics in Silicon On-Wafer RF
           Test-Structures"
           by Ewout P. Vandamme, Dominique M. M.-P. Schreurs,
           and Cees van Dinther
           in *IEEE Transactions on Electron Devices*,
           vol. 48, no. 4, pp. 737-742, 2001

.. [KOL00] "A Four-Step Method for De-Embedding Gigahertz On-Wafer CMOS
           Measurements"
           by Troels Emil Kolding
           in *IEEE Transactions on Electron Devices*,
           vol. 47, no. 4, pp. 734-740, 2000
