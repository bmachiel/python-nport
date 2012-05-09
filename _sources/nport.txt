:mod:`nport` --- n-ports and 2n-ports
=====================================

.. module:: nport
   :synopsis: n-ports and 2n-ports

Signal Representations
----------------------

Depending on the n-port representation (see :ref:`n-port_representations` 
below), the signals on the terminals of an n-port are represented as either
voltages and currents or traveling waves [#f1]_. The conventions for voltages,
currents and voltage waves used in :mod:`nport` are shown in figures
:ref:`2port_VI` and :ref:`2port_WA`.

.. _2port_VI:

.. figure:: images/2port_VI.*
   :align: center
   
   A two-port network descibed by voltages and currents


.. _2port_WA:

.. figure:: images/2port_WA.*
   :align: center
   
   A two-port network descibed by traveling waves


In the voltage-current representation, the currents are flowing into the ports.
In the traveling-wave representation, :math:`a_i` and :math:`b_i` represent
voltage waves flowing into and out of port :math:`i` respectively. The relation 
between voltages and currents and traveling waves is given by

.. math::

    \begin{align}
        v_i &= a_i + b_i\\
        i_i &= \frac{a_i - b_i}{Z_0}
    \end{align}


.. _n-port_representations:

n-port Representations
----------------------

Impedance (Z) parameters
^^^^^^^^^^^^^^^^^^^^^^^^

Pass :data:`Z` or :data:`IMPEDANCE` as :keyword:`type` on object
instantiation.

.. math::

    \begin{bmatrix}
        v_{1} \cr
        v_{2}
    \end{bmatrix} =
    \begin{bmatrix}
        Z_{11} & Z_{12} \cr
        Z_{21} & Z_{22}
    \end{bmatrix}
    \begin{bmatrix}
        i_{1} \cr i_{2}
    \end{bmatrix}

Admittance (Y) parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

Pass :data:`Y` or :data:`ADMITTANCE` as :keyword:`type` on object
instantiation.

.. math::

    \begin{bmatrix}
        i_{1} \cr
        i_{2}
    \end{bmatrix} =
    \begin{bmatrix}
        Y_{11} & Y_{12} \cr
        Y_{21} & Y_{22}
    \end{bmatrix}
    \begin{bmatrix}
        v_{1} \cr v_{2}
    \end{bmatrix}

Scattering (S) parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

Pass :data:`S` or :data:`SCATTERING` as :keyword:`type` on object
instantiation.

.. math::

    \begin{bmatrix}
        b_{1} \cr b_{2}
    \end{bmatrix} =
    \begin{bmatrix}
        S_{11} & S_{12} \cr S_{21} & S_{22}
    \end{bmatrix}
    \begin{bmatrix}
        a_{1} \cr a_{2}
    \end{bmatrix}


2n-port Representations
-----------------------

In addition to the :ref:`n-port_representations`, 2n-ports have two extra
representations.

Transmission (ABCD) parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pass :data:`ABCD` or :data:`TRANSMISSION` as :keyword:`type`
on object instantiation.

.. math::

    \begin{bmatrix}
        v_{1} \cr i_{1}
    \end{bmatrix} =
    \begin{bmatrix}
        A & B \cr C & D
    \end{bmatrix}
    \begin{bmatrix}
        v_{2} \cr - i_{2}
    \end{bmatrix}


Scattering Transfer (T) parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pass :data:`T` or :data:`SCATTERING_TRANSFER` as :keyword:`type`
on object instantiation.

.. math::

    \begin{bmatrix}
        a_{1} \cr b_{1}
    \end{bmatrix} =
    \begin{bmatrix}
        T_{11} & T_{12} \cr T_{21} & T_{22}
    \end{bmatrix}
    \begin{bmatrix}
        b_{2} \cr a_{2}
    \end{bmatrix}


2-port Representations
----------------------

Supported types for 2-port parameters include all of the above and,
additionally:

Hybrid (H) parameters
^^^^^^^^^^^^^^^^^^^^^

Pass :data:`H` or :data:`HYBRID` as :keyword:`type` on object
instantiation.

.. math::

    \begin{bmatrix}
        v_{1} \cr i_{2}
    \end{bmatrix} =
    \begin{bmatrix}
        H_{11} & H_{12} \cr H_{21} & H_{22}
    \end{bmatrix}
    \begin{bmatrix}
        i_{1} \cr v_{2}
    \end{bmatrix}


Inverse Hybrid (G) parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pass :data:`G` or :data:`INVERSE_HYBRID` as :keyword:`type` on
object instantiation.

.. math::

    \begin{bmatrix}
        i_{1} \cr v_{2}
    \end{bmatrix} =
    \begin{bmatrix}
        G_{11} & G_{12} \cr G_{21} & G_{22}
    \end{bmatrix}
    \begin{bmatrix}
        v_{1} \cr i_{2}
    \end{bmatrix}


Frequency-independent n-ports
-----------------------------

A Single n-port Matrix
^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: NPortMatrix
   :members:
   :undoc-members:
   :show-inheritance:


A Single 2n-port Matrix
^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: TwoNPortMatrix
   :members:
   :undoc-members:
   :show-inheritance:


Frequency-dependent n-ports
---------------------------

n-ports
^^^^^^^

.. autoclass:: NPort
   :members:
   :undoc-members:
   :show-inheritance:


2n-ports
^^^^^^^^

.. autoclass:: TwoNPort
   :members:
   :undoc-members:
   :show-inheritance:
   
   
.. [KUR65] "Power Waves and the Scattering Matrix"
           by K. Kurokawa
           in *IEEE Transactions on Microwave Theory and Techniques*,
           vol. 13, no. 2, pp. 194--202, 1965

           
.. rubric:: Footnotes

.. [#f1] Power wave representations, as introduced by Kurokawa [KUR65]_, are
         currently not supported by :mod:`nport`.
