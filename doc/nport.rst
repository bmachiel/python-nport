:mod:`nport` --- n-ports and 2n-ports
=====================================

.. module:: nport
   :synopsis: n-ports and 2n-ports

Convention for voltages, currents and voltage waves.

.. image:: images/2port_VI.*

.. image:: images/2port_WA.*

Supported types for n-port parameters are:

* Z = IMPEDANCE = TYPE_Z = 'Z'

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

* Y = ADMITTANCE = TYPE_Y = 'Y'

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

* S = SCATTERING = TYPE_S = 'S'

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

Supported types for 2n-port parameters are:

* T = SCATTERING_TRANSFER = TYPE_T = 'T'
* ABCD = TRANSMISSION = TYPE_ABCD = 'ABCD'

.. math::

    \begin{bmatrix}
        v_{1} \cr i_{1}
    \end{bmatrix} =
    \begin{bmatrix}
        A & B \cr C & D
    \end{bmatrix}
    \begin{bmatrix}
        v_{2} \cr i_{2}
    \end{bmatrix}




Supported types for 2-port parameters include all of the above and,
additionally:

* H = HYBRID = TYPE_H = 'H'
* G = INVERSE_HYBRID = TYPE_G = 'G'


Classes
-------

.. autoclass:: nport.NPort
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: nport.NPortMatrix
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: nport.TwoNPort
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: nport.TwoNPortMatrix
   :members:
   :undoc-members:
   :show-inheritance:


