:program:`nporttool` --- Command line nport
===========================================

.. program:: nporttool

.. highlightlang:: none


The command line tool :program:`nporttool` enables to perform a subset of the
functions provided by the :mod:`nport` python package.



Functions
---------

The following are the functions supported by :program:`nporttool`.

.. _file_conversion:

File Conversion
^^^^^^^^^^^^^^^

:program:`nporttool` automatically detects the input file format. Using the 
:option:`-f` option, it is possible to specify the format of the output file.
The default output format is Touchstone. Do not specify the extension of the 
output file, it is added automatically.

.. warning::
    
    nporttool will overwrite any existing file without warning!

Example: convert the Touchstone file ``microstrip.s2p`` to CITI format and save
the result to ``microstrip.citi``::

    nporttool -f citi microstrip.s2p microstrip


.. _port_recombination:

Port Recombination
^^^^^^^^^^^^^^^^^^

It is possible to recombine some ports of the n-port described by the
S-parameters. The argument to the :option:`-r` option is a a list of
integers and pairs of integers. An integer specifies the number of the port
that needs to be kept as-is. If the integer is negative, the port's polarity
is reversed. A pair of integers specifies a pair of ports that is to be 
recombined into one port. The second port number in this pair acts  as the 
ground reference to the first port.

Example: recombine the six ports in ``example.s6p`` to four ports::

    nporttool -c "(1,3),(2,4),5,-6" example.s6p recombined

The four ports of ``recombined.s4p`` are defined so that:

* port 1 is original port 1 referenced to port 3
* port 2 is original port 2 referenced to port 4
* port 3 is the same as original port 5
* port 4 is original port 6, but with reversed polarity

If any of the input file's ports are not in the list, these ports will be left
out, practically appearing as open-circuited.

Usage
-----

Invoke :program:`nporttool` like this::

    nporttool [-c portlist] [-f format] <file>

.. cmdoption:: -f <format>
               --format <format>

    Set the output format to *format*. *format* can be one of:
    
    ``tstone``
        Touchstone
    ``citi``
        CITI
    
    .. seealso::
        :ref:`file_conversion`

.. cmdoption:: -r <portlist>
               --recombine <portlist>

    Recombine the ports in the n-port.

    .. seealso::
        :ref:`port_recombination`

.. cmdoption:: -h
               --help

   Print a short description of all command line options.
   
.. cmdoption:: -v
               --verbose

   Be verbose.