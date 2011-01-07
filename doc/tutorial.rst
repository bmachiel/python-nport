.. _tutorial:

Tutorial
========

Basic Usage
-----------

The core of :mod:`nport` are the classes :class:`NPort` and
:class:`TwoNPort`. An :class:`NPort` object stores frequency-dependent n-port
network parameters, such as are typically stored in Touchstone (.s\ *n*\ p where
*n* is a number) or CITI files. These files are typically produced by network
analyzers and field solvers. We can load a Touchstone file into an NPort object
like this::

    >>> import nport
    >>> from nport import touchstone
    
    >>> filter = touchstone.read("filter.s2p")

Once network parameters have been read into an :class:`NPort` object, they can
be manipulated easily by using mathematical operations on it. Suppose we want to
obtain the network parameters of two of the filters connected in cascade. The
cascade connection of two 2-ports is given by the matrix multiplication of the
transfer or ABCD 2-port matrices. An ABCD matrix is a 2n-port representation, so
we first convert the filter :class:`NPort` to a :class:`TwoNPort`::

    >>> filter2n = filter.twonport().convert(nport.ABCD)
    >>> cascade2n = nport.dot(filter2n, filter2n)
    
We can now write the result back to a Touchstone file::

    >>> cascade = cascade2n.nport()
    >>> touchstone.write(cascade, "cascade")
    
If you are working with frequency-independent data, you can use the 
:class:`NPortMatrix` and :class:`TwoNPortMatrix` classes. These classes are used
internally by :class:`NPort` and :class:`TwoNPort` to store the parameter
matrices at each frequency sample. They support most of the same functions as
their frequency-dependent versions::

    >>> z = nport.NPortMatrix([[2.4+6.9j, 5.3+8.1j],
                               [8.4+6.4j, 6.7+2.7j]], nport.Z)
    >>> z.convert(nport.H)
    TwoPortMatrix(type=H)
     [[-1.93591414-6.57060176j,  1.09965504+0.76581066j],
      [-1.40973553-0.3871215j ,  0.12840169-0.05174396j]]
    

    
Extra Functionality
-------------------

Bla::

    filter.ispassive()

    
    filter.ports
    
    
    
    filter.
