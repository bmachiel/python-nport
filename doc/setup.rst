.. _setup:

Getting nport
=============

Requirements
------------

The :mod:`nport` package requires at least version 2.5 of Python.

As for Python packages, :mod:`nport` depends on `NumPy and SciPy`_. It has been
tested with NumPy 1.4 and SciPy 0.8, but older versions of these will likely
work too. Try it, and let me know!

.. _NumPy and SciPy: http://www.scipy.org/

If you don't already have NumPy and SciPy, I recommend installing either the 
`Enthought Python Distribution`_ (Windows and Linux) or `Python(x,y)`_
(Windows). These provide up to date versions of Python, NumPy and SciPy together
with a bunch of other useful packages for scientific use.

.. _Enthought Python Distribution: http://www.enthought.com/
.. _Python(x,y): http://www.pythonxy.com/

Not required but very useful is `IPython`_. This is basically an enhanced 
interactive python shell. IPython is included in the Python distributions
mentioned above.

.. _IPython: http://ipython.scipy.org/


Installation
------------

The most convenient option for getting :mod:`nport` is by using `pip`_ or
`easy_install`_. To automatically download the archive from `PyPI`_ and install
it, run::

    pip install nport
    
or::

    easy_install nport

.. _pip: http://pip.openplans.org/
.. _easy_install: http://pypi.python.org/pypi/setuptools
.. _PyPI: http://pypi.python.org

If you have instead downloaded the source package, you can install it by
unpacking the archive and running::

    python setup.py install


Development
-----------

:mod:`nport` development is coordinated at 
<https://github.com/bmachiel/python-nport>. GitHub provides an `issue tracker`_
where you can report any bugs you might encounter. You will need a GitHub
account to create an issue.

.. _issue tracker: https://github.com/bmachiel/python-nport/issues

If you're interested in contributing to nport development, you can easily
`fork`_ the project and suggest any changes by sending a pull request.

.. _fork: http://help.github.com/forking/
