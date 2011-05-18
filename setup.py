#!/bin/env python

"""

setuptools setup script for nport

"""

from setuptools import setup
from subprocess import Popen, PIPE

# write the git version to nport/version.py
# based on version.py by Douglas Creager <dcreager@dcreager.net>
# http://dcreager.net/2010/02/10/setuptools-git-version-numbers/
try:
    p = Popen(['git', 'describe', '--abbrev=4'],
              stdout=PIPE, stderr=PIPE)
    p.stderr.close()
    line = p.stdout.readlines()[0]
    version = line.strip()[1:]
except OSError as e:
    print("A problem occured while trying to run git: \n" +
          e.strerror + "\n" +
          "Version information is unavailable!")
    version = 'unknown'

version_file = open('nport/version.py', 'w')
version_file.write("__version__ = '%s'\n" % version)
version_file.close()


setup(
    name='nport',
    version=version,
    packages=['nport', 'smith'],
    scripts=['nporttool'],
    requires=['numpy', 'scipy'],
    provides=['nport', 'smith'],
    test_suite='nose.collector',
    
    author="Brecht Machiels",
    author_email="brecht.machiels@esat.kuleuven.be",
    description="Python package for handling n-port data",
    url="https://github.com/bmachiel/python-nport",
    license="GPL",
    keywords="two-port 2n-port s-parameters touchstone citi deembedding smith",
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Electronic Design Automation (EDA)',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ]
)
