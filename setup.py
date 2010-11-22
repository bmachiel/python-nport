#!/bin/env python

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
except:
    print("A problem occured while trying to run git. "
          "Version information is unavailable!")
    version = 'unknown'

version_file = open('nport/version.py', 'w')
version_file.write("__version__ = '%s'\n" % version)
version_file.close()


setup(
    name='nport',
    version=version,
    packages=['nport'],
    scripts=['nporttool'],
    requires=['numpy','scipy'],
    
    author="Brecht Machiels",
    author_email="brecht.machiels@esat.kuleuven.be",
    description="Python package for handling n-port data",
    license="GPL",
    keywords="two-port 2n-port s-parameters touchstone citi deembedding",
    url="https://github.com/bmachiel/python-nport",
)
