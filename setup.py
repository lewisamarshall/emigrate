#!/usr/bin/python
from emigrate import __version__
from setuptools import setup, find_packages
try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except:
    long_description = None

setup(name='emigrate',
      version=__version__,
      author='Lewis A. Marshall',
      author_email='lewis.a.marshall@gmail.com',
      url="https://github.com/lewisamarshall/emigrate",
      classifiers=[
          "Programming Language :: Python",
          "Development Status :: 3 - Alpha",
          "Environment :: Console",
          "Intended Audience :: Science/Research",
          "Operating System :: OS Independent",
          "Topic :: Software Development :: Libraries :: Python Modules",
          "Topic :: Scientific/Engineering :: Chemistry",
          ],
      use_2to3=True,
      description='A package for simulating electrophoresis.',
      packages=find_packages(),
      long_description=long_description,
      requires=['numpy', 'scipy', 'ionize', 'h5py', 'ionize'],
      entry_points={'console_scripts': ['emigrate = emigrate.__main__:main']},
      test_suite="emigrate.tests",
      )
