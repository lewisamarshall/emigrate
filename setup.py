#!/usr/bin/python
from setuptools import setup, find_packages

setup(name='emigrate',
      version='0.5.0',
      author='Lewis A. Marshall',
      author_email='lewis.a.marshall@gmail.com',
      url="https://github.com/lewisamarshall/emigrate",
      classifiers=[
          "Programming Language :: Python",
          "Development Status :: 3 - Alpha",
          "Environment :: Console",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
          "Operating System :: OS Independent",
          "Topic :: Software Development :: Libraries :: Python Modules",
          "Topic :: Scientific/Engineering :: Chemistry",
          ],
      use_2to3=True,
      license='LICENSE',
      description='A package for simulating electrophoresis.',
      long_description=open('README.txt').read(),
      packages=find_packages(),
      requires=['numpy', 'scipy', 'ionize'],
      )
