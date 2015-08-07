emigrate
=====
[![Code Climate](https://codeclimate.com/github/lewisamarshall/emigrate/badges/gpa.svg)](https://codeclimate.com/github/lewisamarshall/emigrate) [![Build Status](https://travis-ci.org/lewisamarshall/emigrate.svg?branch=master)](https://travis-ci.org/lewisamarshall/emigrate) [![Coverage Status](https://coveralls.io/repos/lewisamarshall/emigrate/badge.svg?branch=development&service=github)](https://coveralls.io/github/lewisamarshall/emigrate?branch=development) [![Docs Status](https://readthedocs.org/projects/emigrate/badge/?version=latest)](https://emigrate.readthedocs.org)

A Python package for performing electrophoresis calculations.

**emigrate** calculates equilibrium state and flux of ions under an electric field in a pseudo-1D domain.

The **emigrate** model uses techniques demonstrated by
[Peakmaster][peakmaster], [Spresso][Spresso], and [STEEP][STEEP]. The **emigrate**
model takes into account pH, and will take into account ionic strength, and temperature effects, including
the  most up-to-date temperature model published in STEEP. The **emigrate** object
classes make these techniques directly accessible as a backend for simulations
written in python.

Installation
------------
One-line install using [pip](https://pypi.python.org/pypi/pip):

    pip install emigrate

Tutorial
--------
Want to use **emigrate**? Read the [tutorial][tutorial], written with iPython
Notebook.


[peakmaster]: http://web.natur.cuni.cz/gas/ "Peakmaster"
[Spresso]: http://stanfordspresso.blogspot.com/ "Spresso"
[STEEP]: http://microfluidics.stanford.edu/download/ "STEEP"
[tutorial]: ./tutorial.ipynb  "emigrate Tutorial"
