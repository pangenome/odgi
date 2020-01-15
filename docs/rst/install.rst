Setup
********

=================
Build directions
=================

It is straightforward to build ODGI on a unix-based machine.
First, obtain a copy of the rep and its submodules:

.. code-block:: bash 

   git clone --recursive https://github.com/vgteam/odgi.git
   cd odgi

Then build through cmake:

.. code-block:: bash

   mkdir build
   cd build
   cmake ..
   make

To make a local copy of the documentation:

.. code-block:: bash

   cd docs
   make html

================
Python Usage
================

To import ODGI in python, make sure that the compiled ``lib/odgi.cpython*.so`` file is on your python path and run ``import odgi``.

To add ``lib`` to your python path from within a script:

.. code-block:: python

   import sys
   sys.path.append("./lib")
