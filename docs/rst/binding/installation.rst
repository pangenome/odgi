############
Installation
############

=================
Build directions
=================

First, obtain a copy of the repository and its submodules:

.. code-block:: bash 

   git clone --recursive https://github.com/pangenome/odgi.git
   cd odgi

Then build using ``cmake``:

.. code-block:: bash

   cmake -H. -Bbuild && cmake --build build -- -j 3

To make a local copy of the documentation:

.. code-block:: bash

   cd docs
   make html

================
Python Usage
================

To import ODGI in python, make sure that the compiled ``lib/odgi.cpython*.so`` file is on your `PYTHONPATH` or added to your python path through `sys.path.append` and run ``import odgi``.

For example, assuming that your current working directory is the root of the odgi project:

.. code-block:: python

   import sys
   sys.path.append("./lib")
   import odgi

