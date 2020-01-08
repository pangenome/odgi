Setup
********

=================
Build directions
=================

.. code-block:: bash
   
    mkdir build
    cd build
    cmake ..
    make

To make a local copy of the documentation

.. code-block:: bash
   
   cd docs
   make html

================
Python Usage directions
================

To use ODGI in python, make sure that lib/odgi.cpython*.so file is on your python path and run ``import odgi``.

You can add /lib to your path within python as follows

.. code-block:: python

   import sys
   sys.append.path("./lib")
