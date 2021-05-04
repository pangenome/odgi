.. _installation:

############
Installation
############

.. ========
.. Building
.. ========

It may be necessary to install several system-level libraries to build ``odgi``.
On ``Ubuntu``, these can be installed using ``apt``:

.. code-block:: bash

   sudo apt install build-essential cmake python3-distutils python3-dev libjemalloc-dev

Then, obtain a copy of the repository and its submodules:

.. code-block:: bash 

   git clone --recursive https://github.com/pangenome/odgi.git
   cd odgi

Finally, build ``odgi`` using ``cmake``:

.. code-block:: bash

   cmake -H. -Bbuild && cmake --build build -- -j 2
