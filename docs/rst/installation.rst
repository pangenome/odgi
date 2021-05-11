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

``odgi`` requires a C++ version of 9.3 or higher. You can check your version via:

.. code-block:: bash

    gcc --version
    g++ --version

If this requirement is satisfied, obtain a copy of the repository and its submodules:

.. code-block:: bash 

   git clone --recursive https://github.com/pangenome/odgi.git
   cd odgi

Finally, build ``odgi`` using ``cmake``:

.. code-block:: bash

   cmake -H. -Bbuild && cmake --build build -- -j 2

The ``-j`` argument determines the number of threads used for the compilation process. In the command above it is set to
``2``. As ``odgi`` is a fairly large project, it is recommended to set ``-j`` to the maximum number of available threads. This
can reduce the compilation time significantly.
