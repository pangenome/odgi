######
Usage
######

Install the ``odgi`` software that comes with the python module (see :ref:`installation`).

==========
PYTHONPATH
==========


To import ``odgi`` in ``Python``, make sure that the compiled ``lib/odgi.cpython*.so`` file is in your ``PYTHONPATH``.

Alternatively, add the module to the ``Python`` path through ``sys.path.append``. For example, assuming that your current working directory is the root of the ``odgi`` project, the following code should not lead to errors:

.. code-block:: python

    import sys
    sys.path.append("./lib")
    import odgi

========
Optimise
========

If python3 segfaults on

.. code-block:: bash

    env PYTHONPATH=lib python3 -c 'import odgi'

For improved performance ``odgi`` makes effective use of ``jemalloc`` for in memory data structures. One way is to preload the library, e.g.

.. code-block:: bash

    export LD_PRELOAD=/lib/x86_64-linux-gnu/libjemalloc.so.2 PYTHONPATH=lib python3 -c 'import odgi'

or

.. code-block:: bash

    env LD_PRELOAD=libjemalloc.so.2 PYTHONPATH=lib python3 -c 'import odgi'

that tells the dynamic linker to bind symbols provided by the ``jemalloc`` shared library before loading the other libraries.

==============
Debug with gdb
==============

First compile odgi with debug information. See the README.

To be able to step through odgi code load the python3 interpreter after setting a breakpoint. E.g.

.. code-block:: python

   gdb python3
   b src/odgi.cpp:143
   r
   > Python 3.9.6
   import odgi
   g = odgi.graph()
   g.load("test/DRB1-3123_sorted.og")
   g.get_node_count()

Reached breakpoint 1, odgi::graph_t::get_node_count (this=0x60f0000008b0) at src/odgi.cpp:143
