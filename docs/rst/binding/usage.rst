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

For improved performance ``odgi`` makes effective use of ``jemalloc`` for in memory data structures. One way is to preload the library, e.g.

.. code-block:: bash

    export LD_PRELOAD=/lib/x86_64-linux-gnu/libjemalloc.so.2

that tells the dynamic linker to bind symbols provided by the ``jemalloc`` shared library before loading the other libraries.
