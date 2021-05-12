######
Usage
######

First, make sure you have installed ``odgi`` (see :ref:`installation`).

Then, export the following

.. code-block:: bash

    export LD_PRELOAD=/lib/x86_64-linux-gnu/libjemalloc.so.2

to tell the dynamic linker to bind symbols provided by the ``jemalloc`` shared library before other libraries.

To import ``odgi`` in ``Python``, make sure that the compiled ``lib/odgi.cpython*.so`` file is in your ``PYTHONPATH`` or
added to your ``Python`` path through ``sys.path.append``. For example, assuming that your current working directory is
the root of the ``odgi`` project, the following code should not lead to errors:

.. code-block:: python

    import sys
    sys.path.append("./lib")
    import odgi
