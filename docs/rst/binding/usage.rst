######
Usage
######

.. TODO link to other pages
Make sure to install ``odgi`` (see ::ref rst/installation.rst).

To import ``odgi`` in ``Python``, make sure that the compiled ``lib/odgi.cpython*.so`` file is on your ``PYTHONPATH`` or
added to your ``Python`` path through ``sys.path.append``..

For example, assuming that your current working directory is the root of the ``odgi`` project:

.. code-block:: python

   import sys
   sys.path.append("./lib")
   import odgi
