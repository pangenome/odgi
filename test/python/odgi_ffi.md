# ODGI Python FFI

The odgi toolkit for pangenomics comes with a simple "C" foreign function interface (FFI) that can be used from any computer language.
The header file for the C-API can be found [here](https://github.com/pjotrp/odgi/blob/master/src/odgi-api.h).
In this document we walk through the low-level API using the Python `odgi_ffi` module that comes with odgi.

## Setting it up

First import the module. It may require setting the `PYTHONPATH` to the shared library `odgi_ffi.cpython-39-x86_64-linux-gnu.so`. Also it may be necessary to preload jemalloc with:

    env LD_PRELOAD=libjemalloc.so.2 PYTHONPATH=/odgi/lib python3 -c 'import odgi_ffi'

in a GNU Guix shell prepend to find GLIBC etc.

    LD_LIBRARY_PATH=$LIBRARY_PATH

Now you should be able to use the `odgi_ffi` module and load the graph with 3214 nodes and 12 paths

```python
>>> from odgi_ffi import *

>>> odgi_version()
'f937271'

>>> graph = odgi_load_graph("DRB1-3123_sorted.og")
>>> odgi_get_node_count(graph)
3214

>>> odgi_get_path_count(graph)
12

```
