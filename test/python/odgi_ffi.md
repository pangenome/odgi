% -*- coding: utf-8 -*-

# ODGI Python FFI

The odgi toolkit for pangenomics comes with a simple "C" foreign function interface (FFI) that can be used from any computer language.
The header file for the C-API can be found [here](https://github.com/pjotrp/odgi/blob/master/src/odgi-api.h).
In this document we walk through the low-level API using the Python `odgi_ffi` module that comes with odgi.

Note that odgi also has a higher level Python API. You may also decide to use this low level API to construct your own library.
Feel free to go either way!

The main purpose of this FFI is to allow other languages to explore an in-memory odgi pangenome.

## Setting it up

First import the module. It may require setting the `PYTHONPATH` to the shared library `odgi_ffi.cpython-39-x86_64-linux-gnu.so`. Also it may be necessary to preload jemalloc with:

    env LD_PRELOAD=libjemalloc.so.2 PYTHONPATH=/odgi/lib python3 -c 'import odgi_ffi'

in a GNU Guix shell prepend to find GLIBC etc.

    LD_LIBRARY_PATH=$LIBRARY_PATH

Now you should be able to use the `odgi_ffi` module and load the graph with 3214 nodes and 12 paths

```python
>>> from odgi_ffi import *

>>> graph = odgi_load_graph("DRB1-3123_sorted.og")
>>> odgi_get_node_count(graph)
3214

>>> odgi_get_path_count(graph)
12

```

Note that load graph reads the `.og` file format only. To read a `GFA` file convert that using the odgi tool. E.g..

    odgi build -g DRB1-3123.gfa -o DRB1-3123.og

## Exploring the pagenome

A pangenome is made out of nodes (sequences), edges (connectors) and paths:

![Path through pangenome](../../docs/img/exampleGraphPath.png "Pangenome path")

A path represents one version of a genome - it can represent an individual, or a version of a gene. Anything that is a linked sequence. In the picture the red path represents `CGA TTGG CCGT GT GATAA CGG ACA ATATAAC'.

In the pangenome graph that we loaded with odgi we have 12 paths.
Show the names using a call back


```python

>>> import codecs

>>> def test(p):
...   odgi_get_path_name(graph,p)
...   return True

>>> list = []
>>> odgi_for_each_path_handle2(graph, lambda p: test(p))

```


path_name5 = nil
ODGI.each_path(pangenome) { |path|
  p [path,ODGI.path_name(pangenome,path)]
  path_name5 = ODGI.path_name(pangenome,path) if path==5
}

raise "No path[5]" if not ODGI.has_path(pangenome,path_name5)

path5 = ODGI::get_path(pangenome,path_name5);
p [path5,ODGI.path_name(pangenome,path5)]

ODGI.each_handle(pangenome) { |handle|
  right_edges = []
  ODGI::follow_edges(pangenome,handle,true) { |edge|
    right_edges.push edge
  }
  left_edges = []
  ODGI::follow_edges(pangenome,handle,false) { |edge|
    left_edges.push edge
  }
  p [handle,ODGI::get_id(pangenome,handle),ODGI::get_sequence(pangenome,handle),left_edges,right_edges]
}
