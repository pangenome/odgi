% -*- coding: utf-8 -*-

# ODGI Python FFI

The odgi toolkit for pangenomics comes with a simple "C" foreign function interface (FFI) that can be used from any computer language.
The header file for the C-API can be found [here](https://github.com/pjotrp/odgi/blob/master/src/odgi-api.h).
This C-API is covered by the `odgi-ffi` module.
In this document we walk through the low-level API using the Python `odgi_ffi` module that comes with odgi.

Note that odgi also has an older high-level Python API `import odgi` that is somewhat obsolete. Instead you should probably use below `import odgi_ffi` lower level API to construct your own library.

The main purpose of this FFI is to allow other languages to explore an in-memory odgi pangenome.

## Setting it up

First import the module. It may require setting the `PYTHONPATH` to the shared library `odgi_ffi.cpython-39-x86_64-linux-gnu.so`. Also it may be necessary to preload jemalloc with:

    env LD_PRELOAD=libjemalloc.so.2 PYTHONPATH=/odgi/lib python3 -c 'import odgi_ffi'

in a GNU Guix shell you may prepend `LD_LIBRARY_PATH` to find GLIBC etc.

    LD_LIBRARY_PATH=$LIBRARY_PATH

Now you should be able to use the `odgi_ffi` module and load the graph with 3214 nodes and 12 paths

```python
>>> from odgi_ffi import *

>>> graph = odgi_load_graph("DRB1-3123_sorted.og")
>>> odgi_get_node_count(graph)
3214

>>> odgi_get_path_count(graph)
12

>>> odgi_max_node_id(graph)
3214

```

Note that load graph reads the `.og` file format only. To read a `GFA` file convert that using the odgi tool. E.g..

```sh
    odgi build -g DRB1-3123.gfa -o DRB1-3123.og
```

## API conventions

The function names for odgi_ffi reflect the namings in the ODGI C++ library. This to avoid confusion when digging into internals.
To facilitate different language FFIs we allow for modern C++ constructs (`odgi_` prefix) or use simple C bindings (`odgi_c_` prefix). The modern C++ constructs, including unique_ptr, can be very attractive for type checking and elegant garbage collecting.
Next, by default, data is copied across language barriers. For most odgi handles this is no problem as they are long integers. But if it is text it may be useful not to copy data to speed up processing (or to modify data in place). In this case we use the `odgi_raw_` and `odgi_c_raw_` prefixes. It is ugly, but when defining the API in your favorite language there is chance of renaming or aliasing function names to get a more canonical representation. For Python, following these conventions, to get the path name we can potentially use:

```python
imort odgi_ffi

# These are examples of naming convention:
odgi_get_path_name(graph,p)       # C++ copy handler (immutable & default)
odgi_raw_get_path_name(graph,p)   # C++ string access in place (mutable)
odgi_c_get_path_name(graph,p)     # C copy handler (immutable)
odgi_c_raw_get_path_name(graph,p) # C string access in place (mutable)
```

Note that not all these functions may be implemented. See the reference (FIXME).

Since graph and p represent state we can move it into a class and end up with a nice canonical high-level interface `import odgi_graph` which we will use in the nice API document (FIXME).

In this document we only discuss the native low-level API.

## Exploring the pangenome

A pangenome is made out of nodes (sequences), edges (connectors) and paths:

![Path through pangenome](../../docs/img/exampleGraphPath.png "Pangenome path")

A path represents one version of a genome - it can represent an individual, or a version of a gene. Anything that is a linked sequence. In the picture the red path represents `CGA TTGG CCGT GT GATAA CGG ACA ATATAAC'.

In the DRB1-3123 pangenome graph that we loaded with odgi we have 12 paths.
Show the names using a call back

```python

>>> def do_path(p):
...   print(odgi_get_path_name(graph,p))
...   return True

>>> odgi_for_each_path_handle(graph, lambda p: do_path(p))
gi|568815592:32578768-32589835
gi|568815529:3998044-4011446
gi|568815551:3814534-3830133
gi|568815561:3988942-4004531
gi|568815567:3779003-3792415
gi|568815569:3979127-3993865
gi|345525392:5000-18402
gi|29124352:124254-137656
gi|28212469:126036-137103
gi|28212470:131613-146345
gi|528476637:32549024-32560088
gi|157702218:147985-163915

```

Note that we use a lambda to go through each path handle and that invokes a function because Python does not allow multi-line statements in a lambda (doh!). The `do_path` function prints the output in this case. When test returns False the process stops. This approach is not so functional, so try a more functional:

```python
>>> paths = []
>>> odgi_for_each_path_handle(graph, lambda p: paths.append(odgi_get_path_name(graph,p)))
>>> paths[0:3]
['gi|568815592:32578768-32589835', 'gi|568815529:3998044-4011446', 'gi|568815551:3814534-3830133']

```

This way we can quickly move toward a functional programming interface, see `odgi_graph`. (FIXME).

```python
>>> handles = []
>>> odgi_for_each_handle(graph, lambda h: handles.append([h,odgi_get_sequence(graph,h)]))
True
>>> handles[4:8]
[[8, 'GCTGCCATCAATGCTGGGACTTCAGGCCAA'], [10, 'TGGGAGGCAGGAAGCGTTAGGT'], [12, 'C'], [14, 'AAGATGAGG']]

```

Let's explore the edges of a node. `odgi_follow_edges` returns handles/nodes on the right and left of the edges:

```python

>>> left_edges = []
>>> right_edges = []
>>> for h in handles[4:8]:
...   res1 = odgi_follow_edges(graph,h[0],True,lambda e: right_edges.append(e))
...   res2 = odgi_follow_edges(graph,h[0],False,lambda e: left_edges.append(e))

>>> for n in right_edges:
...   [odgi_get_id(graph,n),odgi_get_sequence(graph,n)]
[7, 'C']
[4, 'G']
[17, 'TACAGATGCA']
[8, 'AAGATGAGG']
[9, 'GTTG']
[3, 'CTTGG']

>>> for n in left_edges:
...   [odgi_get_id(graph,n),odgi_get_sequence(graph,n)]
[14, 'GGGCAG']
[10, 'AGGCAT']
[11, 'AAAGGGGAGCACAAAA']
[5, 'GCTGCCATCAATGCTGGGACTTCAGGCCAA']
[7, 'C']
[4, 'G']

```

Finally we use steps to walk along a path.

```python
>>> path_name = paths[1]
>>> path_name
'gi|568815529:3998044-4011446'

>>> ph = odgi_get_path_handle(graph,path_name)
>>> step = odgi_path_begin(graph,ph)
>>> odgi_get_path_handle_of_step(graph,step) == ph
True

>>> h = odgi_get_handle_of_step(graph,step)
>>> odgi_get_sequence(graph,h)
'AT'

>>> seq = []
>>> while(odgi_has_next_step(graph,step)):
...   h = odgi_get_handle_of_step(graph,step)
...   seq.append(odgi_get_sequence(graph,h))
...   step = odgi_get_next_step(graph,step)

# Show some of the sequences for this path
>>> seq[10:14]
['CCTCT', 'T', 'GTCTCTGC', 'AGGCCACAAGCTATTATGCTTT']

# Number of steps in path
>>> len(seq)
1918

```

## Cleanup

Finally clean up with

```python
>>> odgi_free_graph(graph)

```


## Troubleshooting

### I am getting segfaults

Make sure to set `LD_PRELOAD=libjemalloc.so.2` when running Python scripts.

### I am still crashin

Try running the gdb debugger with something like

```sh
cd test
env LD_PRELOAD=libjemalloc.so.2 LD_LIBRARY_PATH=$GUIX_ENVIRONMENT/lib PYTHONPATH=../lib/ gdb --args python3 -m doctest python/odgi_ffi.md
```

and set a breakpoint on a function, e.g.

```
b odgi_path_begin
r
l
179     const step_handle_i odgi_path_begin(const ograph_t graph, path_handle_i path)
180     {
181       return as_step_handle_i(((graph_t *)graph)->path_begin(as_path_handle(path)));
182     }
p path
$1 = 2
p as_path_handle(path)
$2 = (handlegraph::path_handle_t &) @0x7fffffffcc50: {data = "\002\000\000\000\000\000\000"}
```

etc. See [gdb stepping](https://sourceware.org/gdb/current/onlinedocs/gdb/Continuing-and-Stepping.html#Continuing-and-Stepping) for more info.

### Some odgi system information

The API provides the following low level information

```python
# >>> odgi_version()
''

>>> odgi_long_long_size()
64

>>> odgi_handle_i_size()
64

>>> odgi_step_handle_i_size()
128

>>> hex((0xAB << 32) + 0xCD)
'0xab000000cd'

>>> hex((0xAB << 64) + 0xCD)
'0xab00000000000000cd'

```
