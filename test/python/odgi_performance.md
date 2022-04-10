% -*- coding: utf-8 -*-

# ODGI Python performance

Here we showcase the new `odgi_ffi` module that is faster than the original `odgi` python module. Both are shipped with odgi now.

Python bindings will always be slower than C++, or Rust, because much low level copying is going on and python keeps track of typing. Still python can be useful as a tool for exploring graphs. With the newer `odgi_ffi` module we using values instead of implicit C references.

Let's time loading a graph and traversing it 100 times

```python
>>> import time
>>> from odgi_ffi import *

>>> graph = odgi_load_graph("DRB1-3123_sorted.og")

>>> odgi_get_node_count(graph)
3214

>>> res = []
>>> tic = time.perf_counter()
>>> for x in range(1, 100):
...    res.append(odgi_for_each_handle(graph, lambda h: odgi_get_sequence(graph,h) and True))

>>> toc = time.perf_counter()
>>> print(f"{toc - tic:0.4f} seconds")  # doctest: +SKIP
0.7916 seconds

```

This is twice as fast as the older `odgi` module. That would look like

```python
from odgi import *
import time

gr = graph()

gr.load("DRB1-3123_sorted.og")

gr.get_node_count()

tic = time.perf_counter()
for x in range(1, 100):
    gr.for_each_handle(lambda h: gr.get_sequence(h) and True)
toc = time.perf_counter()
print(f"{toc - tic:0.4f} seconds")

```

Run that in the odgi built source repository with something like

```sh
cd test
env LD_PRELOAD=libjemalloc.so.2 \
  LD_LIBRARY_PATH=$GUIX_ENVIRONMENT/lib:../lib \
  PYTHONPATH=../lib python3 performance2.py
1.4674 seconds
```

For more see the [odgi_ffi interface](odgi_ffi.md).
