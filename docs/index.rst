.. meta::
   :description: odgi: optimized dynamic genome/graph implementation
   :keywords: variation graph, pangenome graph, de-novo assembly

================
Welcome to odgi.
================

Representing large genomic variation graphs with minimal memory overhead requires a careful
encoding of the graph entities. It is possible to build succinct, static data structures to
store queryable graphs, as in xg, but dynamic data structures are more tricky to implement.

``odgi`` (optimized dynamic genome/graph implementation) follows the dynamic GBWT in developing
a byte-packed version of the graph and paths through it. Each node is represented by a byte
array into which we write variable length integers to represent
    1) the node sequence
    2) its edges, and
    3) the paths crossing the node.


.. toctree::
    :maxdepth: 1
    :hidden:

    rst/installation
    rst/quick_start
    rst/use_cases
    rst/commands
    rst/binding

------
Index
------

* :ref:`genindex`
* :ref:`search`