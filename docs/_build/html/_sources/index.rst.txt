.. meta::
   :description: odgi: optimized dynamic genome/graph implementation
   :keywords: variation graph, pangenome graph

=================================
Welcome to the odgi documentation!
=================================

In standard genomic approaches sequences are related to a single linear reference genome introducing reference bias.
`Pangenome graphs <https://pangenome.github.io/>`__ encoded in the variation graph data model describe the all versus all alignment of many sequences.
Representing large pangenome graphs with minimal memory overhead requires a careful
encoding of the graph entities. It is possible to build succinct, static data structures to
store queryable graphs, as in `xg <https://github.com/vgteam/xg>`__, but dynamic data structures are more tricky
to implement.

The optimized dynamic genome/graph implementation ``odgi`` follows the dynamic
`GBWT <https://github.com/jltsiren/gbwt>`__ in developing
a byte-packed version of the graph, edges, and paths through it. The node's id is stored as a ``uint64_t`` and its
sequence is stored as a plain ``std::string``. Bit-compressed dynamic byte arrays, with a local alphabet encoder,
represent the local neighbourhood
of the node:

    1) The node's edges, and
    2) the paths crossing the node.

To ensure minimal memory occupation, only the deltas of the neighbouring steps of a path are hold.

``odgi`` provides a set of tools ranging from graph building, manipulation, layouting, over graph statistics to graph
visualization and gene annotation lift overs.

.. toctree::
    :maxdepth: 1
    :hidden:

    Welcome <self>
    rst/installation
    rst/quick_start
    rst/tutorials
    rst/commands
    rst/binding

------
Index
------

* :ref:`genindex`
* :ref:`search`