.. _odgi:

#########
odgi
#########

dynamic succinct variation graph tool

SYNOPSIS
========
:ref:`odgi build` -g graph.gfa -o graph.og

:ref:`odgi validate` -i graph.og

:ref:`odgi stats` -i graph.og -S

:ref:`odgi degree` -i graph.og -S

:ref:`odgi depth` -i graph.og

:ref:`odgi overlap` -i graph.og -r
path_name

:ref:`odgi cover` -i graph.og -o
graph.paths.og

:ref:`odgi extract` -i graph.og -p
prefix -r path_name:0-17

:ref:`odgi explode` -i graph.og -p
prefix

:ref:`odgi squeeze` -f
input_graphs.txt -o graphs.og

:ref:`odgi position` -i
target_graph.og -g

:ref:`odgi sort` -i graph.og -o
graph.sorted.og -p bSnSnS

:ref:`odgi view` -i graph.og -g

:ref:`odgi kmers` -i graph.og -c -k 23
-e 34 -D 50

:ref:`odgi unitig` -i graph.og -f -t
1324 -l 120

:ref:`odgi viz` -i graph.og -o graph.og.png
-x 1920 -y 1080 -R -t 28

:ref:`odgi draw` -i graph.og -c
coords.lay -p .png -x 1920 -y 1080 -R -t 28

:ref:`odgi paths` -i graph.og -f

:ref:`odgi prune` -i graph.og -o
graph.pruned.og -c 3 -C 345 -T

:ref:`odgi unchop` -i graph.og -o
graph.unchopped.og

:ref:`odgi normalize` -i
graph.og -o graph.normalized.og -I 100 -d

:ref:`odgi bin` -i graph.og -j -w 100 -s -g

:ref:`odgi matrix` -i graph.og -e -d

:ref:`odgi chop` -i graph.og -o
graph.choped.og -c 1000

:ref:`odgi groom` -i graph.og -o
graph.groomed.og

:ref:`odgi layout` -i graph.og -o
graph.og.lay

:ref:`odgi break` -i graph.og -o
graph.broken.og -s 100 -d

:ref:`odgi pathindex` -i graph.og -o graph.xp

:ref:`odgi panpos` -i graph.og -p
Chr1 -n 4

:ref:`odgi server` -i graph.og -p
4000 -ip 192.168.8.9

:ref:`odgi test`

:ref:`odgi version`

DESCRIPTION
===========

**odgi**, the **Optimized Dynamic (genome) Graph Interface**, links a
thrifty dynamic in-memory variation graph data model to a set of
algorithms designed for scalable sorting, pruning, transformation, and
visualization of very large `genome
graphs <https://pangenome.github.io/>`__. **odgi** includes :ref:`python bindings` that can be
used to :ref:`directly interface with its data model <python bindings tutorial>`. This
**odgi** manual provides detailed information about its features and
subcommands, including examples.

COMMANDS
========

Each command has its own man page which can be viewed using e.g. **man
odgi_build.1**. Below we have a brief summary of syntax and subcommand
description.

| **odgi build** [**-g, –gfa**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
  [*OPTION*]…
| The odgi build(1) command constructs a succinct variation graph from a
  GFA. Currently, only GFA1 is supported. For details of the format
  please see https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md.

**odgi validate** [**-i, –input**\ =\ *FILE*] [*OPTION*]… The odgi
validate(1) command validates the graph (currently, it checks if the
paths are consistent with the graph topology).

| **odgi stats** [**-i, –idx**\ =\ *FILE*] [*OPTION*]…
| The odgi stats(1) command produces statistics of a variation graph.
  Among other metrics, it can calculate the #nodes, #edges, #paths and
  the total nucleotide length of the graph. Various histogram summary
  options complement the tool. If [**-B, –bed-multicov**\ =\ *BED*] is
  set, the metrics will be produced for the intervals specified in the
  BED.

**odgi degree** [**-i, –idx**\ =\ *FILE*] [*OPTION*]… The odgi degree(1)
command describes the graph in terms of node degree. For the input
graph, it shows the node.count, edge.count, avg.degree, min.degree, and
max.degree.

**odgi depth** [**-i, –input**\ =\ *FILE*] [*OPTION*]… The odgi depth(1)
command finds the depth of graph as defined by query criteria.

**odgi overlap** [**-i, –input**\ =\ *FILE*] [*OPTION*]… The odgi
overlap(1) command finds the paths touched by the input paths.

| **odgi cover** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
  [*OPTION*]…
| The odgi cover(1) command finds a path cover of a variation graph,
  with a specified number of paths per component.

**odgi extract** [**-f, –input-graphs**\ =\ *FILE*] [**-o,
–out**\ =\ *FILE*] [*OPTION*]… The odgi extract(1) command extracts
parts of the graph as defined by query criteria.

| **odgi explode** [**-i, –idx**\ =\ *FILE*] [**-p,
  –prefix**\ =\ *STRING*] [*OPTION*]…
| The odgi explode(1) command breaks a graph into connected components,
  writing each component in its own file.

**odgi squeeze** [**-f, –input-graphs**\ =\ *FILE*] [**-o,
–out**\ =\ *FILE*] [*OPTION*]… The odgi squeeze(1) command squeezes
multiple graphs into the same file.

**odgi position** [**-i, –target**\ =\ *FILE*] [*OPTION*]… The odgi
position(1) command position parts of the graph as defined by query
criteria.

| **odgi sort** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
  [*OPTION*]…
| The odgi sort(1) command sorts a succinct variation graph. The command
  offers a diverse palette of sorting algorithms to determine the node
  order:

-  A topological sort: A graph can be sorted via `breadth-first search
   (BFS) <https://en.wikipedia.org/wiki/Breadth-first_search>`__ or
   `depth-first search
   (DFS) <https://en.wikipedia.org/wiki/Depth-first_search>`__.
   Optionally, a chunk size specifies how much of the graph to grab at
   once in each topological sorting phase. The sorting algorithm will
   continue the sort from the next node in the prior graph order that
   has not been sorted, yet. The cycle breaking algorithm applies a DFS
   sort until a cycle is found. We break and start a new DFS sort phase
   from where we stopped.

-  A random sort: The graph is randomly sorted. The node order is
   randomly shuffled from `Mersenne Twister
   pseudo-random <http://www.cplusplus.com/reference/random/mt19937/>`__
   generated numbers.

-  A sparse matrix mondriaan sort: We can partition a hypergraph with
   integer weights and uniform hyperedge costs using the
   `Mondriaan <http://www.staff.science.uu.nl/~bisse101/Mondriaan/>`__
   partitioner.

-  A 1D linear SGD sort: Odgi implements a 1D linear, variation graph
   adjusted, multi-threaded version of the `Graph Drawing by Stochastic
   Gradient Descent <https://arxiv.org/abs/1710.04626>`__ algorithm. The
   force-directed graph drawing algorithm minimizes the graph’s energy
   function or stress level. It applies stochastic gradient descent
   (SGD) to move a single pair of nodes at a time.

-  An eades algorithmic sort: Use `Peter Eades’ heuristic for graph
   drawing <http://www.it.usyd.edu.au/~pead6616/old_spring_paper.pdf>`__.

Sorting the paths in a graph my refine the sorting process. For the
users’ convenience, it is possible to specify a whole pipeline of sorts
within one parameter.

| **odgi view** [**-i, –idx**\ =\ *FILE*] [*OPTION*]…
| The odgi view(1) command can convert a graph in odgi format to GFAv1.
  It can reveal a graph’s internal structures for e.g. debugging
  processes.

| **odgi kmers** [**-i, –idx**\ =\ *FILE*] [**-c, –stdout**] [*OPTION*]…
| Given a kmer length, the odgi kmers(1) command can emit all kmers. The
  output can be refined by setting the maximum number of furcations at
  edges or by not considering nodes above a given node degree limit.

| **odgi unitig** [**-i, –idx**\ =\ *FILE*] [*OPTION*]…
| The odgi unitig(1) command can print all unitigs of a given odgi graph
  to standard output in FASTA format. Unitigs can also be emitted in a
  fixed sequence quality FASTQ format. Various parameters can refine the
  unitigs to print.

| **odgi viz** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
  [*OPTION*]…
| The odgi viz(1) command can produce a linear, static visualization of
  an odgi variation graph. It aggregates the pangenome into bins and
  directly renders a raster image. The binning level depends on the
  target width of the PNG to emit. Can be used to produce visualizations
  for gigabase scale pangenomes. For more information about the binning
  process, please refer to :ref:`odgi bin`. If
  reverse coloring was selected, only the bins with a reverse rate of at
  least 0.5 are colored. Currently, there is no parameter to color
  according to the sequence coverage in bins available.

**odgi draw** [**-i, –idx**\ =\ *FILE*] [**-c, –coords-in**\ =\ *FILE*]
[**-p, –png**\ =\ *FILE*] [*OPTION*]… The odgi draw(1) command draws
previously-determined 2D layouts of the graph with diverse annotations.

| **odgi paths** [**-i, –idx**\ =\ *FILE*] [*OPTION*]…
| The odgi paths(1) command allows the investigation of paths of a given
  variation graph. It can calculate overlap statistics of groupings of
  paths.

| **odgi prune** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
  [*OPTION*]…
| The odgi prune(1) command can remove complex parts of a graph. One can
  drop paths, nodes by a certain kind of edge coverage, edges and graph
  tips. Specifying a kmer length and a maximum number of furcations, the
  graph can be broken at edges not fitting into these conditions.

| **odgi unchop** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
  [*OPTION*]…
| The odgi unchop(1) command merges each unitig into a single node.

| **odgi normalize** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
  [*OPTION*]…
| The odgi normalize(1) command
  unchops :ref:`odgi unchop` a given variation graph
  and simplifies redundant furcations.

| **odgi matrix** [**-i, –idx**\ =\ *FILE*] [*OPTION*]…
| The odgi matrix(1) command generates a sparse matrix format out of the
  graph topology of a given variation graph.

| **odgi bin** [**-i, –idx**\ =\ *FILE*] [*OPTION*]…
| The odgi bin(1) command bins a given variation graph. The pangenome
  sequence, the one-time traversal of all nodes from smallest to largest
  node identifier, can be summed up into bins of a specified size. For
  each bin, the path metainformation is summarized. This enables a
  summarized view of gigabase scale graphs. Each step of a path is a bin
  and connected to its next bin via a link. A link has a start bin
  identifier and an end bin identifier.
| The concept of odgi bin is also applied in :ref:`odgi viz`. A demonstration of how the odgi
  bin JSON output can be used for an interactive visualization is
  realized in the `Pantograph <https://graph-genome.github.io/>`__
  project. Per default, odgi bin writes the bins to stdout in a
  tab-delimited format: **path.name**, **path.prefix**, **path.suffix**,
  **bin** (bin identifier), **mean.cov** (mean coverage of the path in
  this bin), **mean.inv** (mean inversion rate of this path in this
  bin), **mean.pos** (mean nucleotide position of this path in this
  bin), **first.nucl** (first nucleotide position of this path in this
  bin), **last.nucl** (last nucleotide position of this path in this
  bin). These nucleotide ranges might span positions that are not
  present in the bin. Example: A range of 1-100 means that the first
  nucleotide has position 1 and the last has position 100, but
  nucleotide 45 could be located in another bin. For an exact positional
  output, please specify [**-j, –json**].

| **odgi chop** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
  [**-c, –chop-to**\ =\ *N*] [*OPTION*]…
| The odgi chop(1) command chops long nodes into short ones while
  preserving the graph topology.

| **odgi layout** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
  [*OPTION*]…
| The odgi layout(1) command computes 2D layouts of the graph using
  stochastic gradient descent (SGD). The input graph must be sorted and
  id-compacted. The algorithm itself is described in `Graph Drawing by
  Stochastic Gradient Descent <https://arxiv.org/abs/1710.04626>`__. The
  force-directed graph drawing algorithm minimizes the graph’s energy
  function or stress level. It applies SGD to move a single pair of
  nodes at a time.

| **odgi flatten** [**-i, –idx**\ =\ *FILE*] [*OPTION*]…
| The odgi flatten(1) command projects the graph sequence and paths into
  FASTA and BED.

| **odgi break** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
  [*OPTION*]…
| The odgi break(1) command finds cycles in a graph via `breadth-first
  search (BFS) <https://en.wikipedia.org/wiki/Breadth-first_search>`__
  and breaks them, also dropping the graph’s paths.

| **odgi pathindex** [**-i, –idx**\ =\ *FILE*] [**-o, –out**\ =\ *FILE*]
  [*OPTION*]…
| The odgi pathindex(1) command generates a path index of a graph. It
  uses succinct data structures to encode the index. The path index
  represents a subset of the features of a fully realized `xg
  index <https://github.com/vgteam/xg>`__. Having a path index, we can
  use :ref:`odgi panpos` to go from
  **path:position** → **pangenome:position** which is important when
  navigating large graphs in an interactive manner like in the
  `Pantograph <https://graph-genome.github.io/>`__ project.

| **odgi panpos** [**-i, –idx**\ =\ *FILE*] [**-p, –path**\ =\ *STRING*]
  [**-n, –nuc-pos**\ =\ *N*] [*OPTION*]…
| The odgi panpos(1) command give a pangenome position for a given path
  and nucleotide position. It requires a path index, which can be
  created with :ref:`odgi pathindex`. Going from
  **path:position** → **pangenome:position** is important when
  navigating large graphs in an interactive manner like in the
  `Pantograph <https://graph-genome.github.io/>`__ project. All input
  and output positions are 1-based.

| **odgi server** [**-i, –idx**\ =\ *FILE*] [**-p, –port**\ =\ *N*]
  [*OPTION*]…
| The odgi server(1) command starts an HTTP server with a given path
  index as input. The idea is that we can go from **path:position** →
  **pangenome:position** via GET requests to the HTTP server. The server
  headers do not block cross origin requests. Example GET request:
  *http://localost:3000/path_name/nucleotide_position*.
| The required path index can be created with :ref:`odgi pathindex`. Going from
  **path:position** → **pangenome:position** is important when
  navigating large graphs in an interactive manner like in the
  `Pantograph <https://graph-genome.github.io/>`__ project. All input
  and output positions are 1-based. If no IP address is specified, the
  server will run on localhost.

| **odgi test** [<TEST NAME|PATTERN|TAGS> …] [*OPTION*]…
| The odgi test(1) command starts all unit tests that are implemented in
  odgi. For targeted testing, a subset of tests can be selected. odgi
  test(1) depends on `Catch2 <https://github.com/catchorg/Catch2>`__. In
  the default setting, all results are printed to stdout.

| **odgi version** [*OPTION*]…
| The odgi version(1) command prints the current git version with tags
  and codename to stdout (like *v-44-g89d022b “back to old ABI”*).
  Optionally, only the release, version or codename can be printed.

BUGS
====

Refer to the **odgi** issue tracker at
https://github.com/pangenome/odgi/issues.

AUTHORS
=======

Erik Garrison from the University of California Santa Cruz wrote the
whole **odgi** tool. Simon Heumos from the Quantitative Biology Center
Tübingen wrote **odgi pathindex**, **odgi panpos**, **odgi server**, and
this documentation. Andrea Guarracino from the University of Rome Tor
Vergata wrote **odgi viz**, **odgi extract**, **odgi cover**, **odgi
explode**, **odgi squeeze**, **odgi depth**, **odgi overlap**, **odgi
validate**, and this documentation.

RESOURCES
=========

**Project web site:** https://github.com/pangenome/odgi

**Git source repository on GitHub:** https://github.com/pangenome/odgi

**GitHub organization:** https://github.com/pangenome

**Discussion list / forum:** https://github.com/pangenome/odgi/issues

COPYING
=======

The MIT License (MIT)

Copyright (c) 2019-2021 Erik Garrison

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
“Software”), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
