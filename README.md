# odgi

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/odgi/README.html)

## optimized dynamic genome/graph implementation

Representing large genomic [variation graphs](https://github.com/vgteam/vg) with minimal memory overhead requires a careful encoding of the graph entities.
It is possible to build succinct, _static_ data structures to store queryable graphs, as in [xg](https://github.com/vgteam/xg), but dynamic data structures are more tricky to implement.

`odgi` follows the dynamic [GBWT](https://github.com/jltsiren/gbwt) in developing a byte-packed version of the graph and paths through it.
Each node is represented by a byte array into which we write variable length integers to represent, 1) the node sequence, 2) its edges, and 3) the paths crossing the node.

The edges and path steps are recorded relativistically, as deltas between the current node id and the target node id, where the node id corresponds to the rank in the global array of nodes.
Graphs built from biological data sets tend to have local partial order, and when sorted the stored deltas will tend to be small.
This allows them to be compressed with a variable length integer representation, resulting in a small in-memory footprint at the cost of packing and unpacking. 

The savings are substantial.  In partially ordered regions of the graph, most deltas will require only a single byte.

<!--- <> The resulting implementation is able to load the whole genome 1000 Genomes Project graph (described in our [publication on vg in Nature Biotechnology](https://www.nature.com/articles/nbt.4227)) in around 20GB of RAM.
Initially, `odgi` has been developed to allow in-memory manipulation of graphs produced by the [seqwish](https://github.com/ekg/seqwish) variation graph inducer.  -->

## building

```
git clone --recursive https://github.com/pangenome/odgi.git
cd odgi
cmake -H. -Bbuild && cmake --build build -- -j 3
```

To build a static executable, use:

```
cmake -DBUILD_STATIC=1 -H. -Bbuild && cmake --build build -- -j 3
```

You'll need to set this flag to 0 or remove and rebuild your build directory if you want to unset this build behavior and get a dynamic binary again.
Static builds are unlikely to be supported on OSX, and require appropriate static libraries on linux.

It may be necessary to install several system-level libraries to build odgi.
On Ubuntu 20.04, these can be installed using apt: `sudo apt install build-essential cmake python3-distutils python3-dev libjemalloc-dev`.

## documentation

`odgi` includes a variety of tools for analyzing and manipulating large pangenome graphs.
Read the full documentation at [https://pangenome.github.io/odgi.github.io/index.html](https://pangenome.github.io/odgi.github.io/index.html).

## tests

Unittests from `vg` have been ported here and are used to validate the behavior of the algorithm.
They can be run via `odgi test`.

## API

`odgi::graph_t` is a `MutablePathDeletableHandleGraph` in the generic variation graph [handle graph](https://github.com/vgteam/libhandlegraph) hierarchical API model.
As such, it is possible to add, delete, and modify nodes, edges, and paths through the graph.
Wherever possible, destructive operations on the graph maintain path validity.

## Documentation
There exists detailed documentation for `odgi`. For an HTML version please click [here](./docs/asciidocs/odgi_docs.html). 
For manpages please click [here](./docs/asciidocs/man).

## Versioning
Each time `odgi` is build, the current version is inferred via `git describe --always --tags`. Assuming, [version.cpp](./src/version.cpp)
is up to date, `odgi version` will not only print out the current tagged version, but its release codename, too. 

## Prepare release
This section is important for developers only. Each time we make a new release, we invoke [prepare_release.sh](./scripts/prepare_release.sh) (`cd` into folder [scripts](./scripts) first!) 
with a new release version and codename. [version.cpp](./src/version.cpp) is updated and the documentation version is bumped up.

## name

`odgi` is a play on the [Italian word "oggi" (/ˈɔd.dʒi/)](https://en.wiktionary.org/wiki/oggi), which means "today".
As of 2019, a standard refrain in genomics is that genome graphs will be useful in _x_ years.
But, if we make them efficient and scalable, they will be useful today.
