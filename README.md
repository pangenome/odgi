# odgi

[![build and test](https://github.com/pangenome/odgi/actions/workflows/build_and_test_on_push.yml/badge.svg)](https://github.com/pangenome/odgi/actions/workflows/build_and_test_on_push.yml)
[![build and python import](https://github.com/pangenome/odgi/actions/workflows/build_and_python_import_on_push.yml/badge.svg)](https://github.com/pangenome/odgi/actions/workflows/build_and_python_import_on_push.yml)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/odgi/README.html)

## optimized dynamic genome/graph implementation (odgi)

`odgi` provides an efficient and succinct dynamic DNA sequence graph model, as well as a host of algorithms that allow the use of such graphs in bioinformatic analyses.

Careful encoding of graph entities allows `odgi` to efficiently compute and transform [pangenomes](https://pangenome.github.io/) with minimal overheads.  `odgi` implements a dynamic data structure that leveraged multi-core CPUs and can be updated on the fly.

The edges and path steps are recorded as deltas between the current node id and the target node id, where the node id corresponds to the rank in the global array of nodes.
Graphs built from biological data sets tend to have local partial order and, when sorted, the deltas be small.
This allows them to be compressed with a variable length integer representation, resulting in a small in-memory footprint at the cost of packing and unpacking.

The RAM and computational savings are substantial.  In partially ordered regions of the graph, most deltas will require only a single byte.

## installation

### building from source

`odgi` requires a C++ version of 9.3 or higher. You can check your version via:

``` bash
gcc --version
g++ --version
```

`odgi` pulls in a host of source repositories as dependencies. It may be necessary to install several system-level 
libraries to build `odgi`. On `Ubuntu 20.04`, these can be installed using `apt`:

```
sudo apt install build-essential cmake python3-distutils python3-dev libjemalloc-dev
```

After installing the required dependencies, clone the `odgi` git repository recursively because of the many submodules
and build with:

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


#### Notes on dependencies

On `Arch Linux`, the `jemalloc` dependency can be installed with:

```
sudo pacman -S jemalloc     # arch linux
```

An alternative way to manage `odgi`'s dependencies is to use `GUIX`. After `sudo apt install guix`, start a GNU Guix
build container with:

```bash
source ./.guix-build
```

### Bioconda

`odgi` recipes for Bioconda are available at https://bioconda.github.io/recipes/odgi/README.html. To install the latest version using `Conda` please execute:

``` bash
conda install -c bioconda odgi
```
### Guix


#### installing via the guix-genomics git repository

First, clone the guix-genomics repository:

``` bash
git clone https://github.com/ekg/guix-genomics
```

And install the `odgi` package to your default GUIX environment:

``` bash
GUIX_PACKAGE_PATH=. guix package -i odgi
```

Now `odgi` is available as a global binary installation.

#### installing via the guix-genomics channel

Add the following to your ~/.config/guix/channels.scm:

``` scm
  (cons*
(channel
  (name 'guix-genomics)
  (url "https://github.com/ekg/guix-genomics.git")
  (branch "master"))
%default-channels)
```

First, pull all the packages, then install odgi to your default GUIX environment:

``` bash
guix pull
guix package -i odgi
```

If you want to build an environment only consisting of the odgi binary, you can do:

``` bash
guix environment --ad-hoc odgi
```

For more details about how to handle Guix channels, please go to https://git.genenetwork.org/guix-bioinformatics/guix-bioinformatics.git.

## documentation

`odgi` includes a variety of tools for analyzing and manipulating large pangenome graphs.
Read the full documentation at [https://odgi.readthedocs.io/](https://odgi.readthedocs.io/).

## multiqc

Since v1.11 [MultiQC](https://multiqc.info/) has an [ODGI module](https://multiqc.info/docs/#odgi). This module can only
work with output from `odgi stats`! For more details take a look at the documentation at [odgi.readthedocs.io/multiqc](https://odgi.readthedocs.io/en/latest/rst/multiqc.html).

## Citation
**Andrea Guarracino\*, Simon Heumos\*, Sven Nahnsen, Pjotr Prins, Erik Garrison**. [ODGI: understanding pangenome graphs](https://www.biorxiv.org/content/10.1101/2021.11.10.467921v1), bioRxiv, 2021\
**\*Shared first authorship**

## tests

Unittests from `vg` have been ported here and are used to validate the behavior of the algorithm.
They can be run via `odgi test`.

## API

`odgi::graph_t` is a `MutablePathDeletableHandleGraph` in the generic variation graph [handle graph](https://github.com/vgteam/libhandlegraph) hierarchical API model.
As such, it is possible to add, delete, and modify nodes, edges, and paths through the graph.
Wherever possible, destructive operations on the graph maintain path validity.

## versioning
Each time `odgi` is build, the current version is inferred via `git describe --always --tags`. Assuming, [version.cpp](./src/version.cpp)
is up to date, `odgi version` will not only print out the current tagged version, but its release codename, too.

## prepare release
This section is important for developers only. Each time we make a new release, we invoke [prepare_release.sh](./scripts/prepare_release.sh) (`cd` into folder [scripts](./scripts) first!)
with a new release version and codename. [version.cpp](./src/version.cpp) is updated and the documentation version is bumped up.

## presentations

[@AndreaGuarracino](https://github.com/andreaguarracino) and [@subwaystation](https://github.com/subwaystation) 
presented `odgi` at the German Bioinformatics Conference 2021: [ODGI - scalable tools for pangenome graphs](https://docs.google.com/presentation/d/1d52kaiOqeH4db4LyMHn7YNjv-mBKvhY2t2zQMNvzgno/edit#slide=id.p).

## name

`odgi` is a play on the [Italian word "oggi" (/ˈɔd.dʒi/)](https://en.wiktionary.org/wiki/oggi), which means "today".
As of 2019, a standard refrain in genomics is that genome graphs will be useful in _x_ years.
But, if we make them efficient and scalable, they will be useful today.
