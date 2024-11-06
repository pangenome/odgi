# odgi

[![build and test](https://github.com/pangenome/odgi/actions/workflows/build_and_test_on_push.yml/badge.svg)](https://github.com/pangenome/odgi/actions/workflows/build_and_test_on_push.yml)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/odgi/README.html)

## optimized dynamic genome/graph implementation (odgi)

`odgi` provides an efficient and succinct dynamic DNA sequence graph model, as well as a host of algorithms that allow the use of such graphs in bioinformatic analyses.

Careful encoding of graph entities allows `odgi` to efficiently compute and transform [pangenomes](https://pangenome.github.io/) with minimal overheads.  `odgi` implements a dynamic data structure that leveraged multi-core CPUs and can be updated on the fly.

The edges and path steps are recorded as deltas between the current node id and the target node id, where the node id corresponds to the rank in the global array of nodes.
Graphs built from biological data sets tend to have local partial order and, when sorted, the deltas be small.
This allows them to be compressed with a variable length integer representation, resulting in a small in-memory footprint at the cost of packing and unpacking.

The RAM and computational savings are substantial. In partially ordered regions of the graph, most deltas will require only a single byte.

## installation

### building from source

`odgi` requires a C++ version of 9.3 or higher. You can check your version via:

``` bash
gcc --version
g++ --version
```

`odgi` pulls in a host of source repositories as dependencies. It may be necessary to install several system-level libraries to build `odgi`. On `Ubuntu 20.04`, these can be installed using `apt`:
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

For more information on optimisations, debugging and GNU Guix builds, see [INSTALL.md](./INSTALL.md) and [CMakeLists.txt](./CMakeLists.txt).

### building with GPU

If you have GPUs and CUDA installed, you can build with GPU to use our GPU-accelerated `odgi-layout`. This will provide significant 57.3x speedup compared to the CPU solution on NVIDIA A100 GPU, reducing execution time from hours to minutes. Check out this [paper](https://arxiv.org/abs/2409.00876) and [repo](https://github.com/tonyjie/odgi) for the detailed performance speedup number. It's going to be presented at [SC'24](https://sc24.conference-program.com/presentation/?id=pap443&sess=sess382)!

Simply build with `-DUSE_GPU=ON` when cmake: 
```
cmake -DUSE_GPU=ON -H. -Bbuild && cmake --build build -- -j 3
```

To run `odgi layout` with GPU, simply add a `--gpu` with the other arguments like: 
```
odgi layout -i ${OG_FILE} -o ${LAY_FILE} --threads ${NUM_THREAD} --gpu
```


### Nix build

If you have `nix`, build and installation in your profile are as simple as:

```
nix-build && nix-env -i ./result
```

#### Notes for distribution

If you need to avoid machine-specific optimizations, use the `CMAKE_BUILD_TYPE=Generic` build type:

```shell
cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Generic && cmake --build build -- -j 3
```

#### Notes on dependencies

On `Arch Linux`, the `jemalloc` dependency can be installed with:

```
sudo pacman -S jemalloc     # arch linux
```

### Bioconda

`odgi` recipes for Bioconda are available at https://bioconda.github.io/recipes/odgi/README.html. To install the latest version using `Conda` please execute:

``` bash
conda install -c bioconda odgi
```


### Docker

To simplify installation and versioning, we have an automated GitHub action that pushes the current docker build to [dockerhub](https://hub.docker.com/r/pangenome/odgi).
To use it, pull the docker image:

```shell
docker pull pangenome/odgi
```

Then, you can run `odgi` with:

```shell
docker run odgi
```


### Guix

An alternative way to manage `odgi`'s dependencies is by using the `GNU GUIX` package manager. We use Guix to develop, test and deploy odgi on our systems.
For more information see [INSTALL](./INSTALL.md).

## FAQs

`odgi` is supported by default by [`PGGB`](https://github.com/pangenome/pggb) and part of the standard construction pipeline. Since `Cactus 2.6.1`, `odgi` has been integrated in the [`minigraph-CACTUS`](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) pangenome build with the additional option `--odgi`; this produces the native `odgi` file format for the full graph, which can later be represented as a 1D visualization or a 2D layout according to the respective parameters — for more information on this users are encouraged to look at [`help guide for minigraph-CACTUS`](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md#odgi).

Regardless, if for any reason the `minigraph-CACTUS` pangenome has been built without the additional `--odgi` option, a simple conversion step can be done to achieve the desired output.
Specifically, the new GFA v1.1 format has to be converted to a GFA v1.0 — the one used by `odgi`. The main difference between the two is whether or not including a W field to encode walks in the pangenome graph, a feature that is not as yet supported in `odgi`; hence why `minigraph-CACTUS` has to go through the following conversion step:

```
vg convert -W -t <n_of_threads> -g <pangenome_v1.1>.gfa -f > <pangenome_v1.0>.gfa
```

After converting the `minigraph-CACTUS` graph in GFA v1.0 format; it is then possible to do all conventional `odgi` analyses. For instance:

```
odgi build -g <pangenome_v1.0>.gfa -o <graph>.og -s -O -t <n_of_threads>
```

which is often useful to speed up `odgi` operations using its native format, especially when part of pipelines or called multiple times. Afterwards, paths can be extracted with `odgi paths`:

```
odgi paths -i <graph>.og -L -t <n_of_threads> | cut -f 1,2 -d '#' | sort | uniq > <prefixes>.tsv
```

and, lasty, any 1D visualization or 2D layout can be produced also for this graph. Below, an example of a 1D viz that colors haplotypes based on their ID:

```
odgi viz -i <graph>.og -o <graph_viz>.png -s '#' -M <prefixes>.tsv -t <n_of_threads>
```

## documentation

`odgi` includes a variety of tools for analyzing and manipulating large pangenome graphs.
Read the full documentation at [https://odgi.readthedocs.io/](https://odgi.readthedocs.io/).

## multiqc

Since v1.11 [MultiQC](https://multiqc.info/) has an [ODGI module](https://multiqc.info/docs/#odgi). This module can only
work with output from `odgi stats`! For more details take a look at the documentation at [odgi.readthedocs.io/multiqc](https://odgi.readthedocs.io/en/latest/rst/multiqc.html).

## Citation
**Andrea Guarracino\*, Simon Heumos\*, Sven Nahnsen, Pjotr Prins, Erik Garrison**. [ODGI: understanding pangenome graphs](https://doi.org/10.1093/bioinformatics/btac308), Bioinformatics, 2022\
**\*Shared first authorship**

**Jiajie Li, Jan-Niklas Schmelzle, Yixiao Du, Simon Heumos, Andrea Guarracino, Giulia Guidi, Pjotr Prins, Erik Garrison, Zhiru Zhang**. [Rapid GPU-Based Pangenome Graph Layout](https://arxiv.org/abs/2409.00876), SC (The International Conference for High Performance Computing, Networking, Storage, and Analysis), 2024

## funding sources

`odgi` has been funded through a variety of mechanisms, including a Wellcome Sanger PhD fellowship and diverse NIH and NSF grants (listed in our paper), as well as funding from the State of Tennessee.
Of particular note is the [contribution of NLnet to the development of a differential privacy model](https://nlnet.nl/project/VariationGraph/), ["privvg"](https://privvg.github.io/), which supported significant maintenance and development effort in the `odgi` toolkit.

## tests

Unittests from `vg` have been ported here and are used to validate the behavior of the algorithm.
They can be run via `odgi test` which is invoked by

```
ctest .
```

## API

`odgi::graph_t` is a `MutablePathDeletableHandleGraph` in the generic variation graph [handle graph](https://github.com/vgteam/libhandlegraph) hierarchical API model.
As such, it is possible to add, delete, and modify nodes, edges, and paths through the graph.
Wherever possible, destructive operations on the graph maintain path validity.

## versioning
Each time `odgi` is build, the current version is inferred via `git describe --always --tags`. Assuming, [version.cpp](./src/version.cpp)
is up to date, `odgi version` will not only print out the current tagged version, but its release codename, too.

## new release (developers only)

- Create a new release [on GitHub](https://github.com/pangenome/odgi/releases/new)
  - `Choose a tag`: v0.X.Y
  - Fill the `Release title`: ODGI v0.X.Y - Miao
  - Fill the `Describe this release` section
  - Tick `This is a pre-release`
  - Click `Publish release`
- Produce a buildable source tarball, containing code for `odgi` and all submodules, and upload it to the release.
    - Execute the following instructions:
    ```shell
    mkdir source-tarball
    cd source-tarball
    git clone --recursive https://github.com/pangenome/odgi
    cd odgi
    git fetch --tags origin
    LATEST_TAG="$(git describe --tags `git rev-list --tags --max-count=1`)"
    git checkout "${LATEST_TAG}"
    git submodule update --init --recursive
    mkdir include
    bash scripts/generate_git_version.sh include
    sed 's/execute_process(COMMAND bash/#execute_process(COMMAND bash/g' CMakeLists.txt -i
    rm -Rf .git
    find deps -name ".git" -exec rm -Rf "{}" \;
    cd ..
    mv odgi "odgi-${LATEST_TAG}"
    tar -czf "odgi-${LATEST_TAG}.tar.gz" "odgi-${LATEST_TAG}"
    rm -Rf "odgi-${LATEST_TAG}"
    ```
  - Open the (pre-)release created earlier
  - Upload the `odgi-v0.X.Y.tar.gz` file
  - Remove the tick on `This is a pre-release`
  - Click `Publish release` (this will trigger the update [on bioconda](http://bioconda.github.io/recipes/odgi/README.html))

[//]: # (This section is important for developers only. Each time we make a new release, we invoke [prepare_release.sh]&#40;./scripts/prepare_release.sh&#41; &#40;`cd` into folder [scripts]&#40;./scripts&#41; first!&#41; with a new release version and codename. [version.cpp]&#40;./src/version.cpp&#41; is updated and the documentation version is bumped up.)

## presentations

[@AndreaGuarracino](https://github.com/andreaguarracino) and [@subwaystation](https://github.com/subwaystation) presented `odgi` at the German Bioinformatics Conference 2021: [ODGI - scalable tools for pangenome graphs](https://docs.google.com/presentation/d/1d52kaiOqeH4db4LyMHn7YNjv-mBKvhY2t2zQMNvzgno/edit#slide=id.p).

## name

`odgi` is a play on the [Italian word "oggi" (/ˈɔd.dʒi/)](https://en.wiktionary.org/wiki/oggi), which means "today".
As of 2019, a standard refrain in genomics is that genome graphs will be useful in _x_ years.
But, if we make them efficient and scalable, they will be useful today.
