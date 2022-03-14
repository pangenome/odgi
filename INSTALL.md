# INSTALL ODGI

This document describes creating a build environment, development and performance tuning options.

## Creating a build environment

### Install and build with GNU Guix

We use GNU Guix to develop, test and deploy odgi on our systems. GNU Guix will run as a package manager on *any* Linux distribution.

On Debian, after installing and updating guix, load the current profile with:

```sh
sudo apt install guix
guix pull
source ~/.config/guix/current/etc/profile
```

That way, with the latest guix in the path start a GNU Guix build container with:

```bash
git clone --recursive https://github.com/pangenome/odgi.git
cd odgi
source ./.guix-build
cd build
cmake ..
make
ctest .
```

Another way of building odgi is with the provided [guix.scm](./guix.scm):

```sh
guix build -f guix.scm
```

And get a development shell with gdb:

```sh
guix shell -C -D -f guix.scm
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
cmake --build . --verbose
make
```

If you get any build issues try removing the CMakeCache files before running cmake

```sh
cd odgi
find -name CMakeCache.txt|xargs rm -v
```

For more information, see [guix.scm](./guix.scm).


#### installing via the guix-genomics git repository

The `guix genomics` git repository contains another package definition of odgi.

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

## Development & debugging

To build odgi with debug information:

```sh
cd build
cmake -DCMAKE_BUILD_TYPE ..
make
ctest . --verbose
```

## ODGI performance tuning

odgi is a high-performance computing tool and some tuning can improve speed greater than 20%.
In this section we describe the different optimization parts we explored.
But first we describe the actions to take to optimize the code output by gcc>11.

### Optimize for preformance

The default odgi build path uses position independent code (PIC) throughout to be able to use the shared libs. So

```
cd build
cmake ..
```

To build the binary without PIC makes it slightly faster

```
cmake -DPIC=OFF ..
```

Note that this only builds the `odgi` binary. For shared libs and FFI we need PIC.
Additionally adding machine flags can improve speed by 10%

```
cmake -DPIC=OFF -DEXTRA_FLAGS="-Ofast -march=native -pipe -msse4.2 -funroll-all-loops" ..
```

To make use of profiler generated optimization (PGO) compile ODGI twice(!) Run the tests (or your favorite job) after compiling with

```
cmake -DPIC=OFF -DEXTRA_FLAGS="-Ofast -march=native -pipe -msse4.2 -funroll-all-loops -fprofile-generate=../pgo" ..
make -j 8
for x in 1 2 3 ; do ctest . ; done
```

With the profiling data in ../pgo we can now ask gcc to optimize compiler output with

```
cmake -DPIC=OFF -DEXTRA_FLAGS="-Ofast -march=native -pipe -msse4.2 -funroll-all-loops -fprofile-use=../pgo" ..
make -j 8
for x in 1 2 3 ; do ctest . ; done
```

This should gain the odgi binary another 10% of speed. To optimize the shared libs run both PGO steps with PIC.

### Position independent code

Normally compilation with position independent code (`-fPIC` option for `gcc`) is  detrimental to performance. In general this is true for odgi too, so to build the binary tool we default.

* On epysode AMD EPYC 3251 8-Core Processor the current master runs at slightly over 25.0s.
* Disabling PIC on odgi-obj drops to 24.6s

### Native compilation

With native compilation we gain a little on the AMD EPYC 3251 8-Core Processor

* Native compilation march=native at 24.5s

### Extra flags

Tested on AMD EPYC 3251 8-Core Processor

* -mavx2 has impact (down to 22.6s)
* -msse4.2 has impact (down to 22.4s)
* -funroll-all-loops has a little impact

The following flags have no impact

* -pipe does not help runtime performance, but helps compile time
* -fomit-frame-pointer
* -ffast-math
* -march=native
* -fforce-addr
* -D_GLIBCXX_PARALLEL

Current recommended build flags

```
cmake -DEXTRA_FLAGS="-Ofast -march=native -pipe -msse4.2 -funroll-all-loops" ..
```

result in CMake generating

```
CXX_FLAGS = -fopenmp -Ofast -march=native -pipe -msse4.2 -funroll-all-loops -fopenmp -fopenmp -DNDEBUG -std=gnu++17
```

### Profile guided optimizations

* Run profile guided optimizations, e.g. https://stackoverflow.com/questions/14492436/g-optimization-beyond-o3-ofast
  + 190s run of odgi test to generate profile
  + PGO goes at 20s. So the overall gain is greater than 20%.

### Static library (libodgi.a) and shared library (libodgi.so)

The optimal build of odgi is with `-DPIC=OFF` that removes the switch `-fPIC` an can not build shared libs.
Therefore the binary always builds with library `libodgi.a`.

To build the shared library use the default, i.e., `-DPIC=ON`. Still, the odgi binary won't use that.
This may change in the future.

### Looking at the libs

After all this the CMake flags are not honoured by libhandlegraph. It builds with

```
CXX_FLAGS =  -O3 -g -fPIC -std=gnu++14
```

Compiling with inlined sources, however, made no dent in performance.

Above metrics were using sdsl-lite which also does not honour global flags with

```
CXX_FLAGS =  -std=c++11 -Wall -Wextra -DNDEBUG -O3 -DNDEBUG
```

this build is also at 22.4s. Let's try adding our flags with the following patch

```patch
diff --git a/CMakeLists.txt b/CMakeLists.txt
index 657949f..97a6863 100644
--- a/CMakeLists.txt
+++ b/CMakeLists.txt
@@ -43,7 +43,10 @@ else()
   message(STATUS "Compiler is recent enough to support C++11.")
 endif()

-if( CMAKE_COMPILER_IS_GNUCXX )
+set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG")
+set (CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -Ofast -pipe -msse4.2 -funroll-all-loops")
+
+if( CMAKE_COMPILER_IS_GNUCXXU )
     append_cxx_compiler_flags("-std=c++11 -Wall -Wextra -DNDEBUG" "GCC" CMAKE_CXX_FLAGS)
     append_cxx_compiler_flags("-O3 -ffast-math -funroll-loops" "GCC" CMAKE_CXX_OPT_FLAGS)
     if ( CODE_COVERAGE )
```

resulting in

```
CXX_FLAGS =  -Ofast -pipe -msse4.2 -funroll-all-loops -DNDEBUG
```

improved speed from 22.7s to 22.4s. That points out that these switches help AND that the native sdsl-lite libs that come with a distro are optimally compiled.

I am not sure what the impact can be of the other libs, but we will find that out after profiling the code.
