# INSTALL ODGI

## Guix

An alternative way to manage `odgi`'s dependencies is by using the `GNU GUIX` package manager. We use Guix to develop, test and deploy odgi on our systems.
For more information see [INSTALL](./INSTALL.md).


After installing and updating guix, load the current profile with:

```sh
sudo apt install guix
guix pull
source ~/.config/guix/current/etc/profile
```

With the latest guix in the path start a GNU Guix build container with:

```bash
git clone --recursive https://github.com/pangenome/odgi.git
cd odgi
source ./.guix-build
cd build
cmake ..
make
ctest .
```

Another way of building odgi is with guix.scm:

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

For more see [guix.scm](./guix.scm).


#### installing via the guix-genomics git repository

guix genomics also contains a package definition of odgi.

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
