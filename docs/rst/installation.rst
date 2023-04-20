.. _installation:

############
Installation
############

Building from source
====================

It may be necessary to install several system-level libraries to build ``odgi``.
On ``Ubuntu``, these can be installed using ``apt``:

.. code-block:: bash

   sudo apt install build-essential cmake python3-distutils python3-dev libjemalloc-dev

Alternatively, after ``sudo apt install guix`` start a GNU Guix build container with

.. code-block:: bash

    source .guix-build

``odgi`` requires a C++ version of 9.3 or higher. You can check your version via:

.. code-block:: bash

    gcc --version
    g++ --version

If this requirement is satisfied, obtain a copy of the repository and its submodules:

.. code-block:: bash 

   git clone --recursive https://github.com/pangenome/odgi.git
   cd odgi

Finally, build ``odgi`` using ``cmake``:

.. code-block:: bash

   cmake -H. -Bbuild && cmake --build build -- -j 2

The ``-j`` argument determines the number of threads used for the compilation process. In the command above it is set to
``2``. As ``odgi`` is a fairly large project, it is recommended to set ``-j`` to the maximum number of available threads. This
can reduce the compilation time significantly.

.. note::

    If you need to avoid machine-specific optimizations, use the ``CMAKE_BUILD_TYPE=Generic`` build type:

    .. code-block:: bash

        cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Generic && cmake --build build -- -j 3


BIOCONDA
========

``odgi`` recipes for Bioconda are available at https://bioconda.github.io/recipes/odgi/README.html. To install the latest version using Conda please execute:

.. code-block:: bash

	conda install -c conda-forge -c bioconda odgi

GUIX
====

Installing via the guix-genomics git repository
-----------------------------------------------

First, clone the `guix-genomics <https://github.com/ekg/guix-genomics>`_ repository:

.. code-block:: bash

    git clone https://github.com/ekg/guix-genomics

And install the ``odgi`` package to your default GUIX environment:

.. code-block:: bash

    GUIX_PACKAGE_PATH=. guix package -i odgi

Now ODGI is available as a global binary installation.

Installing via the guix-genomics channel
----------------------------------------

Add the following to your ``~/.config/guix/channels.scm``:

.. code-block:: scm

        (cons*
      (channel
        (name 'guix-genomics)
        (url "https://github.com/ekg/guix-genomics.git")
        (branch "master"))
      %default-channels)

First, pull all the packages, then install ``odgi`` to your default GUIX environment:

.. code-block:: bash

    guix pull
    guix package -i odgi

If you want to build an environment only consisting of the ``odgi`` binary, you can do:

.. code-block:: bash

    guix environment --ad-hoc odgi

For more details about how to handle Guix channels, please go to
`https://git.genenetwork.org/guix-bioinformatics/guix-bioinformatics.git <https://git.genenetwork.org/guix-bioinformatics/guix-bioinformatics.git#headline-1>`_.

Docker
========

Thanks to our Bioconda recipe https://bioconda.github.io/recipes/odgi/README.html a docker image is generated for free.
Pulling the latest ``odgi`` docker image:

.. code-block:: bash

    docker pull quay.io/biocontainers/odgi:<tag>

Please see `odgi/tags <https://quay.io/repository/biocontainers/odgi?tab=tags>`_ for valid values for ``<tag>``.
