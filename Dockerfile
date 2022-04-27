FROM debian:bullseye-slim AS binary

LABEL description="ODGI: understanding pangenome graphs"
LABEL base_image="debian:bullseye-slim"
LABEL software="ODGI"
LABEL about.home="https://github.com/pangenome/odgi"
LABEL about.license="SPDX:MIT"

# dependencies
RUN apt-get update \
    && apt-get install -y \
                       git \
                       bash \
                       cmake \
                       make \
                       g++ \
                       python3-dev \
                       libatomic-ops-dev \
                       autoconf \
                       libgsl-dev \
                       zlib1g-dev \
                       libzstd-dev \
                       libjemalloc-dev \
                       libhts-dev \
                       build-essential \
                       pkg-config \
  && rm -rf /var/lib/apt/lists/*

RUN git clone --recursive https://github.com/pangenome/odgi.git

RUN cd odgi \
    && cmake -H. -DCMAKE_BUILD_TYPE=Generic -Bbuild \
    && cmake --build build -- -j $(nproc) \
    && cp bin/odgi /usr/local/bin/odgi \
    && rm -rf deps \
    && rm -rf .git \
    && rm -rf build \
    && apt-get clean \
    && apt-get purge  \
    && rm -rf /var/lib/apt/lists/* /tmp/*

RUN chmod 777 /usr/local/bin/odgi

ENTRYPOINT ["odgi"]
