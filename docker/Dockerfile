FROM debian:buster-slim AS compile
RUN apt-get update \
    && apt-get install -y \
                       git \
                       bash \
                       cmake \
                       make \
                       g++ \
                       python3-dev \
    && rm -rf /var/lib/apt/lists/*
RUN git clone --depth=1 --recursive https://github.com/vgteam/odgi.git
RUN cd odgi \
    && cmake -H. -Bbuild \
    && cmake --build build --config Release -- -j $(nproc) \
    && cd build \
    && make install \
    && strip -o /usr/local/bin/odgi-stripped /usr/local/bin/odgi 

FROM debian:buster-slim AS debug
RUN apt-get update \
    && apt-get install -y libgomp1 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=compile /usr/local/bin/odgi /usr/local/bin/
LABEL base_image="debian:buster-slim"
LABEL software="odgi-debug"
LABEL about.home="https://github.com/vgteam/odgi"
LABEL about.license="SPDX:MIT"
ENTRYPOINT ["odgi"]

FROM debian:buster-slim AS binary
RUN apt-get update \
    && apt-get install -y libgomp1 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=compile /usr/local/bin/odgi-stripped /usr/local/bin/odgi
LABEL base_image="debian:buster-slim"
LABEL software="odgi"
LABEL about.home="https://github.com/vgteam/odgi"
LABEL about.license="SPDX:MIT"
ENTRYPOINT ["odgi"]


FROM debian:buster-slim AS python3
RUN apt-get update \
    && apt-get install -y libgomp1 python3 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=compile /usr/local/lib/odgi*.so /usr/local/lib/python3.7/dist-packages/
COPY --from=compile /usr/local/lib/libodgi.so /usr/local/lib/python3.7/dist-packages/
CMD ["bash"]
