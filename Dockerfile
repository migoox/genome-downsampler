FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    git curl build-essential sudo \
    autoconf automake \
    zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-gnutls-dev libssl-dev cmake \
    && rm -rf /var/lib/apt/lists/*

COPY . /app
WORKDIR /app
RUN ./scripts/install_libs.sh install

RUN if [ -d build ]; then rm -rf build; else echo "No build directory"; fi

RUN mkdir build && cd build && \
    cmake --preset gcc-x64-release .. && \
    cmake --build --preset gcc-x64-release .

ENTRYPOINT ["/app/build/src/genome-downsampler"]


