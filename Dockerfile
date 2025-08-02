FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV ORTOOLS_VERSION=or-tools_amd64_ubuntu-22.04_cpp_v9.9.3963
ENV ORTOOLS_DIR_NAME=or-tools_x86_64_Ubuntu-22.04_cpp_v9.9.3963

RUN apt-get update && apt-get install -y \
    git wget curl build-essential \
    autoconf automake make gcc perl \
    zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-gnutls-dev libssl-dev \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# HTSlib
RUN git clone --recurse-submodules https://github.com/samtools/htslib.git && \
    cd htslib && \
    make && \
    make install && \
    cd .. && rm -rf htslib

# OR-Tools
RUN wget https://github.com/google/or-tools/releases/download/v9.9/${ORTOOLS_VERSION}.tar.gz && \
    tar -xvzf ${ORTOOLS_VERSION}.tar.gz && \
    cp -r ${ORTOOLS_DIR_NAME}/bin/* /usr/local/bin/ && \
    cp -r ${ORTOOLS_DIR_NAME}/lib/* /usr/local/lib/ && \
    cp -r ${ORTOOLS_DIR_NAME}/include/* /usr/local/include/ && \
    cp -r ${ORTOOLS_DIR_NAME}/share/* /usr/local/share/ && \
    rm -rf ${ORTOOLS_VERSION}.tar.gz ${ORTOOLS_DIR_NAME}

COPY . /app
WORKDIR /app

RUN if [ -d build ]; then rm -rf build; else echo "No build directory"; fi

RUN mkdir build && cd build && \
    cmake -D WITH_CUDA=OFF .. && \
    cmake --build .

ENTRYPOINT ["/app/build/src/genome-downsampler"]


