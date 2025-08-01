# Genome Downsampler

## Table of Contents
- [Genome Downsampler](#genome-downsampler)
  - [Table of Contents](#table-of-contents)
  - [Command Line Interface (CLI) Options](#command-line-interface-cli-options)
    - [Positional Arguments](#positional-arguments)
    - [Optional Arguments](#optional-arguments)
    - [Usage Examples](#usage-examples)
  - [Running in docker](#running-in-docker)
  - [Installation guide](#installation-guide)
    - [Dependencies](#dependencies)
      - [Common](#common)
      - [HTSlib](#htslib)
      - [OR-Tools](#or-tools)
      - [CUDA (optional)](#cuda-optional)
    - [Install using precompiled binaries](#install-using-precompiled-binaries)
    - [Building from source](#building-from-source)


## Command Line Interface (CLI) Options

This section details the available CLI options for configuring the `genome-downsampler` application. Each option allows you to customize the behavior and output of the downsampling process.

### Positional Arguments

- `INPUT_FILEPATH` : **Input file path**
  - Path to the input .bam file containing DNA sequences. This argument is required.

- `MAX_COVERAGE` : **Maximum coverage**
  - Maximum coverage per reference genome's base pair index. This argument is required.

### Optional Arguments

- `-h`, `--help` : **Help**
  - Print help message and exit.

- `-o`, `--output` `TEXT` : **Output file path**
  - Path to the output .bam file. Default is "output.bam" in the input file's directory.

- `-a`, `--algorithm` `TEXT` : **Algorithm**
  - Algorithm to use for downsampling. Options are: `quasi-mcp-cpu`, `quasi-mcp-cuda`, `mcp-cpu`, `qmcp-cpu`, . Default is `quasi-mcp-cpu`.

- `-b`, `--bed` `TEXT:FILE` : **BED file with amplicon bounds**
  - Path to .bed file specifying amplicon bounds for filtering or prioritization based on the selected algorithm.

- `-t`, `--tsv` `TEXT:FILE` : **TSV file with primer pairings**
  - Path to .tsv file describing pairs of primers from the .bed file to create amplicons.

- `-p`, `--preprocessing-out` `TEXT` : **Output file for preprocessed reads**
  - Path to .bam file for storing reads filtered out during preprocessing. Useful for debugging purposes.

- `-l`, `--min-length` `UINT` : **Minimum sequence length**
  - Minimum length of DNA sequences to retain. Sequences shorter than this length will be filtered out. Default is 90.

- `-q`, `--min-mapq` `UINT` : **Minimum MAPQ value**
  - Minimum Mapping Quality (MAPQ) value of sequences to retain. Sequences with MAPQ lower than this value will be filtered out. Default is 30.

- `-@`, `--threads` `UINT` : **Number of threads**
  - Number of threads for htslib read/write operations. Default is 2.

- `-v`, `--verbose` : **Verbose mode**
  - Enable additional logging for detailed execution information.

### Usage Examples

1. **Basic usage:**
```sh
   genome-downsampler /data/input.bam 100
```

2. **Advanced usage with optional arguments:**
```sh
   genome-downsampler /data/input.bam 100 -a quasi-mcp-cuda -v -l 100 -q 50 -p /data/filtered_out_prep.bam -o /data/output.bam -b /data/primers.bed -t /data/pairs.tsv
```

3. **Verbose mode with preprocessing output:**
```sh
   genome-downsampler /data/input.bam 100 -v -p /data/filtered_out_prep.bam -o /data/output.bam
```

4. **Using amplicon filtering:**
```sh
   genome-downsampler /data/input.bam 100 -v -o /data/output.bam -b /data/primers.bed -t /data/pairs.tsv

```

## Running in docker 

To run the app in `docker`, clone the repository and navigate to it:

```bash
git clone https://github.com/migoox/genome-downsampler
cd genome-downsampler
```   

Build the Docker image using:

```bash
docker build -t genome-downsampler .
```

Now you can run the container with `--help` argument:

```bash
docker run -it genome-downsampler --help 
```

To provide the data and get the output, run the app with a mounted data volume. Assuming you want to work with a file `sample.bam` located in `/home/user/data`, and you'd like the output to appear in the same folder:

```bash
docker run -it -v /home/user/data:/data genome-downsampler /data/sample.bam 100 -o /data/output.bam
```

## Installation guide
### Dependencies
This software only supports GNU/Linux systems, if you are a Windows user, we recommend using WSL. In order to run (or compile), the [**HTSlib**](https://github.com/samtools/htslib) and [**OR-Tools**](https://github.com/google/or-tools) are required to be installed on the your machine.

#### Common
Install the following common dependencies by running:

Debian/Ubuntu/Linux Mint:
```bash
sudo apt install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev
```

Fedora/Red Hat:
```bash
sudo dnf install autoconf automake make gcc perl zlib-devel bzip2-devel xz-devel libcurl-devel openssl-devel
```

OpenSUSE:
```bash
sudo zypper install autoconf automake make gcc perl zlib-devel libbz2-devel xz-devel libcurl-devel libopenssl-devel
```

Arch Linux:
```bash
sudo pacman -S autoconf automake make gcc perl zlib bzip2 xz curl openssl
```

#### HTSlib
If you prefer to install only HTSlib, you can use the following script. For the full installation guide, visit [here](https://github.com/samtools/htslib/blob/develop/INSTALL):
```bash
wget https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
tar -xf htslib-1.20.tar.bz2
cd htslib-1.20
make
sudo make install
```

Alternatively, you can install [**samtools**](http://www.htslib.org/) via package manager, since HTSlib is part of the samtools project.

Debian/Ubuntu/Linux Mint:
```bash
sudo apt install samtools
```

Fedora/Red Hat:
```bash
sudo dnf install samtools
```

OpenSUSE:
```bash
sudo zypper install samtools
```

Arch Linux:
```bash
sudo pacman -S samtools
```

#### OR-Tools
Many package managers does not provide OR-Tools library. To install it on your machine, download the appropriate binaries from [here](https://developers.google.com/optimization/install/cpp/binary_linux) and extract the files. 

Now, supposing that `ORTOOLS_DIR_NAME` represents the path to the extracted directory, use the following commands:

```bash
sudo cp -r ${ORTOOLS_DIR_NAME}/bin/* /usr/local/bin/
sudo cp -r ${ORTOOLS_DIR_NAME}/lib/* /usr/local/lib/
sudo cp -r ${ORTOOLS_DIR_NAME}/include/* /usr/local/include/
sudo cp -r ${ORTOOLS_DIR_NAME}/share/* /usr/local/share/
```

#### CUDA (optional)
The program by default tries to compile all implemented algorithms including CUDA algorithms. That means the default *build* process requires CUDA library (see the [installation guide](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html) and download cuda from [here](https://developer.nvidia.com/cuda-downloads)) and [CUDA capable GPU](https://developer.nvidia.com/cuda-gpus) to *run* the CUDA algorithms. 

Note that if you don't have CUDA capable GPU you can still download and install the CUDA library in order to run default compilation process. However to reduce the binary footprint and make the installation easier for users that may don't want to install the CUDA library, the `WITH_CUDA` flag has been provided (see the [Building from source](#building-from-source) section).

### Install using precompiled binaries
TODO
### Building from source 
1. Install the dependencies.
2. Clone the repository. 
3. Navigate to the repository directory and run 
   - `mkdir build && cd build && cmake .. && cmake --build .` for default build.
   - `mkdir build && cd build && cmake -D WITH_CUDA=OFF .. && cmake --build .` for build without CUDA.
4. The binary file location: `<repository-dir>/build/src/genome-downsampler`.

Available cmake flags:
- `WITH_CUDA`: builds the program with CUDA algorithms. When this option is disabled, the CUDA library is no longer required
- `WITH_TESTS`: builds the program with additional `test` subcommand for testing correctness of the algorithms.



