# Genome Downsampler

## Table of Contents
- [Command Line Interface (CLI) Flags](#command-line-interface-cli-flags)
    - [Required Flags](#required-flags)
    - [Optional Flags](#optional-flags)
    - [Usage Examples](#usage-examples)
- [Installation guide](#installation-guide)
    - [Dependencies](#dependencies)
      - [Common](#common)
      - [HTSlib](#htslib)
      - [OR-Tools](#or-tools)
    - [Install using precompiled binaries](#install-using-precompiled-binaries)
    - [Building from source](#building-from-source)


## Command Line Interface (CLI) Options

This section details the available CLI options for configuring the `cuda-dna-downsampler` application. Each option allows you to customize the behavior and output of the downsampling process.

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
  - Algorithm to use for downsampling. Options are: `cuda-max-flow`, `sequential-cost-scaling`, `sequential-max-flow`. Default is `sequential-max-flow`.

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
   cuda-dna-downsampler /data/input.bam 100
```

2. **Advanced usage with optional arguments:**
```sh
   cuda-dna-downsampler /data/input.bam 100 -a cuda-max-flow -v -l 100 -q 50 -p /data/filtered_out_prep.bam -o /data/output.bam -b /data/primers.bed -t /data/pairs.tsv
```

3. **Verbose mode with preprocessing output:**
```sh
   cuda-dna-downsampler /data/input.bam 100 -v -p /data/filtered_out_prep.bam -o /data/output.bam
```

4. **Using amplicon filtering:**
```sh
   cuda-dna-downsampler /data/input.bam 100 -a sequential-max-flow -v -o /data/output.bam -b /data/primers.bed -t /data/pairs.tsv

```

## Installation guide
### Dependencies
This software supports only GNU/Linux systems, if you are a Windows user, we recommend you to use WSL.
In order to compile the repository you need to have [**HTSlib**](https://github.com/samtools/htslib) and [**OR-Tools**](https://github.com/google/or-tools) libraries installed on your machine.

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
Many package managers does not provide OR-Tools library. In order to install it on your machine, download the appropriate binaries from [here](https://developers.google.com/optimization/install/cpp/binary_linux) and extract the files. 

Now, supposing that `ORTOOLS_DIR_NAME` represents the path to the extracted directory, use the following commands:

```bash
sudo cp -r ${ORTOOLS_DIR_NAME}/bin/* /usr/local/bin/
sudo cp -r ${ORTOOLS_DIR_NAME}/lib/* /usr/local/lib/
sudo cp -r ${ORTOOLS_DIR_NAME}/include/* /usr/local/include/
sudo cp -r ${ORTOOLS_DIR_NAME}/share/* /usr/local/share/
```

### Install using precompiled binaries
TODO
### Building from source 
In order to compile the source:
1. install the dependencies,
2. clone the repository, 
3. navigate to the repository directory and run `mkdir build && cd build && cmake .. && cmake --build .`,
4. the binary file location `<repository-dir>/build/src/cuda-dna-downsampler`.



