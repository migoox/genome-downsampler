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


## Command Line Interface (CLI) Flags

This section details the available CLI flags for configuring the application. Each flag allows you to customize the behavior and output of the cuda-dna-downsampler.

### Required Flags
The following flags must be specified for the cuda-dna-downsampler to run:

- `-i`, `--input` : **Input file path**
  - This flag is required to specify the path to the input file.
  - Example: `-i /path/to/input/file`

- `-a`, `--algorithm` : **Algorithm ID (descriptive string)**
  - This flag is required to specify which algorithm to use. Currently available algorithms are:
    - test - test algorithm just to make sure CLI works properly
    - sequential-cost-scaling - deprecated algorithm of MinCostMaxFlow with cost scaling
    - sequential-max-flow - max flow algorithm based on ortools library
    - cuda-max-flow - CUDA implementation of max flow push-relabel algorithm based on this [paper](https://www.sciencedirect.com/science/article/pii/B9780123859631000058) 
  - Example: `-a "test"`

- `-M`, `--max-coverage` : **Maximum coverage**
  - This flag is required to specify the maximum coverage of reference genome index.
  - Example: `-M 100`

### Optional Flags

- `-v`, `--verbose` : **Verbose mode**
  - Enables verbose mode. The definition of verbosity is determined by each algorithm.
  - Example: `-v`

- `-l`, `--min-length` : **Minimum length of sequence**
  - Sets the minimum length of a sequence. Sequences shorter than this will be filtered out along with their paired sequences before algorithm execution (in the importer). Default is 90.
  - Example: `-l 100`

- `-q`, `--min-mapq` : **Minimum MAPQ of sequence**
  - Sets the minimum MAPQ (Mapping Quality) of a sequence. Sequences with a MAPQ lower than this will be filtered out along with their paired sequences before algorithm execution (in the importer). Default is 30.
  - Example: `-q 50`

- `-o`, `--output` : **Output file path**
  - Sets the path for the output file. If not specified, the default is `output.bam` in the input file's directory.
  - Example: `-o /path/to/output.bam`

- `-p`, `--preprocessing-out` : **Output for reads filtered out during preprocessing file path**
  - Turn on and sets the path for the writing reads which were filtered during preprocessing of input file. It can be useful for debugging.
  - Example: `-p /path/to/preprocessing_filtered_out_output.bam`

- `-b`, `--bed` : **BED file with primers file path**
  - Turn on filtering/grading reads by single-amplicon inclusion based on primers (amplicon bounds) from specified file. When no additional .tsv file is specified it takes every two primers (first with second, third with fourth etc.) as single amplicon bounds.
  - Example: `-b /path/to/primers.bed`

- `-t`, `--tsv` : **TSV file with primers pairing file path**
  - Should be always used with (-b|--bed) flag. It specifies which primers should be paired to create the amplicon.
  - Example: `-t /path/to/pairs.tsv`

- `-@`, `--threads` : **Number of threads**
  - Number of threads for compression/decompression during reading/writing .bam files. Default is set to 2.
  - Example: `-@ 2`

### Usage Examples

1. **Basic usage with required flags:**
   ```sh
   cuda-dna-downsampler -i /data/input.bam -a "sequential-max-flow" -M 100
   ```

2. **Using all optional flags:**
   ```sh
   cuda-dna-downsampler -i /data/input.bam -a "sequential-max-flow" -M 100 -v -l 100 -q 50 -p /data/filtered_out_prep.bam -o /data/output.bam -b /data/primers.bed -t /data/pairs.tsv
   ```

3. **Verbose mode with preprocessing data saved and output path:**
   ```sh
   cuda-dna-downsampler -i /data/input.bam -a "sequential-max-flow" -M 100 -v -p /data/filtered_out_prep.bam -o /data/output.bam
   ```

4. **Verbose mode with amplicon filtering and output path:**
   ```sh
   cuda-dna-downsampler -i /data/input.bam -a "sequential-max-flow" -M 100 -v -o /data/output.bam -b /data/primers.bed -t /data/pairs.tsv
   ```

By combining these flags, you can tailor the behavior and output of the application to meet your specific needs. Ensure that the required flags `-i`, `-a`, and `-M` are always included in your command to avoid errors.

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
3. navigate to the repository directory and run `cmake .`.



