# cuda-dna-downsampling

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

### Branch naming convention
`(feature|bug)/{issue_number}(-{feature_name})?`

ex.
feature/4-add-min-cost-max-flow-library
bug/7

