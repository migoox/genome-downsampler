# cuda-dna-downsampling

## Command Line Interface (CLI) Flags

This section details the available CLI flags for configuring the application. Each flag allows you to customize the behavior and output of the program.

### Required Flags
The following flags must be specified for the program to run:

- `-i`, `--input` : **Input file path**
  - This flag is required to specify the path to the input file.
  - Example: `-i /path/to/input/file`

- `-a`, `--algorithm` : **Algorithm ID (descriptive string)**
  - This flag is required to specify which algorithm to use. Currently available algorithms are:
    - test - test algorithm just to make sure CLI works properly
    - sequential-cost-scaling - deprecated algorithm of MinCostMaxFlow with cost scaling
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

- `-c`, `--csv-history` : **Path to CSV file**
  - Specifies the path to a CSV file where the results will be appended after algorithm execution. If not specified, the default is `./historical_runs.csv` in the current directory.
  - Example: `-c /path/to/history.csv`

- `-o`, `--output` : **Output file path**
  - Sets the path for the output file. If not specified, the default is `output.bam` in the input file's directory.
  - Example: `-o /path/to/output.bam`

### Usage Examples

1. **Basic usage with required flags:**
   ```sh
   program -i /data/input.bam -a "test" -M 100
   ```

2. **Using all optional flags:**
   ```sh
   program -i /data/input.bam -a "test" -M 100 -v -l 100 -q 50 -c /data/history.csv -o /data/output.bam
   ```

3. **Verbose mode with custom CSV history and output path:**
   ```sh
   program -i /data/input.bam -a "test" -M 100 -v -c /data/history.csv -o /data/output.bam
   ```

By combining these flags, you can tailor the behavior and output of the application to meet your specific needs. Ensure that the required flags `-i`, `-a`, and `-M` are always included in your command to avoid errors.

### Branch naming convention
`(feature|bug)/{issue_number}(-{feature_name})?`

ex.
feature/4-add-min-cost-max-flow-library
bug/7

