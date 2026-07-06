# Frame Resolved Assembly for Metagenomics

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Rust 1.70+](https://img.shields.io/badge/rust-1.70+-orange)

## Quick Start

### Installation

Prerequisites:
- Rust 1.70+ ([install](https://www.rust-lang.org/tools/install))
- FragGeneScanRs module (included in repository)

Build FRAME:

```bash
cargo build --release
```

### Basic Usage

```bash
# Run with default parameters (includes read rescue)
./target/release/frame input.fastq

# Specify output directory and k-mer parameters
./target/release/frame input.fastq \
  --output results/ \
  --kmer 31 \
  --min-count 2 \
  --min-length 100

# Use custom HMM model
./target/release/frame input.fastq \
  --hmm-dir ./lib/FragGeneScanRs/train \
  --hmm-model illumina_5

```

## Command-Line Interface

```
Usage: frame [OPTIONS] <INPUT>

ARGUMENTS:
  <INPUT>                Path to input sequencing reads

OPTIONS:
  -o, --output <DIR>              Output directory [default: ./frame_output]
  -k, --kmer <SIZE>               K-mer size (15-63) [default: 31]
  -m, --min-count <COUNT>         Minimum k-mer count threshold [default: 2]
  -l, --min-length <LENGTH>       Minimum unitig length in bp [default: 100]
      --hmm-dir <DIR>             Path to HMM training directory [default: ./lib/FragGeneScanRs/train]
      --hmm-model <NAME>          HMM model name [default: illumina_5]
      --log-level <LEVEL>         Log level: debug, info, warn, error [default: info]
  -h, --help                      Print help
  -V, --version                   Print version
```


### Output Files
- `predictions.gff`: Gene coordinates and annotations (assembly + rescue)
- `proteins.faa`: Translated protein sequences (assembly + rescue)
- `genes.fna`: DNA sequences of predicted genes (assembly + rescue)

### Building from Source

```bash
# Development build with optimizations
cargo build

# Release build with full optimizations
cargo build --release

# Run with specific log level
RUST_LOG=debug ./target/release/frame input.fastq
```
