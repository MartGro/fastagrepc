# FastaGrepC Manual

A fast C implementation for searching patterns in gzipped FASTA files using the Aho-Corasick algorithm.

## Features
- Fast multi-pattern searching using Aho-Corasick algorithm
- Searches both forward and reverse strands
- Case-insensitive search option
- Supports gzipped FASTA files
- Configurable sequence context around matches
- CSV output format

## Installation

```bash
gcc -o fastagrepc main.c -lz
```

## Usage

```bash
./fastagrepc <fasta_file> <patterns_file> [context] [sequence_only] [ignore_case]
```

### Parameters

1. `fasta_file`: Input FASTA file (can be gzipped)
2. `patterns_file`: CSV file containing patterns to search for
3. `context`: Number of bases to include before and after match (default: 0)
4. `sequence_only`: Search only in sequences, not headers (1=yes, 0=no, default: 0)
5. `ignore_case`: Case-insensitive search (1=yes, 0=no, default: 0)

### Pattern File Format

The pattern file should be a CSV with header row and two columns:
```csv
name,sequence
Pattern1,ATCG
Pattern2,GCTA
```

### Output Format

CSV format with the following columns:
- header: FASTA sequence header
- pattern_name: Name of the matched pattern
- pattern_sequence: Sequence of the matched pattern
- position: Position in the sequence (0-based)
- strand: "forward" or "reverse"
- context: Sequence context including the match

## Examples

### Basic Search
Search for patterns in a FASTA file:

```bash
./fastagrepc input.fa.gz patterns.csv > results.csv
```

### With Context
Include 10 bases before and after each match:

```bash
./fastagrepc input.fa.gz patterns.csv 10 > results.csv
```

### Case-Insensitive Search
Search patterns case-insensitively with 5 bases context:

```bash
./fastagrepc input.fa.gz patterns.csv 5 0 1 > results.csv
```

### Sample patterns.csv
```csv
name,sequence
PAM1,NGG
Pattern1,ATCGATCG
sgRNA1,GUUUUAGAGCUA
```

### Example Output
```csv
header,pattern_name,pattern_sequence,position,strand,context
chr1,PAM1,NGG,145,forward,ACCNGGTAT
chr1,sgRNA1,GUUUUAGAGCUA,2010,reverse,TAGCTCTAAAAC
```

## Performance Considerations
- Uses Aho-Corasick for efficient multi-pattern matching
- Memory usage scales with:
  - Input sequence length
  - Number of patterns
  - Context size
- 1MB buffer size for reading compressed files

## Notes
- Patterns containing commas will have them replaced with semicolons in output
- Reverse strand positions are reported in forward strand coordinates
- N in patterns matches any base
