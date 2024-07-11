# MutChecker

MutChecker is a tool for checking the mutation status of a given genome range. 


## Installation

```bash
pip install mutchecker
```

## Usage

```bash
mutchecker -r accession.bam -r chr1:1000-2000 -o output_prefix
```

```bash
mutchecker -p accession.genome.fasta -q query.fasta -o output_prefix
```