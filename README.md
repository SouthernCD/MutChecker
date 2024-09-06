# MutChecker

MutChecker is a tool for checking the mutation status of a given genome range. 


## Installation

```bash
mamba create -n mutchecker -c bioconda -y python==3.9 biopython==1.80 freebayes whatshap bcftools samtools htslib numexpr miniprot blast
conda activate mutchecker
pip install mutchecker
```

## Usage

1. `bamstat`: Stat the depth and coverage of all CDS in a bam file

```bash
usage: mutchecker bamstat [-h] gff_file bam_file

Stat the depth and coverage of all CDS in a bam file

positional arguments:
  gff_file    reference gff file
  bam_file    sorted and markdup bam file

optional arguments:
  -h, --help  show this help message and exit
```

2. `reseq`: Check the mutation of a specific gene by resequencing data

```bash
usage: mutchecker reseq [-h] [-i GENE_ID] [-o OUTPUT_DIR] [-e EXON_EXTEND] [-c] genome_file gff_file bam_file

Check the mutation of a specific gene by resequencing data

positional arguments:
  genome_file           reference genome file
  gff_file              reference gff file
  bam_file              sorted and markdup bam file

optional arguments:
  -h, --help            show this help message and exit
  -i GENE_ID, --gene_id GENE_ID
                        gene id, default=None
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        output directory, default="mutchecker_output"
  -e EXON_EXTEND, --exon_extend EXON_EXTEND
                        extend length of cds, default=500
  -c, --clean           clean the intermediate files, default False
```