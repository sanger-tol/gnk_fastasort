# gnk_fastasort
A package to create a preferred order of a FASTA scaffolds, useful for PretextSnapshot and SAMTOOLS reordering.

This will reorganise a fasta file into groupings, based on chromosome naming. So unlocs of SUPER_1 will be re-ordered alongside it in size order. Other unlocs will be arranged in size order, after everything else.

## Installation
```
git clone https://github.com/sanger-tol/gnk_fastasort.git

cd gnk_fastasort

pip install ./

fastasort -h
```

## Usage 
This script can be used in two ways, both however rely on the underlying naming of the fasta being SUPER_* style. In neither case is a fasta file required, it is designed to work in conjunction with SAMTOOLS faidx which will can re-organise a fasta given a tsv file of names. This tool emits a tsv of `names`, `original naming`, `parent molecule` and `length` or just `names`.

For example, a fasta fai file where scaffolds are named SUPER_* or a GCA accessioned fasta which was originally named SUPER_* (this is visible from the sequence report on ncbi.).

For genomes stored on NCBI and are accessioned. You can use the `--gca_accession` to call the API and download a sequence report.
```
fastasort --gca_accession GCA_964017025.1
```

For local fasta names where headers are `SUPER_1`, `SUPER_2`, `SUPER_2_unloc_1` e.g. SUPER_* style, you can run:
```
fastasort --index {ASSEMBLY}.fa.fai
```
