# Ingest pipeline (adapted) to create minimal Nextcalde datasets

## Requirements

### required files

- subdirectory data/pathogen_name/ must have: an annotation.gff3 file, pathogen.json, reference.fasta, and reference.gbk files

data folder with virus name, inside contains the annotation.gff3, reference.gbk and reference.fasta created previously with script bin/generate from genbank.py. Additionally, a pathogen.json file must be present in each folder.

should also add coordinates of special protein of interest in extend-metadata.py

the config file must contain for each virus: the taxon id, minimum length to filter on, and the reference (called root) under refine.

## outputs

in dataset folder:

- the Nextclade dataset for each configured pathogen

## running

- from a terminal, launch Snakemake, using the environment provided in envs/ingest-genetic-dataset.yaml
