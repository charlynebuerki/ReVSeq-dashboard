# Data curation overview

This repository formats the data inputs for the dashboard in 2 pipelines:

- Nextstrain-pipeline: for phylogeny generation and metadata formatting
- Extract-pileup: for pileup extraction

See the ReadMe files in each sub-directory for implementation.

## Requirements

- Nextstrain-pipeline:
  - consensus genomes of selected virus
  - associated metadata
  - configurations detailed in the ReadMe file

- Extract-pileup:
  - .bam file for each sample
  - configurations detailed in the ReadMe file

## Using the outputs

Once the pipelines have been generated, copy over the formatted data to the data/ directory of the dashboard
