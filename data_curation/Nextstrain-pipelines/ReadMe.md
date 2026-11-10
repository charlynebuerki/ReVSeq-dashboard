# Pre-processing data for dashboard integration

This pipeline formats the sequences and metadata required for dashboard integration.

## inputs

local data under folder: data/local_data/strain_name containing 2 files:

- metadata.tsv with minimal column names SampleID, date, virus_identified, location, Match_PCR if PRC results
- consensus genomes contained in a multi-fasta format of the same virus.

Only sequences of samples listed in metadata.tsv go into the tree: a sequence belongs to a sample when its name, up to the first `|`, equals the SampleID. Each sample must have exactly one sequence; the pipeline stops otherwise. Left-out sequences are listed in logs/select_sequences/strain_name.log.

You can parametrize which strains to compute in the configfile.yaml. This supports both existing Nextclade datasets and custom datasets. Custom datasets must be configured in the ingest/ folder. See the ReadMe file for that directory for more details.

## running

You can launch the workflow using Snakemake, using the envrionment provided in envs/.

## output

requested strain datasets will be found under results/auspice_trees/virus_name.json . The dashboard metadata formatted will be under results/dashboard-metadata/, including metadata_co_infection_sequencing.tsv: the samples with more than one virus, for the co-infection page.
