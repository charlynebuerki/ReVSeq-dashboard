"""
This part of the workflow handles sorting downloaded sequences and metadata
by aligning them to reference sequences.

It produces output files as

    metadata = "data/{strain}/metadata.tsv"
    sequences = "data/{strain}/sequences.fasta"

"""


rule sort:
    input:
        sequences=rules.curate.output.sequences,
    output:
        sequences="results/{strain}/sequences_sorted.fasta",
    shell:
        """
        seqkit rmdup {input.sequences} > {output}
        """


rule metadata:
    input:
        metadata=rules.subset_metadata.output.subset_metadata,
        sequences="results/{strain}/sequences_sorted.fasta",
    output:
        metadata="data/{strain}/metadata_raw.tsv",
    run:
        import pandas as pd
        from Bio import SeqIO

        strains = [s.id for s in SeqIO.parse(input.sequences, "fasta")]
        d = (
            pd.read_csv(input.metadata, sep="\t", index_col="accession")
            .loc[strains]
            .drop_duplicates()
        )
        d.to_csv(output.metadata, sep="\t")


rule nextclade:
    input:
        sequences="results/{strain}/sequences_sorted.fasta",
        ref="data/{strain}/reference.fasta",
    output:
        nextclade="results/{strain}/nextclade.tsv",
    conda:
        "../envs/ingest-genetic-datasets.yaml"
    params:
        dataset="data/{strain}/",
        #output_columns="seqName clade qc.overallScore qc.overallStatus alignmentScore  alignmentStart  alignmentEnd  coverage dynamic",
    threads: 8
    shell:
        """
        nextclade run -D {params.dataset}  -j {threads} \
                          --output-tsv {output.nextclade} \
                          {input.sequences}
        """


rule extend_metadata:
    input:
        nextclade="results/{strain}/nextclade.tsv",
        metadata="data/{strain}/metadata_raw.tsv",
    output:
        metadata="results/{strain}/metadata.tsv",
    conda:
        "../envs/ingest-genetic-datasets.yaml"
    shell:
        """
        python3 bin/extend-metadata.py --metadata {input.metadata} \
                                       --id-field accession \
                                       --virus-type {wildcards.strain} \
                                       --nextclade {input.nextclade} \
                                       --output {output.metadata}
        """


# keep the root of the tree (reference sequence)
def get_root(wildcards):
    import os

    root = config["refine"]["root"][wildcards.strain]
    root_file = f"results/{wildcards.strain}/root.txt"

    # Create directory if it doesn't exist
    os.makedirs(f"results/{wildcards.strain}", exist_ok=True)

    # Create the file if it doesn't exist
    if not os.path.exists(root_file):
        with open(root_file, "w") as f:
            f.write(root)

    return root_file


# filter to only high-quality sequences & representative sequences for each country
# don't include swiss sequences from our project
rule filter:
    message:
        """
        Filtering recent sequences to: 
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - from {params.min_date} to {params.max_date}
        """
    input:
        sequences="results/{strain}/sequences_sorted.fasta",
        metadata="results/{strain}/metadata.tsv",
    output:
        sequences="results/{strain}/filtered_initial.fasta",
        metadata="results/{strain}/filtered_metadata_initial.tsv",
        log="logs/{strain}/filtered.log",
    conda:
        "../envs/ingest-genetic-datasets.yaml"
    log:
        "logs/filtering/{strain}.log",
    params:
        group_by=config["filter"]["group_by"],
        root=lambda wildcards: get_root(wildcards),
        sequences_per_group=config["filter"]["subsample_max_sequences"]["recent"],
        strain_id_field="accession",
        min_date=lambda wildcards: config["filter"]["resolution"]["min_date"],
        max_date=lambda wildcards: config["filter"]["resolution"]["max_date"],
        min_length=lambda wildcards: config["filter"]["min_length"][wildcards.strain],
        min_coverage=f"genome_coverage>{config['filter']['min_coverage']}",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --include {params.root} \
            --max-date {params.max_date} \
            --min-length {params.min_length} \
            --query '{params.min_coverage} & bioproject_accession != "PRJEB83635" & (`qc.overallStatus` == "good" | `qc.overallStatus` == "mediocre") ' \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --output-log {output.log} \
             > {log} 2>&1
        """
