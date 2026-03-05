# pipeline to clean the metadata associated with each virus

VIRUSES = list(config.get("viruses", []))


rule clean_metadata:
    input:
        metadata="data/local_data/{strain}/metadata.tsv",
    output:
        metadata_sequencing="tmp/{strain}/metadata_sequencing.tsv",
        metadata_pcr="tmp/{strain}/metadata_pcr.tsv",
    log:
        "logs/clean_metadata/{strain}.log",
    script:
        "scripts/clean_metadata.py"


rule assemble_metadata:
    input:
        metadata_sequencing=expand(
            "tmp/{strain}/metadata_sequencing.tsv", strain=VIRUSES
        ),
        metadata_pcr=expand("tmp/{strain}/metadata_pcr.tsv", strain=VIRUSES),
    output:
        metadata_sequencing="results/dashboard-metadata/metadata_sequencing.tsv",
        metadata_pcr="results/dashboard-metadata/metadata_pcr.tsv",
    run:
        # concatenate metadata files for all strains into a single file for the dashboard
        import pandas as pd

        all_metadata_sequencing = pd.concat(
            [pd.read_csv(f, sep="\t") for f in input.metadata_sequencing],
            ignore_index=True,
        )
        all_metadata_sequencing.to_csv(
            output.metadata_sequencing, sep="\t", index=False
        )

        all_metadata_pcr = pd.concat(
            [pd.read_csv(f, sep="\t") for f in input.metadata_pcr],
            ignore_index=True,
        )
        all_metadata_pcr.to_csv(output.metadata_pcr, sep="\t", index=False)

        # delete folder tmp
        import shutil

        shutil.rmtree("tmp")
