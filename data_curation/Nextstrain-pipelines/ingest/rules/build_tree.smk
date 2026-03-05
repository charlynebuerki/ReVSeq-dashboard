# workflow to create a nextclade tree


# get the genes to translate for virus
def _get_genes_for_virus(wildcard):
    # read the gff annotaiton file for the virus and select the gene, return a list of genes,
    # gff file is expect to have a "Name=X" on each line for the gene name, where X is gene name
    import pandas as pd

    gff_file = f"results/nextclade/{wildcard}/annotation.gff3"
    gff = pd.read_csv(
        gff_file,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ],
    )
    # Prefer CDS features since this GFF uses CDS entries with Name=...
    gene_source = gff[gff["type"] == "CDS"]
    if gene_source.empty:
        gene_source = gff[gff["type"] == "gene"]
    genes = (
        gene_source["attributes"]
        .str.extract(r"Name=([^;]+)")[0]
        .dropna()
        .unique()
        .tolist()
    )
    print(f"Genes found in {gff_file}: {genes}")
    return genes


def _get_genes_string(wildcards):
    genes = _get_genes_for_virus(wildcards.strain)
    print(f"Genes for {wildcards.strain}: {genes}")
    return " ".join(genes)


# copy everything to nextclade directory in results
rule copy_to_nextclade:
    input:
        sequences="results/{strain}/filtered_initial.fasta",
        metadata="results/{strain}/filtered_metadata_initial.tsv",
        pathogen_json="defaults/pathogen.json",
        genome_annotation="data/{strain}/annotation.gff3",
        genbank_file="data/{strain}/reference.gbk",
        reference="data/{strain}/reference.fasta",
    output:
        sequences="results/nextclade/{strain}/filtered_pre.fasta",
        metadata="results/nextclade/{strain}/metadata_pre.tsv",
        pathogen_json="results/nextclade/{strain}/pathogen.json",
        genome_annotation="results/nextclade/{strain}/annotation.gff3",
        genbank_file="results/nextclade/{strain}/reference.gbk",
        reference="results/nextclade/{strain}/reference.fasta",
        confirm_ready="results/nextclade/{strain}/.dataset_ready",
    shell:
        """
        mkdir -p results/nextclade/{wildcards.strain} 
        cp {input.sequences} {output.sequences}
        cp {input.metadata} {output.metadata}
        cp {input.pathogen_json} {output.pathogen_json}
        cp {input.genome_annotation} {output.genome_annotation}
        cp {input.reference} {output.reference}
        cp {input.genbank_file} {output.genbank_file}
        touch {output.confirm_ready}
        """


rule align:
    input:
        sequences="results/nextclade/{strain}/filtered_pre.fasta",
        reference="results/nextclade/{strain}/reference.fasta",
        annotation="results/nextclade/{strain}/annotation.gff3",
    output:
        alignment="results/nextclade/{strain}/aligned.fasta",
        tsv="results/nextclade/{strain}/nextclade.tsv",
    params:
        translation_template=lambda w: f"results/nextclade/{w.strain}/translations/cds_{{cds}}.translation.fasta",
    conda:
        "../envs/ingest-genetic-datasets.yaml"
    log:
        "logs/align/{strain}.log",
    shell:
        """
        nextclade run \
            {input.sequences} \
            --input-ref {input.reference} \
            --input-annotation {input.annotation} \
            --output-translations {params.translation_template} \
            --output-tsv {output.tsv} \
            --output-fasta {output.alignment} \
            &> {log}
        """


rule get_outliers:
    """
    Automatically identify sequences with >{params.allowed_divergence} substitutions
    (likely to be sequencing errors or low quality/misannotated sequences) and put them in outliers.txt
    """
    input:
        nextclade="results/nextclade/{strain}/nextclade.tsv",
    output:
        outliers="results/nextclade/{strain}/outliers.txt",
        tmp="tmp/{strain}/outliers.txt",
    conda:
        "../envs/ingest-genetic-datasets.yaml"
    params:
        allowed_divergence=config["filter"]["allowed_divergence"],
    shell:
        """
        tsv-filter -H -v --is-numeric totalSubstitutions {input.nextclade} \
        > {output.tmp}
        tsv-filter -H \
            --is-numeric totalSubstitutions \
            --gt totalSubstitutions:{params.allowed_divergence} \
            {input.nextclade} \
        | tail -n +2 >> {output.tmp}
        cat {output.tmp} \
        | tsv-select -H -f seqName \
        | tail -n +2 > {output.outliers} 
        """


rule exclude:
    """
    Rule to allow for manual and automatic exclusion of sequences
    without triggering a new subsampling that could
    surface new bad sequences resulting in an infinite loop
    """
    input:
        sequences="results/nextclade/{strain}/aligned.fasta",
        metadata="results/nextclade/{strain}/metadata_pre.tsv",
        outliers="results/nextclade/{strain}/outliers.txt",
    output:
        filtered_sequences="results/nextclade/{strain}/filtered.fasta",
        filtered_metadata="results/nextclade/{strain}/filtered_metadata.tsv",
        strains="results/nextclade/{strain}/tree_strains.txt",
    conda:
        "../envs/ingest-genetic-datasets.yaml"
    log:
        "logs/exclude/{strain}.log",
    params:
        strain_id_field="accession",
        root_strain=lambda wildcards: get_root(wildcards),
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata-id-columns {params.strain_id_field} \
            --metadata {input.metadata} \
            --exclude  {input.outliers} \
            --include {params.root_strain} \
            --output-sequences {output.filtered_sequences} \
            --output-metadata {output.filtered_metadata} \
            --output-strains {output.strains} \
            > {log} 2>&1
        """


rule tree:
    input:
        alignment="results/nextclade/{strain}/filtered.fasta",
    output:
        tree="results/nextclade/{strain}/tree_raw.nwk",
    conda:
        "../envs/ingest-genetic-datasets.yaml"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
        """


# root using dates in treetime, use 15000 as sequence length (good enough, doesn't matter)
rule root:
    input:
        tree=rules.tree.output.tree,
        metadata=rules.exclude.output.filtered_metadata,
    output:
        tree="results/nextclade/{strain}/tree_rooted.nwk",
    params:
        outdir="results/nextclade/{strain}/tt_out",
    conda:
        "../envs/ingest-genetic-datasets.yaml"
    shell:
        """
        treetime clock \
            --tree {input.tree} \
            --sequence-length 15000 \
            --dates {input.metadata} \
            --name-column accession \
            --clock-filter 5 \
            --clock-filter-method residual \
            --outdir {params.outdir}
        cp {params.outdir}/rerooted.newick {output.tree}
        """


rule prune_outliers:
    input:
        tree=rules.root.output.tree,
    output:
        tree="results/nextclade/{strain}/tree_rooted_pruned.nwk",
    params:
        outliers="results/nextclade/{strain}/tt_out/outliers.tsv",
    run:
        import pandas as pd
        from Bio import Phylo
        import os
        import shutil

        if not os.path.exists(params.outliers):
            # If outliers file doesn't exist, just copy the tree as-is
            shutil.copy(input.tree, output.tree)
        else:
            # Prune outliers if the file exists
            outliers = pd.read_csv(params.outliers, sep="\t", index_col=0)
            T = Phylo.read(input.tree, "newick")

            for n in outliers.index:
                if outliers.loc[n, "given_date"] > 1980 and ("-egg" not in n):
                    print("prune", n)
                    T.prune(n)
            Phylo.write(T, output.tree, "newick")


rule refine:
    input:
        tree="results/nextclade/{strain}/tree_rooted_pruned.nwk",
        alignment="results/nextclade/{strain}/filtered.fasta",
        metadata="results/nextclade/{strain}/filtered_metadata.tsv",
    output:
        tree="results/nextclade/{strain}/tree_refined.nwk",
        node_data="results/nextclade/{strain}/branch_lengths.json",
    params:
        #root=lambda wildcards: config["refine"]["root"][wildcards.strain],
        root="mid_point",
        strain_id_field="accession",
    log:
        "logs/refine/{strain}.log",
    conda:
        "../envs/ingest-genetic-datasets.yaml"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --root {params.root} \
            --metadata-id-columns {params.strain_id_field} \
            --metadata {input.metadata} \
            --keep-polytomies \
            --keep-root \
            --divergence-units mutations-per-site \
            --output-node-data {output.node_data} \
            --output-tree {output.tree} \
             > {log} 2>&1
        """


rule ancestral:
    input:
        tree="results/nextclade/{strain}/tree_refined.nwk",
        alignment="results/nextclade/{strain}/filtered.fasta",
        reference="results/nextclade/{strain}/reference.fasta",
        annotation="results/nextclade/{strain}/reference.gbk",
        gff_file="results/nextclade/{strain}/annotation.gff3",
    output:
        node_data="results/nextclade/{strain}/muts.json",
    params:
        translation_template=r"results/nextclade/{strain}/translations/cds_%GENE.translation.fasta",
        output_translation_template=r"results/nextclade/{strain}/translations/cds_%GENE.ancestral.fasta",
        genes=_get_genes_string,
        inference="joint",
    conda:
        "../envs/ingest-genetic-datasets.yaml"
    log:
        "logs/ancestral/{strain}.log",
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --inference {params.inference} \
            --infer-ambiguous \
            --annotation {input.annotation} \
            --genes {params.genes} \
            --translations {params.translation_template} \
            --output-node-data {output.node_data} \
            --output-translations {params.output_translation_template} \
            --root-sequence {input.reference} \
             > {log} 2>&1 
        """


rule dummy_clades:
    """
    Nextclade requires clade membership to be specified for each node
    in the tree. This rule creates a dummy clade membership for each node
    """
    input:
        "results/nextclade/{strain}/branch_lengths.json",
    output:
        "results/nextclade/{strain}/dummy_clades.json",
    shell:
        """
        jq '.nodes |= map_values({{"clade_membership": "dummy"}})' {input} > {output}
        """


rule export:
    input:
        tree="results/nextclade/{strain}/tree_refined.nwk",
        metadata="results/nextclade/{strain}/filtered_metadata.tsv",
        mutations="results/nextclade/{strain}/muts.json",
        branch_lengths="results/nextclade/{strain}/branch_lengths.json",
        clades="results/nextclade/{strain}/dummy_clades.json",
        auspice_config="defaults/auspice_config.json",
    output:
        auspice="results/nextclade/{strain}/auspice.json",
    params:
        strain_id_field="accession",
    log:
        "logs/export/{strain}.log",
    conda:
        "../envs/ingest-genetic-datasets.yaml"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --auspice-config {input.auspice_config} \
            --node-data {input.mutations} {input.branch_lengths} {input.clades} \
            --output {output.auspice} \
             > {log} 2>&1
        """


rule rename_tips:
    input:
        auspice_json="results/nextclade/{strain}/auspice.json",
        metadata="results/nextclade/{strain}/filtered_metadata.tsv",
    output:
        auspice_json="results/nextclade/{strain}/auspice_renamed.json",
    run:
        import pandas as pd

        name_map = {
            d.accession: d.strain
            for d in pd.read_csv(input.metadata, sep="\t").itertuples()
        }

        import json

        with open(input.auspice_json) as fh:
            data = json.load(fh)


        def rename(n, name_map):
            n["name"] = name_map.get(n["name"], n["name"])
            if "children" in n:
                for c in n["children"]:
                    rename(c, name_map)


        rename(data["tree"], name_map)

        with open(output.auspice_json, "w") as fh:
            json.dump(data, fh)


rule subsample_example_sequences:
    input:
        all_sequences="results/nextclade/{strain}/filtered_pre.fasta",
        tree_strains="results/nextclade/{strain}/tree_strains.txt",
    output:
        example_sequences="results/nextclade/{strain}/example_sequences.fasta",
    conda:
        "../envs/ingest-genetic-datasets.yaml"
    shell:
        """
        # Exclude tree sequences from all sequences
        seqkit grep -v -f {input.tree_strains} {input.all_sequences} \
        | seqkit sample -n 30 -s 42 > {output.example_sequences}
        """


rule assemble_dataset:
    input:
        tree="results/nextclade/{strain}/auspice_renamed.json",
        reference="results/nextclade/{strain}/reference.fasta",
        annotation="results/nextclade/{strain}/annotation.gff3",
        sequences="results/nextclade/{strain}/example_sequences.fasta",
        pathogen="results/nextclade/{strain}/pathogen.json",
    output:
        tree="dataset/{strain}/tree.json",
        reference="dataset/{strain}/reference.fasta",
        annotation="dataset/{strain}/annotation.gff3",
        sequences="dataset/{strain}/sequences.fasta",
        readme="dataset/{strain}/README.md",
        changelog="dataset/{strain}/CHANGELOG.md",
        pathogen="dataset/{strain}/pathogen.json",
        dataset_zip="dataset/dataset_{strain}.zip",
    log:
        "logs/assemble_dataset/{strain}.log",
    shell:
        """
        cp {input.tree} {output.tree}
        cp {input.reference} {output.reference}
        cp {input.annotation} {output.annotation}
        cp {input.sequences} {output.sequences}
        cp {input.pathogen} {output.pathogen}
        touch {output.readme}
        touch {output.changelog}
        zip -rj {output.dataset_zip}  dataset/{wildcards.strain}/* > {log} 2>&1

        """


rule test:
    input:
        dataset="dataset/dataset_{strain}.zip",
        sequences="dataset/{strain}/sequences.fasta",
    output:
        output=directory("test_out/{strain}"),
        output_confirm="test_out/{strain}/.test_ready",
    conda:
        "../envs/ingest-genetic-datasets.yaml"
    log:
        "logs/test/{strain}.log",
    shell:
        """
        nextclade run \
            --input-dataset {input.dataset} \
            --output-all {output.output} \
            {input.sequences} \
            &> {log}
        touch {output.output_confirm}
        """
