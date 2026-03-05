# pipeline to download background nextstrain dataset


def _get_virus_types(wildcard):
    dataset = config["nextclade_datasets"].get(wildcard, "")
    # Treat only Nextstrain dataset names as nextclade; paths/URLs are custom zips
    if isinstance(dataset, str) and dataset.startswith("nextstrain/"):
        return "nextclade"
    return "other"


rule fetch_nextclade_datasets:
    conda:
        "../envs/python-genetic-data.yaml"
    # logging handled inside shell to avoid wildcard mismatch issues
    # log:
    #     "logs/download_nextclade_datasets/{virus}.log",
    # Produce only a sentinel file; downstream rules will depend on this
    # sentinel and reference the dataset directory path directly. This
    # avoids `dir(...)` outputs which can trigger wildcard-matching issues
    # across outputs/logs for some Snakemake versions.
    output:
        sentinel="data/nextclade/{virus}/.dataset_ready",
    params:
        dataset_name=lambda wildcards: config["nextclade_datasets"][wildcards.virus],
        virus_type=lambda wildcards: _get_virus_types(wildcards.virus),
    shell:
        """
        mkdir -p data/nextclade/{wildcards.virus}
        if [ "{params.virus_type}" = "nextclade" ]; then
            nextclade dataset get \
                --name {params.dataset_name} \
                --output-dir data/nextclade/{wildcards.virus} 
        else
            # for other datasets, we assume the dataset_name is a URL to a custom dataset zip file
            if [ -f "{params.dataset_name}" ]; then
                cp "{params.dataset_name}" data/nextclade/{wildcards.virus}/dataset.zip
            elif [[ "{params.dataset_name}" == file://* ]]; then
                cp "${{params.dataset_name#file://}}" data/nextclade/{wildcards.virus}/dataset.zip
            else
                curl -L "{params.dataset_name}" -o data/nextclade/{wildcards.virus}/dataset.zip
            fi
            unzip -o data/nextclade/{wildcards.virus}/dataset.zip -d data/nextclade/{wildcards.virus}
        fi
        touch {output.sentinel}
        """


rule calculate_tree:
    conda:
        "../envs/python-genetic-data.yaml"
    log:
        "logs/calculate_nextclade_tree/{virus}.log",
    input:
        dataset_ready="data/nextclade/{virus}/.dataset_ready",
        local_sequences="data/local_data/{virus}/sequences.fasta",
    output:
        fasta="results/nextclade/{virus}/sequences.fasta",
        metadata="results/nextclade/{virus}/metadata.tsv",
        tree="results/nextclade/{virus}/tree.json",
    shell:
        """
        nextclade run \
            --input-dataset data/nextclade/{wildcards.virus} \
            --output-tree {output.tree} \
            --output-fasta {output.fasta} \
            --output-tsv {output.metadata} \
            {input.local_sequences} \
            &> {log}
        """


# to filter to only new sequences add this to url: f_Node%20type=New
rule export_auspice_trees:
    conda:
        "../envs/python-genetic-data.yaml"
    log:
        "logs/export_auspice_trees/{virus}.log",
    input:
        tree="results/nextclade/{virus}/tree.json",
    output:
        sentinel_ready="results/auspice_trees/.export_ready_{virus}",
    params:
        nice_virus_name=lambda wildcards: config["nice_virus_names"].get(
            wildcards.virus
        ),
    shell:
        """
        python -c "
import json
with open('{input.tree}', 'r') as f:
    data = json.load(f)
if 'meta' not in data:
    data['meta'] = {{}}
if 'display_defaults' not in data['meta']:
    data['meta']['display_defaults'] = {{}}
data['meta']['display_defaults']['sidebar'] = 'closed'
with open('results/auspice_trees/{params.nice_virus_name}.json', 'w') as f:
    json.dump(data, f, indent=2)
"
        touch {output.sentinel_ready}
        """
