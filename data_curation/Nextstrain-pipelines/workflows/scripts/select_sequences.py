"""Keep the consensus genomes of the samples listed in a virus's metadata.

A sequence belongs to a sample when the part of its name before the first "|" equals the
sample's SampleID (so "S1|2025-01-20|BE" matches "S1"). Sequences of other samples are left
out and listed in the log. A sample with several sequences stops the pipeline, since it is
unclear which one belongs in the tree.
"""

from collections import Counter

import pandas as pd


def sample_id(name: str) -> str:
    return str(name).split("|")[0].strip()


def read_fasta(path: str) -> list[tuple[str, str]]:
    records = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                records.append((line[1:].strip(), []))
            elif line.strip():
                records[-1][1].append(line.strip())
    return [(name, "".join(seq)) for name, seq in records]


ids = set(pd.read_csv(snakemake.input.metadata, sep="\t", dtype=str)["SampleID"].map(sample_id))
records = read_fasta(snakemake.input.sequences)
kept = [(name, seq) for name, seq in records if sample_id(name) in ids]
found = Counter(sample_id(name) for name, _ in kept)

with open(snakemake.log[0], "w") as log:
    log.write(f"kept {len(kept)} of {len(records)} sequences\n")
    log.writelines(f"left out, not in metadata: {name}\n" for name, _ in records if sample_id(name) not in ids)
    log.writelines(f"no sequence for: {i}\n" for i in sorted(ids - set(found)))

repeated = sorted(i for i, n in found.items() if n > 1)
if repeated:
    raise ValueError(
        f"{snakemake.input.sequences}: several sequences for {', '.join(repeated)}; keep one per sample."
    )

with open(snakemake.output.sequences, "w") as out:
    out.writelines(f">{name}\n{seq}\n" for name, seq in kept)
