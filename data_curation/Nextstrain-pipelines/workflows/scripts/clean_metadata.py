import pandas as pd


def _normalize_virus_text(value: object) -> set[str]:
    if pd.isna(value):
        return set()
    text = str(value).lower()
    # drop parenthetical subtype/strain details
    text = text.replace("(", " ").replace(")", " ")
    # remove common qualifiers
    for token in ["human ", "virus ", "viruses "]:
        text = text.replace(token, "")
    # normalize separators
    for sep in [";", "/", "|", "+"]:
        text = text.replace(sep, ",")

    # should remove anything in between parantheses, e.g. "coronavirus 229e (hCoV-229E)" -> "coronavirus 229e"
    import re

    text = re.sub(r"\(.*?\)", "", text)
    # normalize whitespace
    text = " ".join(text.split())
    # split into items
    parts = [p.strip() for p in text.split(",") if p.strip()]
    # normalize common synonyms
    synonyms = {
        "influenza a": "influenza a",
        "influenza a h3n2": "influenza a",
        "influenza a h1n1": "influenza a",
        "influenza a h5n1": "influenza a",
        "influenza a h7n9": "influenza a",
        "influenza b": "influenza b",
        "coronavirus 229e": "coronavirus 229e",
        "coronavirus hk u1": "coronavirus hku1",
        "coronavirus oc43": "coronavirus oc43",
        "coronavirus nl63": "coronavirus nl63",
        "human coronavirus 229e": "coronavirus 229e",
        "human coronavirus hku1": "coronavirus hku1",
        "human coronavirus oc43": "coronavirus oc43",
        "human coronavirus nl63": "coronavirus nl63",
        "parainfluenza 1": "parainfluenza 1",
        "parainfluenza 2": "parainfluenza 2",
        "parainfluenza 3": "parainfluenza 3",
        "parainfluenza 4": "parainfluenza 4",
        "rsv": "rsv",
        "respiratory syncytial": "rsv",
        "respiratory syncytial virus": "rsv",
        "Influenza A": "influenza a",
        "Influenza B": "influenza b",
    }

    normalized = set()
    if "influenza a" in text:
        normalized.add("influenza a")
    if "influenza b" in text:
        normalized.add("influenza b")
    for p in parts:
        p = p.replace("-", " ")
        p = " ".join(p.split())
        if p.startswith("influenza a "):
            p = "influenza a"
        if p.startswith("influenza b "):
            p = "influenza b"
        # Handle RSV variants: "rsv a", "rsv b", "rsv type a", etc.
        if p.startswith("rsv ") or p.startswith("respiratory syncytial "):
            p = "rsv"
        normalized.add(synonyms.get(p, p))

    return normalized


def _match_virus_values(row) -> bool:
    left = _normalize_virus_text(row.get("virus_identified"))
    right = _normalize_virus_text(row.get("Virus_PCR"))
    if not left or not right:
        return False
    return not left.isdisjoint(right)


###############################################

input_path = snakemake.input["metadata"]
output_path_sequencing = snakemake.output["metadata_sequencing"]
output_path_pcr = snakemake.output["metadata_pcr"]


d = pd.read_csv(input_path, sep="\t")

# make Sampled_Date/date column consistent in yyyy-mm-dd format
if "Sampled_Date" in d.columns and "date" in d.columns:
    # Both exist: prefer Sampled_Date, fallback to date
    date_source = d["Sampled_Date"].fillna(d["date"])
elif "Sampled_Date" in d.columns:
    date_source = d["Sampled_Date"]
elif "date" in d.columns:
    date_source = d["date"]
else:
    date_source = pd.Series([pd.NA] * len(d))

date_source = date_source.astype(str).str.strip()
date_source = date_source.replace({"": pd.NA, "nan": pd.NA, "NaT": pd.NA})

parsed_dates = pd.to_datetime(date_source, format="%Y-%m-%d", errors="coerce")
parsed_dates = parsed_dates.fillna(
    pd.to_datetime(date_source, format="%d.%m.%Y", errors="coerce")
)
parsed_dates = parsed_dates.fillna(pd.to_datetime(date_source, errors="coerce"))

# ensure final column is named "date"
d = d.drop(columns=[c for c in ["Sampled_Date", "date"] if c in d.columns])
d["date"] = parsed_dates.dt.strftime("%Y-%m-%d")

# determine whether virus_identified and Virus_PCR match, output match in new column Match_PCR
d["Match_PCR"] = d.apply(_match_virus_values, axis=1)

# if SampleID exists, we create a variable called strain, which is of the format {SampleID}|{Sampled_Date}|{Canton}
if "SampleID" in d.columns:
    # Use Canton if available, fallback to location
    location_value = (
        d["Canton"] if "Canton" in d.columns else pd.Series([pd.NA] * len(d))
    )
    if "location" in d.columns:
        location_value = location_value.fillna(d["location"])

    d["strain"] = (
        d["SampleID"].astype(str)
        + "|"
        + parsed_dates.dt.strftime("%Y-%m-%d").fillna("")
        + "|"
        + location_value.astype(str)
    )

# Rename columns to standardized names
if "virus_subtype" in d.columns and "virus_identified" not in d.columns:
    d.rename(columns={"virus_subtype": "virus_identified"}, inplace=True)
if "Canton" in d.columns and "location" not in d.columns:
    d.rename(columns={"Canton": "location"}, inplace=True)

# select only relevant columns: SampleID, virus_identified, date, location, percentage_of_genome_covered
cols = [
    "strain",
    "virus_identified",
    "date",
    "location",
    # "percentage_of_genome_covered",
    "Match_PCR",
]

# select only columns that exist
sequencing_cols = [col for col in cols if col in d.columns]
d_sequencing = d[sequencing_cols]

d_sequencing.to_csv(output_path_sequencing, sep="\t", index=False)


#################################
# PCR metadata processing
#################################

if "Virus_PCR" in d.columns and d["Virus_PCR"].notna().any():
    pcr_cols = [
        col
        for col in ["strain", "date", "location", "virus_identified"]
        if col in d.columns
    ]
    d_pcr = d[pcr_cols].copy()

    # Unlist Virus_PCR column: split by "," and create separate rows
    d_pcr["virus_identified_pcr"] = d["Virus_PCR"].fillna("")
    d_pcr["virus_identified_pcr"] = d_pcr["virus_identified_pcr"].str.split(",")

    # Explode to create one row per virus in Virus_PCR
    d_pcr = d_pcr.explode("virus_identified_pcr", ignore_index=True)

    # Clean up the virus_identified_pcr values (strip whitespace)
    d_pcr["virus_identified_pcr"] = d_pcr["virus_identified_pcr"].str.strip()

    # Remove empty rows
    d_pcr = d_pcr[d_pcr["virus_identified_pcr"] != ""]

    # Compare virus_identified_pcr with virus_identified using normalization logic
    def _match_sequencing(row) -> bool:
        left = _normalize_virus_text(row.get("virus_identified"))
        right = _normalize_virus_text(row.get("virus_identified_pcr"))
        if not left or not right:
            return False
        return not left.isdisjoint(right)

    d_pcr["Match_Sequencing"] = d_pcr.apply(_match_sequencing, axis=1)

    # Select and reorder columns
    d_pcr = d_pcr[
        [
            col
            for col in [
                "strain",
                "date",
                "location",
                "Match_Sequencing",
                "virus_identified_pcr",
            ]
            if col in d_pcr.columns
        ]
    ]
else:
    # Create empty PCR metadata if Virus_PCR doesn't exist
    d_pcr = pd.DataFrame(
        columns=[
            "strain",
            "date",
            "location",
            "Match_Sequencing",
            "virus_identified_pcr",
        ]
    )

d_pcr.to_csv(output_path_pcr, sep="\t", index=False)
