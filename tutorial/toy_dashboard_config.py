# Dashboard configuration for the tutorial's toy dataset.
# Activate it with:  export REVSEQ_CONFIG=tutorial/toy_dashboard_config.py
# Any setting left out keeps its default from dashboard/config.py.

SITE_TITLE = "ReVSeq tutorial dashboard"
REGION_NAME = "Magic Land (Switzerland)"

# Modules shown on every virus page unless a page overrides them.
DEFAULT_MODULES = {
    "barplot_pcr": True,
    "barplot_sequencing": True,
    "map": True,
    "pileup": True,
    "tree": True,
}

# One entry per virus page. `data_name` must match the harmonized virus name (see STRAIN_ALIASES).
STRAIN_CONFIG = {
    "RSV": {
        "label": "RSV - A/B",
        "data_name": "RSV - A/B",
        # Substrains (RSV-A, RSV-B) are found from rsv-a_avg.json / rsv-b_avg.json.
        "trees": [
            {"title": "RSV A", "dataset": "RSV-A"},
            {"title": "RSV B", "dataset": "RSV-B"},
        ],
    },
    "Influenza_A": {
        "label": "Influenza A",
        "data_name": "Influenza A",
        "pileup_data_prefix": "flu-a",  # finds flu-a-h3n2_avg.json (subtype H3N2)
        "pileup_segments": [
            {"value": segment, "label": segment, "annotation": f"Influenza_A_{segment}.gb"}
            for segment in ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
        ],
        "pileup_default_segment": "HA",
        "trees": [{"title": "H3N2 segment 4 (HA)", "dataset": "Influenza-A-H3N2-HA"}],
    },
    "Parainfluenza_3": {
        "label": "Parainfluenza 3",
        "data_name": "Parainfluenza 3",
        "pileup_data_prefix": "hpiv-3",
        "trees": [{"title": "", "dataset": "HPIV-3"}],
    },
}

COLOR_BY_STRAIN = {
    "RSV - A/B": "#7B61C9",
    "Influenza A": "#2A9D8F",
    "Parainfluenza 3": "#E9A23B",
}

# Extra spellings of virus names in your metadata -> harmonized name (added to the defaults).
STRAIN_ALIASES = {}  # e.g. {"hRSV-A": "RSV - A/B"}

# Map: the toy samples are Swiss, so we keep the bundled canton map.
# For another country, point `geojson` to your file and `id_property` to the feature
# property holding the codes used in your metadata `location` column.
MAP_CONFIG = {
    "geojson": "dashboard/static/swiss_cantons.geojson",
    "id_property": "canton",
    "name_property": "NAME",
    "center": {"lat": 46.8, "lon": 8.3},
    "zoom": 5.5,
    "style": "open-street-map",
}

# Extra location spellings -> region code (added to the defaults). Case and accents are ignored.
LOCATION_ALIASES = {"Kanton Bern": "BE"}
