from copy import deepcopy
import re
import unicodedata

DEFAULT_MODULES = {
    "barplot_pcr": True,
    "barplot_sequencing": True,
    "map": True,
    "pileup": True,
    "tree": False,
}

DEFAULT_PILEUP_LEVELS = ["all", "substrain", "individual"]
DEFAULT_PILEUP_MAX_INDIVIDUAL_TRACES = 30
MIXED_PAGE_ENABLED = True

STRAIN_CONFIG = {
    "Influenza_A": {
        "label": "Influenza A",
        "data_name": "Influenza A",
        "modules": {"tree": True},
        "pileup_data_prefix": "flu-a",
        "pileup_segments": [
            {"value": "PB2", "label": "PB2", "annotation": "Influenza_A_PB2.gb"},
            {"value": "PB1", "label": "PB1", "annotation": "Influenza_A_PB1.gb"},
            {"value": "PA", "label": "PA", "annotation": "Influenza_A_PA.gb"},
            {"value": "HA", "label": "HA", "annotation": "Influenza_A_HA.gb"},
            {"value": "NP", "label": "NP", "annotation": "Influenza_A_NP.gb"},
            {"value": "NA", "label": "NA", "annotation": "Influenza_A_NA.gb"},
            {"value": "MP", "label": "MP", "annotation": "Influenza_A_MP.gb"},
            {"value": "NS", "label": "NS", "annotation": "Influenza_A_NS.gb"},
        ],
        "pileup_default_segment": "HA",
        "trees": [
            {"title": "H1N1 segment 4 (HA)", "dataset": "Influenza-A-H1N1-HA"},
            {"title": "H3N2 segment 4 (HA)", "dataset": "Influenza-A-H3N2-HA"},
        ],
    },
    "Influenza_B": {
        "label": "Influenza B",
        "data_name": "Influenza B",
        "modules": {"tree": True},
        "pileup_data_prefix": "flu-b",
        "pileup_segments": [
            {"value": "PB2", "label": "PB2", "annotation": "Influenza_B_PB2.gb"},
            {"value": "PB1", "label": "PB1", "annotation": "Influenza_B_PB1.gb"},
            {"value": "HA", "label": "HA", "annotation": "Influenza_B_HA.gb"},
            {"value": "NP", "label": "NP", "annotation": "Influenza_B_NP.gb"},
            {"value": "NA", "label": "NA", "annotation": "Influenza_B_NA.gb"},
            {"value": "MP", "label": "MP", "annotation": "Influenza_B_MP.gb"},
            {"value": "NS", "label": "NS", "annotation": "Influenza_B_NS.gb"},
        ],
        "pileup_default_segment": "HA",
        "trees": [
            {"title": "Influenza B segment 4 (HA)", "dataset": "Influenza-B"},
        ],
    },
    "Metapneumovirus": {
        "label": "Metapneumovirus",
        "data_name": "Metapneumovirus",
        "modules": {"tree": True},
        "pileup_data_prefix": "hmpv",
        "trees": [
            {"title": "", "dataset": "HMPV"},
        ],
    },
    "Parainfluenza_1": {
        "label": "Parainfluenza 1",
        "data_name": "Parainfluenza 1",
        "modules": {"tree": True},
        "pileup_data_prefix": "hpiv-1",
        "trees": [{"title": "", "dataset": "HPIV-1"}],
    },
    "Parainfluenza_2": {
        "label": "Parainfluenza 2",
        "data_name": "Parainfluenza 2",
        "modules": {"tree": True},
        "pileup_data_prefix": "hpiv-2",
        "trees": [{"title": "", "dataset": "HPIV-2"}],
    },
    "Parainfluenza_3": {
        "label": "Parainfluenza 3",
        "data_name": "Parainfluenza 3",
        "modules": {"tree": True},
        "pileup_data_prefix": "hpiv-3",
        "trees": [{"title": "", "dataset": "HPIV-3"}],
    },
    "Parainfluenza_4a": {
        "label": "Parainfluenza 4a",
        "data_name": "Parainfluenza 4a",
        "modules": {"tree": True},
        "pileup_data_prefix": "hpiv-4a",
        "trees": [{"title": "", "dataset": "HPIV-4a"}],
    },
    "RSV": {
        "label": "RSV - A/B",
        "data_name": "RSV - A/B",
        "modules": {"tree": True},
        "trees": [
            {"title": "RSV A", "dataset": "RSV-A"},
            {"title": "RSV B", "dataset": "RSV-B"},
        ],
    },
    "SARS-CoV-2": {
        "label": "SARS-CoV-2",
        "data_name": "SARS-CoV-2",
        "modules": {"tree": True},
        "pileup_data_prefix": "sars-cov-2",
        "trees": [{"title": "", "dataset": "SARS-CoV-2"}],
    },
    "coronavirus_229E": {
        "label": "coronavirus 229E",
        "data_name": "coronavirus 229E",
        "modules": {"tree": True},
        "pileup_data_prefix": "cov-229e",
        "trees": [{"title": "", "dataset": "Coronavirus-229E"}],
    },
    "coronavirus_HKU1": {
        "label": "coronavirus HKU1",
        "data_name": "coronavirus HKU1",
        "modules": {"tree": True},
        "pileup_data_prefix": "cov-hku1",
        "trees": [{"title": "", "dataset": "Coronavirus-HKU1"}],
    },
    "coronavirus_NL63": {
        "label": "coronavirus NL63",
        "data_name": "coronavirus NL63",
        "modules": {"tree": True},
        "pileup_data_prefix": "cov-nl63",
        "trees": [{"title": "", "dataset": "Coronavirus-NL63"}],
    },
    "coronavirus_OC43": {
        "label": "coronavirus OC43",
        "data_name": "coronavirus OC43",
        "modules": {"tree": True},
        "pileup_data_prefix": "cov-oc43",
        "trees": [{"title": "", "dataset": "Coronavirus-OC43"}],
    },
}

STRAIN_ORDER = list(STRAIN_CONFIG.keys())
DASHBOARD_STRAIN_NAMES = {item["data_name"] for item in STRAIN_CONFIG.values()}
SLUG_BY_DATA_NAME = {item["data_name"]: slug for slug, item in STRAIN_CONFIG.items()}

COLOR_BY_STRAIN = {
    "Influenza A": "#6875D9",
    "Influenza B": "#E158D3",
    "Metapneumovirus": "royalblue",
    "Parainfluenza 2": "beige",
    "Parainfluenza 3": "#848877",
    "Parainfluenza 4a": "#FDBF6F",
    "Polyomavirus": "#D9664D",
    "RSV - A/B": "thistle",
    "SARS-CoV-2": "#E9D588",
    "Rhino- / Enterovirus": "#74457B",
    "coronavirus 229E": "pink",
    "coronavirus HKU1": "#BAE6B0",
    "coronavirus NL63": "#A968E1",
    "coronavirus OC43": "#65E0C8",
}

# Source-specific normalization maps.
_BASE_CANONICAL_BY_KEY = {
    "adenovirus": "Adenovirus",
    "bocavirus": "Bocavirus",
    "metapneumovirus": "Metapneumovirus",
    "human metapneumovirus": "Metapneumovirus",
    "human metapneumovirus a": "Metapneumovirus",
    "human metapneumovirus b": "Metapneumovirus",
    "parainfluenza 1": "Parainfluenza 1",
    "human parainfluenza 1": "Parainfluenza 1",
    "human parainfluenza virus 1": "Parainfluenza 1",
    "parainfluenza 2": "Parainfluenza 2",
    "parainfluenza 3": "Parainfluenza 3",
    "parainfluenza 4a": "Parainfluenza 4a",
    "human parainfluenza 2": "Parainfluenza 2",
    "human parainfluenza 3": "Parainfluenza 3",
    "human parainfluenza virus 2": "Parainfluenza 2",
    "human parainfluenza virus 3": "Parainfluenza 3",
    "human parainfluenza virus 4a": "Parainfluenza 4a",
    "Parainfluenza 4": "Parainfluenza 4a",
    "influenza a": "Influenza A",
    "influenza a virus": "Influenza A",
    "influenza b": "Influenza B",
    "influenza b virus": "Influenza B",
    "influenza b virus (b/brisbane/60/2008)": "Influenza B",
    "rsv - a/b": "RSV - A/B",
    "rsv a": "RSV - A/B",
    "rsv b": "RSV - A/B",
    "respiratory syncytial virus (type a)": "RSV - A/B",
    "human respiratory syncytial virus 9320 (type b)": "RSV - A/B",
    "rhino- / enterovirus": "Rhino- / Enterovirus",
    "rhino-/enterovirus": "Rhino- / Enterovirus",
    "sars-cov-2": "SARS-CoV-2",
    "coronavirus 229e": "coronavirus 229E",
    "coronavirus hku1": "coronavirus HKU1",
    "coronavirus nl63": "coronavirus NL63",
    "coronavirus oc43": "coronavirus OC43",
    "human coronavirus 229e": "coronavirus 229E",
    "human coronavirus hku1": "coronavirus HKU1",
    "human coronavirus nl63": "coronavirus NL63",
    "human coronavirus oc43": "coronavirus OC43",
}


def _normalize_key(value):
    if value is None:
        return ""
    value = str(value).strip()
    value = re.sub(r"\s+", " ", value)
    return value.lower()


def _normalize_location_key(value):
    if value is None:
        return ""
    value = str(value).strip()
    value = (
        unicodedata.normalize("NFKD", value).encode("ascii", "ignore").decode("ascii")
    )
    value = re.sub(r"\s+", " ", value)
    return value.lower()


_CANTON_BY_KEY = {
    "ag": "AG",
    "ai": "AI",
    "ar": "AR",
    "be": "BE",
    "bl": "BL",
    "bs": "BS",
    "fr": "FR",
    "ge": "GE",
    "gl": "GL",
    "gr": "GR",
    "ju": "JU",
    "lu": "LU",
    "ne": "NE",
    "nw": "NW",
    "ow": "OW",
    "sg": "SG",
    "sh": "SH",
    "so": "SO",
    "sz": "SZ",
    "tg": "TG",
    "ti": "TI",
    "ur": "UR",
    "vd": "VD",
    "vs": "VS",
    "zg": "ZG",
    "zh": "ZH",
    "aargau": "AG",
    "appenzell innerrhoden": "AI",
    "appenzell ausserrhoden": "AR",
    "bern": "BE",
    "berne": "BE",
    "basel-landschaft": "BL",
    "basel landschaft": "BL",
    "basel-stadt": "BS",
    "basel stadt": "BS",
    "fribourg": "FR",
    "freiburg": "FR",
    "geneve": "GE",
    "genf": "GE",
    "glarus": "GL",
    "graubunden": "GR",
    "grisons": "GR",
    "jura": "JU",
    "lucerne": "LU",
    "luzern": "LU",
    "neuchatel": "NE",
    "nidwalden": "NW",
    "obwalden": "OW",
    "st gallen": "SG",
    "st. gallen": "SG",
    "schaffhausen": "SH",
    "solothurn": "SO",
    "schwyz": "SZ",
    "thurgau": "TG",
    "ticino": "TI",
    "uri": "UR",
    "vaud": "VD",
    "valais": "VS",
    "wallis": "VS",
    "zug": "ZG",
    "zurich": "ZH",
    "zurich canton": "ZH",
}


def harmonize_detected_strain(raw_value, source: str):
    # Returns: (canonical_strain, display_label). Unknown values are dropped upstream.
    key = _normalize_key(raw_value)
    if not key:
        return None, None

    canonical = _BASE_CANONICAL_BY_KEY.get(key)
    if canonical is None:
        if "influenza a" in key and ("h1n1" in key or "h3n2" in key):
            canonical = "Influenza A"
        else:
            return None, None

    if canonical not in DASHBOARD_STRAIN_NAMES:
        return None, None

    if source == "pcr":
        # PCR legend should always use canonical labels.
        return canonical, canonical

    # Sequencing keeps selected subtype display labels.
    if canonical == "RSV - A/B":
        if key in {"rsv a", "respiratory syncytial virus (type a)"}:
            return canonical, "RSV-A"
        if key in {"rsv b", "human respiratory syncytial virus 9320 (type b)"}:
            return canonical, "RSV-B"
    if canonical == "Influenza A":
        if "h1n1" in key:
            return canonical, "H1N1"
        if "h3n2" in key:
            return canonical, "H3N2"
    if canonical == "Parainfluenza 2":
        if "human parainfluenza" in key:
            return canonical, "Human parainfluenza 2"
    if canonical == "Parainfluenza 3":
        if "human parainfluenza" in key:
            return canonical, "Human parainfluenza 3"
    if canonical == "Parainfluenza 4a":
        if "human parainfluenza" in key:
            return canonical, "Human parainfluenza virus 4a"

    return canonical, canonical


def harmonize_canton(raw_value):
    key = _normalize_location_key(raw_value)
    if not key:
        return None
    return _CANTON_BY_KEY.get(key, key.upper() if len(key) == 2 else None)


def get_strain_config(slug):
    config = STRAIN_CONFIG.get(slug)
    if config is None:
        return None

    result = deepcopy(config)
    result["slug"] = slug
    raw_modules = result.get("modules", {})
    modules = DEFAULT_MODULES.copy()
    modules.update(raw_modules)
    result["modules"] = modules
    result["pileup_levels"] = result.get("pileup_levels", DEFAULT_PILEUP_LEVELS.copy())
    result["pileup_max_individual_traces"] = result.get(
        "pileup_max_individual_traces", DEFAULT_PILEUP_MAX_INDIVIDUAL_TRACES
    )
    result["pileup_annotation"] = result.get("pileup_annotation", f"{slug}.gb")
    result["pileup_data_prefix"] = result.get("pileup_data_prefix", slug)
    return result


def get_strain_options():
    options = [{"value": "/home", "label": "All strains"}]
    options.extend(
        {
            "value": f"/strain/{slug}",
            "label": STRAIN_CONFIG[slug]["label"],
        }
        for slug in STRAIN_ORDER
    )
    if MIXED_PAGE_ENABLED:
        options.append({"value": "/mixed", "label": "Mixed"})
    return options


def is_mixed_page_enabled() -> bool:
    return bool(MIXED_PAGE_ENABLED)
