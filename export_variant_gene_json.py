import json
from pathlib import Path

import pandas as pd

BASE_DIR = Path(__file__).resolve().parent
INPUT_CSV = BASE_DIR / "variant_gene_edges.csv"
OUTPUT_JSON = BASE_DIR / "variant_gene.json"


def clean_value(value):
    if pd.isna(value):
        return None
    if isinstance(value, str):
        value = value.strip()
        return value if value else None
    return value


def parse_pipe_list(value):
    value = clean_value(value)
    if value is None:
        return []

    parts = [item.strip() for item in str(value).split("|") if item.strip()]
    return parts


def parse_bool(value):
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    return text in {"true", "1", "yes", "y"}


def build_node_id(node_type, raw_id):
    return f"{node_type}:{str(raw_id).strip()}"


def main():
    if not INPUT_CSV.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_CSV}")

    df = pd.read_csv(INPUT_CSV)

    required_columns = {
        "source", "source_type",
        "target", "target_type",
        "variant_id", "variant_type", "rsid", "haplotype",
        "hgnc_name", "gene_family", "primary_gene_family",
        "gene_family_id", "gene_family_match_source",
        "in_clinical_variants_tsv", "in_interactions_tsv", "in_relationships_tsv"
    }

    missing = required_columns - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    nodes = {}
    links = []

    for i, row in df.iterrows():
        source = clean_value(row["source"])
        source_type = clean_value(row["source_type"])
        target = clean_value(row["target"])
        target_type = clean_value(row["target_type"])

        if source is None or target is None or source_type is None or target_type is None:
            continue

        source_type = str(source_type).lower()
        target_type = str(target_type).lower()

        source_id = build_node_id(source_type, source)
        target_id = build_node_id(target_type, target)

        rsid = clean_value(row["rsid"])
        haplotype = clean_value(row["haplotype"])
        variant_id = clean_value(row["variant_id"])
        variant_type = clean_value(row["variant_type"])
        hgnc_name = clean_value(row["hgnc_name"])
        gene_family = parse_pipe_list(row["gene_family"])
        primary_gene_family = clean_value(row["primary_gene_family"])
        gene_family_ids = parse_pipe_list(row["gene_family_id"])
        gene_family_match_source = clean_value(row["gene_family_match_source"])

        in_clinical_variants_tsv = parse_bool(row["in_clinical_variants_tsv"])
        in_interactions_tsv = parse_bool(row["in_interactions_tsv"])
        in_relationships_tsv = parse_bool(row["in_relationships_tsv"])

        # -----------------------------
        # Variant node
        # -----------------------------
        if source_type == "variant" and source_id not in nodes:
            nodes[source_id] = {
                "id": source_id,
                "label": rsid if rsid is not None else str(source),
                "type": "variant",
                "group": "variant",
                "variant_key": source,
                "variant_id": variant_id,
                "variant_type": variant_type,
                "rsid": rsid,
                "haplotype": haplotype,
                "in_clinical_variants_tsv": in_clinical_variants_tsv,
                "in_interactions_tsv": in_interactions_tsv,
                "in_relationships_tsv": in_relationships_tsv,
            }

        # -----------------------------
        # Gene node
        # -----------------------------
        if target_type == "gene" and target_id not in nodes:
            nodes[target_id] = {
                "id": target_id,
                "label": str(target),
                "type": "gene",
                "group": "gene",
                "gene_symbol": str(target),
                "hgnc_name": hgnc_name,
                "gene_family": gene_family,
                "primary_gene_family": primary_gene_family,
                "gene_family_id": gene_family_ids,
                "gene_family_match_source": gene_family_match_source,
            }

        # fallback if file direction is reversed in future
        if source_type == "gene" and source_id not in nodes:
            nodes[source_id] = {
                "id": source_id,
                "label": str(source),
                "type": "gene",
                "group": "gene",
                "gene_symbol": str(source),
                "hgnc_name": hgnc_name,
                "gene_family": gene_family,
                "primary_gene_family": primary_gene_family,
                "gene_family_id": gene_family_ids,
                "gene_family_match_source": gene_family_match_source,
            }

        if target_type == "variant" and target_id not in nodes:
            nodes[target_id] = {
                "id": target_id,
                "label": rsid if rsid is not None else str(target),
                "type": "variant",
                "group": "variant",
                "variant_key": target,
                "variant_id": variant_id,
                "variant_type": variant_type,
                "rsid": rsid,
                "haplotype": haplotype,
                "in_clinical_variants_tsv": in_clinical_variants_tsv,
                "in_interactions_tsv": in_interactions_tsv,
                "in_relationships_tsv": in_relationships_tsv,
            }

        # -----------------------------
        # Link
        # -----------------------------
        links.append({
            "id": f"link_{i}",
            "source": source_id,
            "target": target_id,
            "relation": f"{source_type}-{target_type}",
            "source_type": source_type,
            "target_type": target_type,
            "variant_id": variant_id,
            "variant_type": variant_type,
            "rsid": rsid,
            "haplotype": haplotype,
            "gene_symbol": str(target) if target_type == "gene" else (str(source) if source_type == "gene" else None),
            "hgnc_name": hgnc_name,
            "gene_family": gene_family,
            "primary_gene_family": primary_gene_family,
            "gene_family_id": gene_family_ids,
            "gene_family_match_source": gene_family_match_source,
            "in_clinical_variants_tsv": in_clinical_variants_tsv,
            "in_interactions_tsv": in_interactions_tsv,
            "in_relationships_tsv": in_relationships_tsv,
        })

    graph_json = {
        "metadata": {
            "graph_type": "variant-gene knowledge graph",
            "input_file": INPUT_CSV.name,
            "node_count": len(nodes),
            "link_count": len(links),
        },
        "nodes": list(nodes.values()),
        "links": links,
    }

    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump(graph_json, f, indent=2, ensure_ascii=False)

    print(f"Saved JSON to: {OUTPUT_JSON}")
    print(f"Nodes: {len(nodes)}")
    print(f"Links: {len(links)}")


if __name__ == "__main__":
    main()