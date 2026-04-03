from pathlib import Path

import pandas as pd
from sqlalchemy import create_engine


DB_CONFIG = {
    "host": "localhost",
    "database": "pharmgraph",
    "user": "maywelkin",
    "password": "",
    "port": 5432,
}


def get_engine():
    conn_string = (
        f"postgresql://{DB_CONFIG['user']}:{DB_CONFIG['password']}"
        f"@{DB_CONFIG['host']}:{DB_CONFIG['port']}/{DB_CONFIG['database']}"
    )
    return create_engine(conn_string)


def normalize_symbol(value):
    if pd.isna(value):
        return None
    value = str(value).strip()
    return value if value else None


def split_pipe_values(value):
    if pd.isna(value):
        return []
    return [item.strip() for item in str(value).split("|") if item and item.strip()]


def load_hgnc_family_lookup(hgnc_path):
    usecols = [
        "symbol",
        "name",
        "gene_group",
        "gene_group_id",
        "alias_symbol",
        "prev_symbol",
    ]

    hgnc = pd.read_csv(hgnc_path, sep="\t", dtype=str, usecols=usecols).fillna("")

    exact = hgnc[["symbol", "name", "gene_group", "gene_group_id"]].copy()
    exact["match_symbol"] = exact["symbol"].str.strip()
    exact["match_source"] = "symbol"
    exact["priority"] = 0

    prev = hgnc[["name", "gene_group", "gene_group_id", "prev_symbol"]].copy()
    prev["match_symbol"] = prev["prev_symbol"].str.split("|")
    prev = prev.explode("match_symbol")
    prev["match_symbol"] = prev["match_symbol"].fillna("").str.strip()
    prev["match_source"] = "prev_symbol"
    prev["priority"] = 1

    alias = hgnc[["name", "gene_group", "gene_group_id", "alias_symbol"]].copy()
    alias["match_symbol"] = alias["alias_symbol"].str.split("|")
    alias = alias.explode("match_symbol")
    alias["match_symbol"] = alias["match_symbol"].fillna("").str.strip()
    alias["match_source"] = "alias_symbol"
    alias["priority"] = 2

    prev = prev.rename(columns={"prev_symbol": "symbol"})
    alias = alias.rename(columns={"alias_symbol": "symbol"})

    lookup = pd.concat([exact, prev, alias], ignore_index=True, sort=False)
    lookup = lookup[lookup["match_symbol"] != ""].copy()

    lookup = (
        lookup.sort_values(["match_symbol", "priority"])
        .drop_duplicates(subset=["match_symbol"], keep="first")
        .rename(
            columns={
                "name": "hgnc_name",
                "gene_group": "gene_family",
                "gene_group_id": "gene_family_id",
            }
        )
    )

    return lookup[
        [
            "match_symbol",
            "match_source",
            "hgnc_name",
            "gene_family",
            "gene_family_id",
        ]
    ]


def load_reference_gene_sets(base_dir):
    gene_sets = {
        "in_clinical_variants_tsv": set(),
        "in_interactions_tsv": set(),
        "in_relationships_tsv": set(),
    }

    clinical_path = base_dir / "clinicalVariants.tsv"
    if clinical_path.exists():
        clinical_df = pd.read_csv(clinical_path, sep="\t", dtype=str, usecols=["gene"])
        gene_sets["in_clinical_variants_tsv"] = {
            normalize_symbol(x) for x in clinical_df["gene"] if normalize_symbol(x)
        }

    interactions_path = base_dir / "interactions.tsv"
    if interactions_path.exists():
        interactions_df = pd.read_csv(
            interactions_path,
            sep="\t",
            dtype=str,
            usecols=["gene_name", "gene_claim_name"],
        )
        interaction_genes = set()
        for col in ["gene_name", "gene_claim_name"]:
            interaction_genes.update(
                normalize_symbol(x)
                for x in interactions_df[col]
                if normalize_symbol(x)
            )
        gene_sets["in_interactions_tsv"] = interaction_genes

    relationships_path = base_dir / "relationships.tsv"
    if relationships_path.exists():
        relationships_df = pd.read_csv(
            relationships_path,
            sep="\t",
            dtype=str,
            usecols=["Entity1_name", "Entity1_type"],
        )
        relationships_df = relationships_df[
            relationships_df["Entity1_type"].fillna("").str.strip().str.lower() == "gene"
        ]
        gene_sets["in_relationships_tsv"] = {
            normalize_symbol(x)
            for x in relationships_df["Entity1_name"]
            if normalize_symbol(x)
        }

    return gene_sets


def main():
    base_dir = Path(__file__).resolve().parent
    hgnc_path = base_dir / "hgnc_complete_set.tsv"

    if not hgnc_path.exists():
        raise FileNotFoundError(f"Could not find HGNC file: {hgnc_path}")

    engine = get_engine()

    query = """
    SELECT
        COALESCE(v.rsid, v.haplotype) AS source,
        'variant' AS source_type,
        g.gene_symbol AS target,
        'gene' AS target_type,
        v.variant_id,
        v.variant_type,
        v.rsid,
        v.haplotype
    FROM variants v
    JOIN genes g ON v.gene_id = g.gene_id
    WHERE v.rsid IS NOT NULL OR v.haplotype IS NOT NULL
    """

    edges_df = pd.read_sql(query, engine)

    print("Number of variant-gene edges:", len(edges_df))
    print(edges_df.head())

    edges_df["target"] = edges_df["target"].apply(normalize_symbol)
    edges_df["source"] = edges_df["source"].apply(normalize_symbol)

    hgnc_lookup = load_hgnc_family_lookup(hgnc_path)

    edges_df = edges_df.merge(
        hgnc_lookup,
        left_on="target",
        right_on="match_symbol",
        how="left",
    )

    edges_df["gene_family"] = edges_df["gene_family"].fillna("Unknown")
    edges_df["gene_family_id"] = edges_df["gene_family_id"].fillna("")
    edges_df["hgnc_name"] = edges_df["hgnc_name"].fillna("")
    edges_df["gene_family_match_source"] = edges_df["match_source"].fillna("unmatched")

    edges_df["primary_gene_family"] = edges_df["gene_family"].apply(
        lambda x: split_pipe_values(x)[0] if x != "Unknown" and split_pipe_values(x) else "Unknown"
    )

    reference_gene_sets = load_reference_gene_sets(base_dir)
    for col_name, gene_set in reference_gene_sets.items():
        edges_df[col_name] = edges_df["target"].isin(gene_set)

    edges_df = edges_df[
        [
            "source",
            "source_type",
            "target",
            "target_type",
            "variant_id",
            "variant_type",
            "rsid",
            "haplotype",
            "hgnc_name",
            "gene_family",
            "primary_gene_family",
            "gene_family_id",
            "gene_family_match_source",
            "in_clinical_variants_tsv",
            "in_interactions_tsv",
            "in_relationships_tsv",
        ]
    ]

    output_file = base_dir / "variant_gene_edges.csv"
    edges_df.to_csv(output_file, index=False)

    print(f"Saved to {output_file}")
    print("Number of unique variants:", edges_df["source"].nunique())
    print("Number of unique genes:", edges_df["target"].nunique())
    print("Number of nodes:", edges_df["source"].nunique() + edges_df["target"].nunique())
    print("Number of edges in graph:", len(edges_df))

    matched_rows = (edges_df["gene_family_match_source"] != "unmatched").sum()
    print("Rows matched to HGNC gene family:", matched_rows)
    print("Rows unmatched to HGNC gene family:", len(edges_df) - matched_rows)

    print("\nTop 10 genes with most variants:")
    top_genes = edges_df["target"].value_counts().head(10)
    for gene, deg in top_genes.items():
        print(f"{gene}: {deg}")

    print("\nTop 10 primary gene families:")
    top_families = edges_df["primary_gene_family"].value_counts().head(10)
    for family, count in top_families.items():
        print(f"{family}: {count}")


if __name__ == "__main__":
    main()