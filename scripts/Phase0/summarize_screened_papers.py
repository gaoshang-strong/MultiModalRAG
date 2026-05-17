"""
Generate a markdown summary of screened paper results from screen_similar_data_papers.py.
No API calls — reads CSVs from literature/00_similar_data_paper_review/ directly.

Usage:
    python summarize_screened_papers.py
    python summarize_screened_papers.py --min-score 5 --top-n 30
    python summarize_screened_papers.py --min-if 10 --top-n-if 10
"""

import argparse
import os
from collections import Counter
from datetime import date

import pandas as pd

INPUT_DIR = "/ShangGaoAIProjects/Gastric/genome/literature/00_similar_data_paper_review"
OUTPUT_DIR = INPUT_DIR
JIF_CSV = "/ShangGaoAIProjects/Gastric/genome/literature/Journal_impact_factor.csv"


def parse_args():
    parser = argparse.ArgumentParser(description="Summarize screened paper results as markdown.")
    parser.add_argument("--min-score", type=float, default=4.0,
                        help="Minimum similarity_score to include in top paper list (default: 4.0)")
    parser.add_argument("--top-n", type=int, default=30,
                        help="Maximum number of top papers to list by similarity score (default: 30)")
    parser.add_argument("--min-if", type=float, default=5.0,
                        help="Minimum journal impact factor for the high-IF section (default: 5.0)")
    parser.add_argument("--top-n-if", type=int, default=20,
                        help="Maximum papers in the high-IF section (default: 20)")
    return parser.parse_args()


def load_combined(input_dir: str) -> pd.DataFrame:
    combined_path = os.path.join(input_dir, "combined_ranked.csv")
    if os.path.isfile(combined_path):
        return pd.read_csv(combined_path, dtype=str).fillna("")
    # Fall back to merging individual screened CSVs
    dfs = []
    for fname in os.listdir(input_dir):
        if fname.endswith("_screened.csv"):
            dfs.append(pd.read_csv(os.path.join(input_dir, fname), dtype=str).fillna(""))
    if not dfs:
        raise FileNotFoundError(f"No screened CSVs found in {input_dir}")
    combined = pd.concat(dfs, ignore_index=True).drop_duplicates(subset=["pmid"])
    combined["similarity_score"] = pd.to_numeric(combined["similarity_score"], errors="coerce").fillna(0)
    return combined.sort_values("similarity_score", ascending=False).reset_index(drop=True)


def load_jif_map(jif_csv: str) -> dict[str, float]:
    """Load journal -> 2024 JIF from Journal_impact_factor.csv. Returns empty dict if file missing."""
    if not os.path.isfile(jif_csv):
        return {}
    df = pd.read_csv(jif_csv, dtype=str).fillna("")
    result = {}
    for _, row in df.iterrows():
        journal = row.get("journal", "").strip()
        val = row.get("impact_factor_2024", "")
        try:
            result[journal] = float(val)
        except (ValueError, TypeError):
            pass
    return result


def safe_int(val):
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return None


def fmt_list(val: str, sep: str = "; ") -> str:
    if not val or val.strip() == "":
        return "_not reported_"
    items = [v.strip() for v in val.split(sep) if v.strip()]
    return " · ".join(items) if items else "_not reported_"


def score_bar(score: float, max_score: float = 10.0, width: int = 10) -> str:
    filled = round(score / max_score * width)
    return "█" * filled + "░" * (width - filled)


def generate_summary(df: pd.DataFrame, min_score: float, top_n: int,
                     jif_map: dict, min_if: float, top_n_if: int) -> str:
    df["similarity_score"] = pd.to_numeric(df["similarity_score"], errors="coerce").fillna(0)
    df["sample_size_int"] = df["sample_size"].apply(safe_int)

    total = len(df)
    similar = (df["similar_to_scenario"].str.lower() == "true").sum()
    has_germline = (df["has_germline"].str.lower() == "true").sum()
    has_matched = (df["has_matched_normal"].str.lower() == "true").sum()
    has_both = (
        (df["has_germline"].str.lower() == "true") &
        (df["has_matched_normal"].str.lower() == "true")
    ).sum()

    score_bins = pd.cut(
        df["similarity_score"],
        bins=[0, 2, 4, 6, 8, 10],
        labels=["0–2", "2–4", "4–6", "6–8", "8–10"],
        include_lowest=True,
    )
    score_dist = score_bins.value_counts().sort_index()

    cancer_counts = Counter(
        v.strip().lower()
        for v in df["cancer_type"].dropna()
        if v.strip()
    )
    seq_counts = Counter(
        v.strip().lower()
        for v in df["sequencing_type"].dropna()
        if v.strip()
    )

    top_df = df[df["similarity_score"] >= min_score].head(top_n)

    lines = [
        f"# Screened Paper Summary — Similar Data Cohort Review",
        f"",
        f"Generated: {date.today().isoformat()}  ",
        f"Source: `literature/00_similar_data_paper_review/`",
        f"",
        f"---",
        f"",
        f"## Overview",
        f"",
        f"| Metric | Count |",
        f"|---|---|",
        f"| Total papers screened | {total} |",
        f"| Similar to scenario | {similar} ({similar/total*100:.0f}%) |",
        f"| Has germline analysis | {has_germline} ({has_germline/total*100:.0f}%) |",
        f"| Has matched tumor-normal | {has_matched} ({has_matched/total*100:.0f}%) |",
        f"| Has both germline + matched normal | {has_both} ({has_both/total*100:.0f}%) |",
        f"",
        f"---",
        f"",
        f"## Similarity Score Distribution",
        f"",
        f"| Score range | Papers |",
        f"|---|---|",
    ]
    for label, count in score_dist.items():
        lines.append(f"| {label} | {count} |")

    lines += [
        f"",
        f"---",
        f"",
        f"## Cancer Type Breakdown",
        f"",
        f"| Cancer type | Papers |",
        f"|---|---|",
    ]
    for ct, n in cancer_counts.most_common(15):
        lines.append(f"| {ct} | {n} |")

    lines += [
        f"",
        f"---",
        f"",
        f"## Sequencing Type Breakdown",
        f"",
        f"| Sequencing type | Papers |",
        f"|---|---|",
    ]
    for st, n in seq_counts.most_common():
        lines.append(f"| {st} | {n} |")

    lines += [
        f"",
        f"---",
        f"",
        f"## Top {len(top_df)} Papers by Similarity Score",
        f"_(score ≥ {min_score})_",
        f"",
    ]

    for _, row in top_df.iterrows():
        score = row["similarity_score"]
        bar = score_bar(score)
        pmid = row.get("pmid", "")
        title = row.get("title", "Untitled")
        journal = row.get("journal", "")
        year = row.get("year", "")
        doi_url = row.get("doi_url", "")
        pubmed_url = row.get("pubmed_url", f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/")
        cancer = row.get("cancer_type", "")
        seq = row.get("sequencing_type", "")
        n = safe_int(row.get("sample_size", ""))
        sample_str = f"n={n}" if n else "n=NR"
        reason = row.get("reason", "")
        findings = fmt_list(row.get("key_findings", ""))
        analyses = fmt_list(row.get("key_analyses", ""))
        endpoints = fmt_list(row.get("clinical_endpoints", ""))
        germ = "✓" if row.get("has_germline", "").lower() == "true" else "✗"
        matched = "✓" if row.get("has_matched_normal", "").lower() == "true" else "✗"
        link = doi_url if doi_url else pubmed_url

        lines += [
            f"### {score:.1f}/10 {bar} — {title}",
            f"",
            f"**{journal}** ({year}) | {cancer} | {seq} | {sample_str} | "
            f"Germline: {germ} | Matched normal: {matched}",
            f"",
            f"[PubMed](https://pubmed.ncbi.nlm.nih.gov/{pmid}/) · [DOI]({link})" if link else f"PMID: {pmid}",
            f"",
            f"**Why similar:** {reason}",
            f"",
            f"**Key findings:** {findings}",
            f"",
            f"**Analyses:** {analyses}",
            f"",
            f"**Clinical endpoints:** {endpoints}",
            f"",
            f"---",
            f"",
        ]

    # --- High impact factor section ---
    if jif_map:
        df["_jif"] = df["journal"].map(lambda j: jif_map.get(j.strip()))
        high_if_df = (
            df[df["_jif"].notna() & (df["_jif"] >= min_if)]
            .sort_values("_jif", ascending=False)
            .head(top_n_if)
        )

        lines += [
            f"---",
            f"",
            f"## Top {len(high_if_df)} Papers by Journal Impact Factor",
            f"_(JIF ≥ {min_if}, verified 2024 Clarivate JCR, max {top_n_if} papers)_",
            f"",
        ]

        for _, row in high_if_df.iterrows():
            jif = row["_jif"]
            sim = row["similarity_score"]
            pmid = row.get("pmid", "")
            title = row.get("title", "Untitled")
            journal = row.get("journal", "")
            year = row.get("year", "")
            doi_url = row.get("doi_url", "")
            cancer = row.get("cancer_type", "")
            seq = row.get("sequencing_type", "")
            n = safe_int(row.get("sample_size", ""))
            sample_str = f"n={n}" if n else "n=NR"
            findings = fmt_list(row.get("key_findings", ""))
            analyses = fmt_list(row.get("key_analyses", ""))
            germ = "✓" if row.get("has_germline", "").lower() == "true" else "✗"
            matched = "✓" if row.get("has_matched_normal", "").lower() == "true" else "✗"
            link = doi_url if doi_url else f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

            lines += [
                f"### JIF {jif:.1f} — {title}",
                f"",
                f"**{journal}** ({year}) | {cancer} | {seq} | {sample_str} | "
                f"Germline: {germ} | Matched normal: {matched} | Similarity: {sim:.1f}/10",
                f"",
                f"[PubMed](https://pubmed.ncbi.nlm.nih.gov/{pmid}/) · [DOI]({link})",
                f"",
                f"**Key findings:** {findings}",
                f"",
                f"**Analyses:** {analyses}",
                f"",
                f"---",
                f"",
            ]
    else:
        lines += [
            f"---",
            f"",
            f"## Top Papers by Journal Impact Factor",
            f"",
            f"_Journal_impact_factor.csv not found — run the journal IF lookup first._",
            f"",
        ]

    return "\n".join(lines)


def main():
    args = parse_args()

    df = load_combined(INPUT_DIR)
    print(f"Loaded {len(df)} papers from {INPUT_DIR}")

    jif_map = load_jif_map(JIF_CSV)
    if jif_map:
        print(f"Loaded {len(jif_map)} journal impact factors from {JIF_CSV}")
    else:
        print(f"Warning: {JIF_CSV} not found — high-IF section will be skipped.")

    md = generate_summary(df, args.min_score, args.top_n, jif_map, args.min_if, args.top_n_if)

    out_path = os.path.join(OUTPUT_DIR, f"summary_{date.today().isoformat()}.md")
    with open(out_path, "w") as f:
        f.write(md)

    print(f"Summary written → {out_path}")


if __name__ == "__main__":
    main()
