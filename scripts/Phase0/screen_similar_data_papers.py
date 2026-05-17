"""
Screen PubMed search results for papers similar to the FDU gastric cancer cohort study.

Reads CSVs from Paper_searching/results/, calls DeepSeek to extract structured data
from each abstract, computes a weighted similarity score, and saves ranked results
to literature/00_similar_data_paper_review/.

Raw DeepSeek responses are saved incrementally to raw_{stem}.json after each batch,
allowing crash recovery and schema-free re-parsing without re-calling the API.

Usage:
    # Normal run — calls DeepSeek, saves raw JSON + screened CSV
    python screen_similar_data_papers.py \
        --csvs cancer_targeted_panel_cohort_clinical_outcomes_2026-05-16.csv \
               cancer_matched_germline_somatic_panel_2026-05-16.csv

    # Re-parse from saved JSON without any API calls
    python screen_similar_data_papers.py \
        --csvs cancer_targeted_panel_cohort_clinical_outcomes_2026-05-16.csv \
        --reparse

Environment:
    DEEPSEEK_API_KEY   required (not needed with --reparse)
"""

import argparse
import json
import os
import sys
import time

import pandas as pd
from openai import OpenAI

DEEPSEEK_BASE_URL = "https://api.deepseek.com"
DEFAULT_MODEL = "deepseek-v4-flash"

PAPER_SEARCH_DIR = "/ShangGaoAIProjects/Paper_searching/results"
OUTPUT_DIR = "/ShangGaoAIProjects/Gastric/genome/literature/00_similar_data_paper_review"

# --- Scenario and prompt --------------------------------------------------------

SCENARIO = """\
A hospital cohort of ~500 gastric cancer patients with targeted panel DNA sequencing.
Each patient has matched tumor and blood samples, providing both somatic mutations
(tumor vs. blood paired calling) and germline variants (blood only). Clinical data
includes MMR status (pMMR/dMMR), HER2 IHC, TMB, MSI status, pathological stage
(pT/pN), overall survival, and IHC markers (claudin18.2, EBER, CDX2, MUC6, KI-67).
A subset of ~45 patients received neoadjuvant chemotherapy ± PD-1 immunotherapy
with tumor regression grade (TRG) scores. Key planned analyses: significantly mutated
gene (SMG) detection, mutational signature decomposition, TMB/MSI biomarker validation,
molecular subtyping (EBV/MSI/GS/CIN), treatment response prediction, survival analysis,
germline pathogenic variant discovery, biallelic tumor suppressor inactivation
(double-hit using matched germline+somatic), pathway-level multi-hit analysis, and
cross-cohort comparison with TCGA-STAD.\
"""

FIELD_DEFINITIONS = """\
For each paper, extract these fields from the title and abstract:

cancer_type       (string)  Cancer type studied, e.g. "gastric", "lung", "pan-cancer", "colorectal".
sequencing_type   (string)  Primary sequencing technology — one of:
                            "targeted panel", "WES", "WGS", "RNA-seq", "multi-omics", "other".
sample_size       (int|null) Number of patients or samples; null if not stated.
has_matched_normal (bool)   True if tumor-normal paired sequencing is described.
has_germline      (bool)    True if germline variants are analyzed.
has_somatic       (bool)    True if somatic mutations are analyzed.
clinical_endpoints (array)  Clinical variables used, e.g. ["overall survival", "TMB", "MSI",
                            "treatment response", "TRG", "HER2"].
key_analyses      (array)   Main analyses performed, e.g. ["biomarker discovery",
                            "survival analysis", "molecular subtyping", "mutation landscape",
                            "treatment response prediction", "mutational signatures"].
key_findings      (array)   2-3 key results from the abstract focused on genomic or clinical
                            findings, e.g. ["TMB-high patients had significantly better OS
                            (HR 0.42, p=0.003)", "Panel-derived TMB concordant with WES at r=0.91"].

score_sequencing  (int 0-3)
    0 = no DNA sequencing (RNA-seq, epigenomics, or proteomics only)
    1 = WGS or WES, no targeted panel
    2 = WES or WGS with strong clinical annotation
    3 = targeted panel sequencing with clinical outcomes data

score_cancer      (int 0-3)
    0 = cancer type unrelated to gastric or GI tract
    1 = any GI cancer (colorectal, esophageal, pancreatic, hepatic, etc.)
    2 = gastric or stomach cancer
    3 = gastric cancer with treatment or immunotherapy context

score_design      (int 0-3)
    0 = no matched or paired sequencing design
    1 = somatic mutations only, no matched normal
    2 = tumor-normal paired sequencing (somatic + matched normal)
    3 = paired sequencing with germline variant analysis AND clinical outcome data

score_analyses    (int 0-3)
    0 = no overlap with the planned analyses listed in the research scenario
    1 = 1 overlapping analysis type
    2 = 2-3 overlapping analysis types
    3 = 4 or more overlapping analysis types

similar_to_scenario (bool)  True if this paper's data structure or analytical questions
                            are meaningfully similar to the research scenario.
reason            (string)  One sentence explaining why similar_to_scenario is true or false.\
"""

SYSTEM_PROMPT = (
    "You are a biomedical literature screening assistant specializing in cancer genomics. "
    "You must return only valid JSON — no prose, no markdown fences."
)


def build_prompt(batch: list[dict]) -> str:
    lines = [
        "Research scenario:",
        SCENARIO,
        "",
        FIELD_DEFINITIONS,
        "",
        'Return a JSON object with a single key "papers" whose value is an array'
        " — one object per paper in the exact schema above, in the same order as the input.",
        "",
        "Papers:",
    ]
    for i, paper in enumerate(batch, 1):
        lines.append(
            f"[{i}] PMID: {paper['pmid']} | Title: {paper['title']} | Abstract: {paper['abstract']}"
        )
    return "\n".join(lines)


# --- Scoring --------------------------------------------------------------------

# Weights: sequencing×2, cancer×2, design×3, analyses×1 → max=24 → normalized to 0–10
WEIGHTS = {"score_sequencing": 2, "score_cancer": 2, "score_design": 3, "score_analyses": 1}
MAX_SCORE = sum(v * 3 for v in WEIGHTS.values())  # 24


def compute_similarity_score(extracted: dict) -> float:
    try:
        raw = sum(
            int(extracted.get(k, 0) or 0) * w for k, w in WEIGHTS.items()
        )
        return round(raw / MAX_SCORE * 10, 2)
    except (TypeError, ValueError):
        return 0.0


# --- API calls ------------------------------------------------------------------

def call_api(client: OpenAI, model: str, batch: list[dict]) -> str:
    response = client.chat.completions.create(
        model=model,
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": build_prompt(batch)},
        ],
        response_format={"type": "json_object"},
    )
    return response.choices[0].message.content


def parse_response(raw: str, batch: list[dict]) -> dict:
    """Map response items to pmids by position — DeepSeek does not echo pmid in output."""
    data = json.loads(raw)
    if not isinstance(data, dict) or "papers" not in data:
        raise ValueError(f"Expected JSON object with 'papers' key, got: {type(data).__name__} "
                         f"with keys {list(data.keys()) if isinstance(data, dict) else 'N/A'}")
    items = data["papers"]
    if not isinstance(items, list):
        raise ValueError(f"'papers' value must be a list, got {type(items).__name__}")
    if len(items) != len(batch):
        raise ValueError(f"Response length {len(items)} != batch length {len(batch)}")
    return {str(paper["pmid"]).strip(): item for paper, item in zip(batch, items)}


def process_batch(client: OpenAI, model: str, batch: list[dict], max_retries: int = 3) -> dict:
    last_exc = None
    for attempt in range(max_retries):
        try:
            raw = call_api(client, model, batch)
            try:
                return parse_response(raw, batch)
            except (json.JSONDecodeError, ValueError):
                raw2 = call_api(client, model, batch)
                return parse_response(raw2, batch)
        except (json.JSONDecodeError, ValueError) as exc:
            last_exc = exc
            print(f"  Parse error (attempt {attempt + 1}): {exc}")
            break
        except Exception as exc:
            last_exc = exc
            wait = 2 ** attempt
            print(f"  API error (attempt {attempt + 1}): {exc}. Retrying in {wait}s...")
            time.sleep(wait)
    print(f"  Batch failed: {last_exc}")
    return {}


# --- Raw JSON persistence -------------------------------------------------------

def load_raw_json(path: str) -> dict:
    if os.path.isfile(path):
        with open(path) as f:
            return json.load(f)
    return {}


def save_raw_json(path: str, results: dict) -> None:
    with open(path, "w") as f:
        json.dump(results, f, indent=2)


# --- Main screening logic -------------------------------------------------------

EXTRACTED_COLS = [
    "cancer_type", "sequencing_type", "sample_size",
    "has_matched_normal", "has_germline", "has_somatic",
    "clinical_endpoints", "key_analyses", "key_findings",
    "score_sequencing", "score_cancer", "score_design", "score_analyses",
    "similar_to_scenario", "reason",
]
LIST_COLS = {"clinical_endpoints", "key_analyses", "key_findings"}


def screen_csv(
    client: OpenAI,
    model: str,
    csv_path: str,
    batch_size: int,
    raw_json_path: str,
    reparse: bool,
) -> pd.DataFrame:
    df = pd.read_csv(csv_path, dtype=str).fillna("")
    has_abstract = df["abstract"].str.strip() != ""
    df_with = df[has_abstract].copy()
    df_without = df[~has_abstract].copy()

    print(f"\n{os.path.basename(csv_path)}: {len(df)} papers "
          f"({len(df_with)} with abstracts, {len(df_without)} skipped)")

    # Load any previously saved results (enables crash recovery and --reparse)
    results: dict[str, dict] = load_raw_json(raw_json_path)
    if results:
        print(f"  Loaded {len(results)} existing results from {os.path.basename(raw_json_path)}")

    if reparse:
        print("  --reparse: skipping API calls, using saved JSON only.")
    else:
        records = df_with[["pmid", "title", "abstract"]].to_dict(orient="records")
        pending = [r for r in records if str(r["pmid"]).strip() not in results]

        if not pending:
            print("  All papers already processed — nothing to call.")
        else:
            print(f"  {len(pending)} papers to process ({len(records) - len(pending)} already cached).")
            batches = [pending[i: i + batch_size] for i in range(0, len(pending), batch_size)]
            for idx, batch in enumerate(batches):
                start = idx * batch_size + 1
                end = min(start + batch_size - 1, len(pending))
                print(f"  Batch {idx + 1}/{len(batches)} — papers {start}–{end} ...")
                batch_results = process_batch(client, model, batch)
                for paper in batch:
                    pmid = str(paper["pmid"]).strip()
                    results[pmid] = batch_results.get(pmid, {"reason": "parse_error"})
                # Save after every batch so a crash loses at most one batch
                save_raw_json(raw_json_path, results)
            print(f"  Raw JSON saved → {raw_json_path}")

    for col in EXTRACTED_COLS:
        df_with[col] = df_with["pmid"].map(
            lambda p, c=col: results.get(str(p).strip(), {}).get(c)
        )

    for col in LIST_COLS:
        df_with[col] = df_with[col].apply(
            lambda v: "; ".join(v) if isinstance(v, list) else (v or "")
        )

    df_with["similarity_score"] = df_with["pmid"].apply(
        lambda p: compute_similarity_score(results.get(str(p).strip(), {}))
    )

    for col in EXTRACTED_COLS:
        df_without[col] = None
    df_without["similarity_score"] = 0.0

    return pd.concat([df_with, df_without], ignore_index=True)


# --- CLI ------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Screen papers for similarity to FDU gastric cancer cohort study."
    )
    parser.add_argument(
        "--csvs", nargs="+", required=True,
        help="CSV filenames inside Paper_searching/results/",
    )
    parser.add_argument("--model", default=DEFAULT_MODEL, help=f"DeepSeek model (default: {DEFAULT_MODEL})")
    parser.add_argument("--batch-size", type=int, default=10, help="Papers per API call (default: 10)")
    parser.add_argument("--reparse", action="store_true",
                        help="Re-parse from saved raw JSON without calling the API")
    return parser.parse_args()


def main():
    args = parse_args()

    api_key = os.environ.get("DEEPSEEK_API_KEY", "")
    if not api_key and not args.reparse:
        sys.exit("Error: DEEPSEEK_API_KEY environment variable is not set.")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    client = OpenAI(api_key=api_key or "none", base_url=DEEPSEEK_BASE_URL)

    all_dfs = []

    for csv_name in args.csvs:
        csv_path = os.path.join(PAPER_SEARCH_DIR, csv_name)
        if not os.path.isfile(csv_path):
            print(f"Warning: {csv_path} not found, skipping.")
            continue

        stem = os.path.splitext(csv_name)[0]
        raw_json_path = os.path.join(OUTPUT_DIR, f"raw_{stem}.json")

        df = screen_csv(client, args.model, csv_path, args.batch_size, raw_json_path, args.reparse)

        out_path = os.path.join(OUTPUT_DIR, f"{stem}_screened.csv")
        df.to_csv(out_path, index=False)
        print(f"  Screened CSV → {out_path}")

        all_dfs.append(df)

    if len(all_dfs) > 1:
        combined = pd.concat(all_dfs, ignore_index=True)
        combined = combined.drop_duplicates(subset=["pmid"])
        combined = combined.sort_values("similarity_score", ascending=False).reset_index(drop=True)
        combined_path = os.path.join(OUTPUT_DIR, "combined_ranked.csv")
        combined.to_csv(combined_path, index=False)
        print(f"\nCombined ranked → {combined_path} ({len(combined)} unique papers)")


if __name__ == "__main__":
    main()
