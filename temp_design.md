# Multimodal RAG — POC Plan

### Goal

Given a genomic finding from our own FDU cohort, automatically build an evidence chain across external omics data and structured biological knowledge, then synthesize a coherent narrative. Literature is a secondary layer — only consulted when the data story is strong and needs mechanistic naming.

**End-to-end: Finding → Evidence Plan → Multi-omics retrieval → Narrative → (optional) Literature**

---

## Core Architecture

```
FDU somatic finding (seed)
        ↓
[LLM] Evidence Plan
        ↓
Outside data — evidence chain (primary)
  ├── TCGA-STAD somatic   → does the finding replicate?
  ├── TCGA-STAD RNA-seq   → is expression altered when mutated?
  ├── TCGA-STAD CNV       → is the gene also amplified/deleted?
  ├── TCGA-STAD methylation → is it epigenetically dysregulated?
  ├── TCGA-STAD clinical  → survival / treatment response correlation
  └── Neo4j KB            → pathway membership, GO terms, BioGRID interactions
        ↓
[LLM] Is the evidence chain strong?
        ↓ yes                         ↓ no
[Literature layer]               Stop / revise finding
  ├── Generate PubMed query from Neo4j context + finding
  ├── Fetch top abstracts (PubMed E-utilities API)
  ├── LLM: extract mechanism from abstracts
  └── If promising → full article card
        ↓
[LLM] Narrative synthesis
```

---

## Phase 0 — Foundation

**Step 0.1 — Generate Finding Cards from FDU somatic MAF**
Run a minimal analysis: mutation frequency per gene, top mutated genes, variant classification breakdown, pre-treatment vs post-treatment differences. Serialize results into structured Finding Cards. These are the seeds for the whole pipeline.

**Step 0.2 — Audit available TCGA-STAD data layers**
Inventory what is actually usable in each omics folder (`somatic_mutation`, `rnaseq`, `copy_number`, `methylation`, `clinical`). Map sample overlap across layers — the evidence chain is only as strong as the sample intersection.

**Step 0.3 — Audit Neo4j KB**
For a representative set of top FDU genes, query Neo4j and map what is returned: pathways, GO terms, interactors, disease associations. This defines what the Evidence Plan can ask Neo4j for.

**Step 0.4 — Design the Evidence Plan schema**
Based on 0.1–0.3, define what an Evidence Plan looks like: which layers to query, what question each layer answers, what a "strong" result looks like per layer.

---

## Phase 1 — Multi-omics Retrieval

**Step 1.1 — TCGA somatic replication**
For each top FDU gene/finding, query TCGA-STAD somatic mutations. Report frequency, variant spectrum, co-mutation patterns.

**Step 1.2 — Expression consequence**
In TCGA RNA-seq, compare expression of the mutated gene (and pathway members) between mutant vs wildtype samples.

**Step 1.3 — CNV and methylation context**
Check if the gene shows copy number loss or promoter methylation in the same samples. Multi-hit evidence strengthens the story.

**Step 1.4 — Clinical correlation**
Correlate mutation status with survival (OS/DFS) and treatment response in TCGA clinical data.

**Step 1.5 — Neo4j biological context**
Query Neo4j for the gene's pathway memberships, GO terms, and top interactors. This frames the biological interpretation of the data findings.

---

## Phase 2 — Evidence Assessment and Literature (conditional)

**Step 2.1 — LLM evidence assessment**
Given all retrieved omics + Neo4j evidence, the LLM evaluates whether the chain is coherent and complete. If yes, proceed to narrative. If a mechanistic gap remains, trigger literature.

**Step 2.2 — Literature (on demand)**
Using the Neo4j context + finding, generate a targeted PubMed search string. Fetch top abstracts. LLM extracts mechanism. Only if a finding is novel or unexpected, escalate to a full article card.

**Step 2.3 — Narrative synthesis**
LLM synthesizes all evidence layers into a coherent biological story grounded in the FDU data, validated by TCGA, contextualized by Neo4j, and (if triggered) named by literature.

---

## Phase 3 — End-to-End Test

Run the full pipeline on 2–3 Finding Cards. Evaluate manually:
- Does the Evidence Plan correctly identify which omics layers are relevant?
- Does each layer return meaningful signal?
- Is the evidence chain coherent without literature?
- Does the narrative accurately reflect the data?

---

## Open Decisions

| # | Decision | Options |
|---|---|---|
| 1 | **Evidence Plan format** | Fixed template per layer vs LLM-generated free-form plan |
| 2 | **"Strong evidence" threshold** | How many layers must show signal before proceeding to narrative? |
| 3 | **TCGA sample overlap** | Use only samples present in all layers vs layer-specific sample sets |
| 4 | **Literature trigger** | Always run literature vs only when LLM flags a mechanistic gap |
| 5 | **What does "done" look like?** | A working notebook? A CLI script? An API? |
