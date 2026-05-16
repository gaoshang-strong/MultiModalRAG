# Biomedical Multimodal RAG — Gastric Cancer Multi-Omics

## Goal

Build a multimodal RAG system over TCGA gastric cancer data (RNAseq + somatic mutations) as a vehicle for developing expertise as a biomedical AI agent scientist.

---

## Core Innovation

**The problem with multi-omics integration:** When combining data from different labs, cohorts, or platforms, batch effects and patient population differences make direct data-level integration statistically fragile — sometimes impossible.

**The insight:** Don't integrate the data. Integrate the findings.

Analyze each omics dataset independently using rigorous, well-established methods. Then use your own findings as seeds to query a knowledge base built from external omics, literature, and public databases. The RAG system becomes the cross-omics reasoning layer — it identifies convergent evidence, surfaces contradictions, and assembles a coherent biological story.

This mirrors how good scientists actually think: start from what your own data tells you, then ask whether the world agrees.

```
Own RNAseq findings ──┐
                       ├──► Seed queries ──► Knowledge Base ──► Story
Own mutation findings ─┘                     (lit + public omics)
```

The output is not a merged matrix. It is a structured narrative: gene X (mutated in your cohort) → pathway Y (dysregulated in your RNAseq) → supported by N independent studies → clinically annotated in COSMIC.

---

## Skills This Project Develops

### For a Senior AI Agent Scientist role at a biotech company:

**AI/ML Engineering**
- RAG architecture: chunking, embedding, hybrid retrieval, reranking
- Agent design: tool use, multi-step reasoning, LangGraph/LangChain
- Vector databases and knowledge graph construction
- Structured LLM output and schema design
- RAG evaluation: faithfulness, hallucination detection, domain QA

**Biomedical Domain**
- Multi-omics data types and their analytical idioms (DEG, MAF, GSEA)
- Pathway and ontology databases: KEGG, Reactome, MSigDB, GO
- Clinical/genomic databases: TCGA, COSMIC, ClinVar, GEO
- Literature retrieval and semantic search over biomedical papers

**Systems**
- ML workflow orchestration
- Model serving and structured output parsing
- Evaluation frameworks for domain-specific RAG

---

## Data Sources

- **Own data:** TCGA-STAD RNAseq and somatic mutation (MAF) files
- **Knowledge base:** PubMed literature, COSMIC, MSigDB gene sets, pathway databases, gene information cards

---

## High-Level Architecture

```
[TCGA RNAseq]  [TCGA Somatic Mut.]
      │                │
      ▼                ▼
 Independent       Independent
  Analysis          Analysis
      │                │
      └────── Finding Cards ──────┐
                                  ▼
                         [Knowledge Base]
                    literature + public omics
                       + pathway databases
                                  │
                                  ▼
                        Story-Building Agent
                     (seed → retrieve → narrate)
                                  │
                                  ▼
                    Structured Biological Narrative
                    with citations and evidence score
```

---

## Design Principles

- Analyze each dataset with its own appropriate statistical method
- Never force data harmonization across cohorts
- Use the RAG layer for reasoning, not for data fusion
- Every claim in the output narrative must be traceable to a source
- Evaluate the system on biological coherence, not just retrieval metrics
