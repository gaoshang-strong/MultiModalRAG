# Project Rules

## Analysis Skills — MANDATORY

The following analyses have standardized skills. When the user asks for these analyses, ALWAYS invoke the corresponding skill. Never write new code or reimplement the analysis from scratch.

| Analysis | Skill to invoke | Script |
|----------|----------------|--------|
| DepMap CRISPR dependency / synthetic lethality | `/depmap-synleth` | `scripts/depmap/depmap_synleth.py` |
| mRNA differential expression (TCGA-STAD, TP53-mut vs WT) | `/mrna-expression` | `scripts/mrna/mrna_expression.R` |
| Neo4j pathway query + DeepSeek interpretation | `/neo4j-query` | (Cypher, no fixed script) |
| Pathway network figure from JSON config | (called within `/neo4j-query` step 7) | `scripts/neo4j/draw_pathway.py` |

## Before writing any new script

Before writing or modifying any analysis script, list all ambiguous parameter choices and methodological decisions, and wait for user confirmation. Do not write a single line of code until the user has confirmed every decision point.
