# Project Rules

## Analysis Skills — MANDATORY

The following analyses have standardized skills. When the user asks for these analyses, ALWAYS invoke the corresponding skill. Never write new code or reimplement the analysis from scratch.

| Analysis | Skill to invoke | Script |
|----------|----------------|--------|
| DepMap CRISPR dependency / synthetic lethality | `/depmap-synleth` | `scripts/depmap/depmap_synleth.py` |

## Before writing any new script

Before writing or modifying any analysis script, list all ambiguous parameter choices and methodological decisions, and wait for user confirmation. Do not write a single line of code until the user has confirmed every decision point.
