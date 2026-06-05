# DepMap Synthetic Lethality Analysis

**Date:** 2026-06-04 18:46

## Parameters

| Parameter | Value |
|-----------|-------|
| genes | ATM |
| mut_gene | TP53 |
| cancer | Esophagus/Stomach |
| pancancer | False |
| n_mut | 76 |
| n_wt | 14 |

## Analysis Decisions

| Decision | Value |
|----------|-------|
| Cancer filter field | OncotreeLineage |
| Mutant definition | damaging mutation count > 0 OR CN < 0.3 |
| Min sample size | mut ≥ 5, WT ≥ 3 |
| Test (n ≥ 20 both) | Mann-Whitney U (two-sided) |
| Test (n < 20 either) | Permutation test (10,000 iterations) |
| CI | Bootstrap 95%, 2,000 iterations, delta = mean(mut) − mean(WT) |

## Results

| Gene | N_mut | N_WT | Median_mut | Median_WT | Delta | CI 95% | p | Test | Sig |
|------|-------|------|-----------|----------|-------|--------|---|------|-----|
| ATM | 60 | 9 | -0.0716 | 0.0538 | -0.102 | [-0.205, 0.005] | 0.093 | permutation | ns |

## Figures

### ATM
- **stripbox**: `ATM_stripbox.png`
- **violin**: `ATM_violin.png`
