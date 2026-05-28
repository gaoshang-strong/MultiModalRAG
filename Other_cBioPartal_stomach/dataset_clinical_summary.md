# cBioPortal Stomach Cancer Datasets — Clinical Summary

## Overview

| Dataset | N Patients | Cancer Type | Sequencing | Source |
|---|---|---|---|---|
| egc_tmucih_2015 | 78 | Esophagogastric Adenocarcinoma | WXS | TMUCIH (China) |
| stad_oncosg_2018 | ~146 | Stomach Adenocarcinoma | WXS | Multi-cohort: ICGC-SG, HK, SG |
| stad_pfizer_uhongkong | 100 | Stomach Adenocarcinoma | WXS | Hong Kong (Pfizer-funded) |
| stad_uhongkong | 22 | Stomach Adenocarcinoma | WXS | Hong Kong (subset with survival) |
| stad_utokyo | 30 | Stomach Adenocarcinoma | WXS | Tokyo (Japan), diffuse-only |

---

## 1. egc_tmucih_2015

**N = 78 patients** | Esophagogastric Adenocarcinoma | WXS, matched tumor-normal

### Patient-level variables
| Variable | Description |
|---|---|
| AGE | Age at diagnosis (range: ~25–82) |
| CLIN_T/N/M_STAGE | Clinical TNM staging |
| EBV_STATUS | Positive / Negative |
| EBV_PRESENT | EBV copy number (numeric) |
| LOCATION | Body / Antrum / Cardia |
| OS_MONTHS | Overall survival in months |
| OS_STATUS | 0:LIVING / 1:DECEASED |
| SEX | Male / Female |
| STAGE | UICC stage (1–4) |
| SUBTYPE | I (intestinal) / D (diffuse) |

### Sample-level variables
| Variable | Description |
|---|---|
| ARID1A_MUTATION | Wt / Mut |
| CDH1_STATUS | CDH1_wt / CDH1_mut |
| CLONAL_LABEL | Ultraclonal / Non-Ultraclonal |
| DNA_REPAIR_STATUS | Def (deficient) / NonDef |
| GRADE | Tumor grade (1–3, partially missing) |
| MUTATION_COUNT | Somatic mutation count |
| TMB_NONSYNONYMOUS | Nonsynonymous TMB (mut/Mb) |
| TP53_MUTATION | Wt / Mut |
| PLATFORM | WXS |
| SOMATIC_STATUS | Matched |

### Notes
- Oncotree code: EGC (Esophagogastric Cancer)
- Survival data (OS) available for most patients
- EBV copy number provided alongside binary EBV_STATUS
- About half the cohort is EBV-positive

---

## 2. stad_oncosg_2018

**N ≈ 146 patients** | Stomach Adenocarcinoma | WXS, matched tumor-normal  
Multi-cohort: ICGC-Singapore (~24), Hong Kong Pfizer (~91), Singapore Apollo/Tan (~31)

### Patient-level variables
| Variable | Description |
|---|---|
| AGE | Age at diagnosis (HK and SG cohorts; missing for some ICGC) |
| COHORT | ICGC / HK / SG |
| LAURENS_CLASSIFICATION | Intestinal / Diffuse / Mixed |
| MOLECULAR_SUBTYPE | EBV / MSI / CIN / GS |
| STAGE / STAGE_T/N/M | Pathologic TNM staging |
| TUMOR_SITE | Antrum / Body / Cardia / Lesser Curve / Greater Curve / Pylorus / Incisura |
| SEX | Male / Female |

### Sample-level variables
| Variable | Description |
|---|---|
| MUTATION_COUNT | Somatic mutation count |
| TMB_NONSYNONYMOUS | Nonsynonymous TMB (mut/Mb) |
| ONCOTREE_CODE | STAD |
| SOMATIC_STATUS | Matched |
| XCELL_* | xCell immune cell deconvolution scores (partial; available for SG and some HK samples) |

### Notes
- Oncotree code: STAD
- TCGA-like molecular subtypes (EBV, MSI, CIN, GS) available for HK and SG subcohorts; missing for ICGC
- Lauren's classification available for HK and SG subcohorts; missing for ICGC
- xCell immune microenvironment deconvolution scores available for a subset (~30 samples)
- Tumor site uses mixed granularity (some: Antrum/Body/Cardia; others: Lesser Curve/Greater Curve etc.)

---

## 3. stad_pfizer_uhongkong

**N = 100 patients** | Stomach Adenocarcinoma | WXS, matched tumor-normal  
Hong Kong cohort, Pfizer-funded

### Patient-level variables
| Variable | Description |
|---|---|
| AGE | Age at diagnosis (range: ~32–89) |
| EBV_STATUS | 0 (negative) / 1 (positive) |
| H_PYLORI_INFECTION | 0 / 1 |
| M_STAGE | Distant metastasis (0/1) |
| N_STAGE | Nodal stage (0–3) |
| NEOADJUVANT_CHEMO | 0 / 1 |
| PATHOLOGY_LAUREN | intestinal / diffuse / mixed |
| TUMOR_SITE | antrum / body / cardia / diffuse/site unknown |
| UICC_TUMOR_STAGE | IA / IB / II / IIIA / IIIB / IV |
| SEX | Male / Female |

### Sample-level variables
| Variable | Description |
|---|---|
| MSI_STATUS | microsatellite stable or low-level MSI / high-level MSI |
| MUTATION_COUNT | Somatic mutation count |
| TMB_NONSYNONYMOUS | Nonsynonymous TMB (mut/Mb) |
| TUMOR_CONTENT | Estimated tumor cellularity (%) |
| TUMOR_DIFFERENTIATION | well / moderate / poor |
| T_STAGE | T1–T4 |
| ONCOTREE_CODE | STAD |
| SOMATIC_STATUS | Matched |

### Notes
- **No survival data** (no OS or DFS columns)
- H. pylori infection status and neoadjuvant chemotherapy recorded
- MSI status phrased as full text (not MSS/MSI-H shorthand)
- Patient IDs overlap with stad_uhongkong (same pfg numbering); this is the larger, Pfizer-sponsored version without survival data

---

## 4. stad_uhongkong

**N = 22 patients** | Stomach Adenocarcinoma | WXS, matched tumor-normal  
Hong Kong cohort (earlier / smaller subset with survival outcomes)

### Patient-level variables
| Variable | Description |
|---|---|
| AGE | Age at diagnosis (range: ~37–87) |
| DFS_MONTHS | Disease-free survival in months |
| DFS_STATUS | 0:DiseaseFree / 1:Recurred/Progressed |
| EBV_STATUS | Negative / Positive |
| H_PYLORI_INFECTION | 0 / 1 |
| M_STAGE | Distant metastasis (0/1) |
| N_STAGE | Nodal stage (0–3) |
| OS_MONTHS | Overall survival in months |
| OS_STATUS | 0:LIVING / 1:DECEASED |
| PATHOLOGY_LAUREN | intestinal / diffuse (signet ring cell) / mixed |
| TUMOR_SITE | antrum / body / cardia |
| UICC_TUMOR_STAGE | 0 / IA / IB / II / IIIA / IIIB / IV |
| VITAL_STATUS | Alive or censored / Died of disease |
| SEX | Male / Female |

### Sample-level variables
| Variable | Description |
|---|---|
| MSI_STATUS | MSS / MSI |
| MUTATION_COUNT | Somatic mutation count |
| TMB_NONSYNONYMOUS | Nonsynonymous TMB (mut/Mb) |
| TUMOR_CONTENT | Estimated tumor cellularity (%) |
| TUMOR_DIFFERENTIATION | well / moderate / poor |
| TYPE_OF_SURGERY | curative / palliative |
| T_STAGE | Tis / T1–T4 |
| ONCOTREE_CODE | STAD |
| SOMATIC_STATUS | Matched |

### Notes
- **Both OS and DFS** available — the only HK dataset with disease-free survival
- Includes type of surgery (curative vs palliative)
- Lauren's classification distinguishes "diffuse (signet ring cell)" explicitly
- Patient IDs overlap with stad_pfizer_uhongkong; this is a smaller, earlier subset
- One patient (pfg282T) has stage 0 (Tis); one (pfg016T) is stage IV but alive (long survivor)

---

## 5. stad_utokyo

**N = 30 patients** | Stomach Adenocarcinoma | WXS, matched tumor-normal  
University of Tokyo, Japan; diffuse-type gastric cancer focus

### Patient-level variables
| Variable | Description |
|---|---|
| AGE_GROUP | Age range (e.g., 40–80, 60–70, 70–80, 80–90); ~half missing |
| BONE_SCAN_RESULT | Integer score (0–5); ~half missing |
| LAUREN_CLASS | All **Diffuse** |
| M_STAGE | M0 / M1 |
| N_STAGE | N0 / N1 / N2 / N3a / N3b |
| TUMOR_SITE | U (upper) / M (middle) / L (lower) |
| UICC_TUMOR_STAGE | IB / IIA / IIB / IIIA / IIIB / IIIC / IV; ~half missing |
| SEX | Male / Female |

### Sample-level variables
| Variable | Description |
|---|---|
| GRADE | Poor / Poor (signet); partial |
| HER2_IHC_SCORE | 0 / 1; partial |
| MUTATION_COUNT | Somatic mutation count |
| TMB_NONSYNONYMOUS | Nonsynonymous TMB (mut/Mb) |
| RHOA_MUTATION | Specific RHOA variant (Y42C, L22R, R5W, Y74D, R5Q); partial |
| ONCOTREE_CODE | STAD |
| SOMATIC_STATUS | Matched |
| T_STAGE | T1b / T2 / T3 / T4a; partial |

### Notes
- **All diffuse-type gastric cancer** — dataset specifically designed to study diffuse/signet-ring histology
- **RHOA mutations** explicitly annotated; known driver in diffuse gastric cancer
- **No survival data** (no OS or DFS)
- Age given as grouped ranges, not exact values; ~half of patients have missing clinical annotations
- Bone scan result included (unusual variable; possibly for bone metastasis assessment)
- HER2 IHC score available for partially annotated cases
- Tumor site uses anatomical thirds (U/M/L) rather than named regions

---

## Cross-Dataset Comparison

| Feature | egc_tmucih_2015 | stad_oncosg_2018 | stad_pfizer_uhongkong | stad_uhongkong | stad_utokyo |
|---|:---:|:---:|:---:|:---:|:---:|
| N patients | 78 | ~146 | 100 | 22 | 30 |
| OS data | Yes | No | No | Yes | No |
| DFS data | No | No | No | Yes | No |
| Lauren's classification | Indirect (Subtype I/D) | Yes | Yes | Yes | Yes (all Diffuse) |
| Molecular subtype | Partial (EBV, DNA repair) | Yes (EBV/MSI/CIN/GS) | MSI status | MSI status | No |
| H. pylori | No | No | Yes | Yes | No |
| EBV status | Yes | Partial | Yes | Yes | No |
| TMB | Yes | Yes | Yes | Yes | Yes |
| RHOA mutation | No | No | No | No | Yes |
| HER2 | No | No | No | No | Partial |
| Immune deconvolution | No | Partial (xCell) | No | No | No |
| Surgery type | No | No | No | Yes | No |
| Population | Chinese (Beijing) | Multi (SG/HK) | Chinese (HK) | Chinese (HK) | Japanese |
