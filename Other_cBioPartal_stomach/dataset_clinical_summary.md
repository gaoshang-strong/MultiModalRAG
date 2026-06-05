# cBioPortal Stomach Cancer Datasets — Clinical Summary

## Overview


| Dataset               | N Patients        | Cancer Type                    | Sequencing    | Source                           |
| --------------------- | ----------------- | ------------------------------ | ------------- | -------------------------------- |
| egc_tmucih_2015       | 78                | Esophagogastric Adenocarcinoma | WXS           | TMUCIH (China)                   |
| stad_oncosg_2018      | ~146              | Stomach Adenocarcinoma         | WXS           | Multi-cohort: ICGC-SG, HK, SG    |
| stad_pfizer_uhongkong | 100               | Stomach Adenocarcinoma         | WXS           | Hong Kong (Pfizer-funded)        |
| stad_uhongkong        | 22                | Stomach Adenocarcinoma         | WXS           | Hong Kong (subset with survival) |
| stad_utokyo           | 30                | Stomach Adenocarcinoma         | WXS           | Tokyo (Japan), diffuse-only      |
| egc_msk_2017          | 305 (341 samples) | EGC (STAD + GEJ + ESCA)        | MSK-IMPACT341 | MSK, metastatic                  |
| egc_mskcc_2020        | 487               | EGC (STAD + GEJ + ESCA)        | MSK-IMPACT    | MSK                              |
| egc_msk_tp53_ccr_2022 | 237               | EGC (STAD + GEJ + ESCA)        | MSK-IMPACT    | MSK, TP53-focused                |
| egc_msk_2023          | 902               | EGC (STAD + GEJ + ESCA)        | MSK-IMPACT341 | MSK, early vs average onset      |


---

## 1. egc_tmucih_2015

**N = 78 patients** | Esophagogastric Adenocarcinoma | WXS, matched tumor-normal

### Patient-level variables


| Variable         | Description                      |
| ---------------- | -------------------------------- |
| AGE              | Age at diagnosis (range: ~25–82) |
| CLIN_T/N/M_STAGE | Clinical TNM staging             |
| EBV_STATUS       | Positive / Negative              |
| EBV_PRESENT      | EBV copy number (numeric)        |
| LOCATION         | Body / Antrum / Cardia           |
| OS_MONTHS        | Overall survival in months       |
| OS_STATUS        | 0:LIVING / 1:DECEASED            |
| SEX              | Male / Female                    |
| STAGE            | UICC stage (1–4)                 |
| SUBTYPE          | I (intestinal) / D (diffuse)     |


### Sample-level variables


| Variable          | Description                          |
| ----------------- | ------------------------------------ |
| ARID1A_MUTATION   | Wt / Mut                             |
| CDH1_STATUS       | CDH1_wt / CDH1_mut                   |
| CLONAL_LABEL      | Ultraclonal / Non-Ultraclonal        |
| DNA_REPAIR_STATUS | Def (deficient) / NonDef             |
| GRADE             | Tumor grade (1–3, partially missing) |
| MUTATION_COUNT    | Somatic mutation count               |
| TMB_NONSYNONYMOUS | Nonsynonymous TMB (mut/Mb)           |
| TP53_MUTATION     | Wt / Mut                             |
| PLATFORM          | WXS                                  |
| SOMATIC_STATUS    | Matched                              |


### Notes

- Oncotree code: EGC (Esophagogastric Cancer)
- Survival data (OS) available for most patients
- EBV copy number provided alongside binary EBV_STATUS
- About half the cohort is EBV-positive

---

## 2. stad_oncosg_2018

**N ≈ 146 patients** | Stomach Adenocarcinoma | WXS, matched tumor-normal  
Multi-cohort: ICGC-Singapore (~~24), Hong Kong Pfizer (~~91), Singapore Apollo/Tan (~31)

### Patient-level variables


| Variable               | Description                                                                |
| ---------------------- | -------------------------------------------------------------------------- |
| AGE                    | Age at diagnosis (HK and SG cohorts; missing for some ICGC)                |
| COHORT                 | ICGC / HK / SG                                                             |
| LAURENS_CLASSIFICATION | Intestinal / Diffuse / Mixed                                               |
| MOLECULAR_SUBTYPE      | EBV / MSI / CIN / GS                                                       |
| STAGE / STAGE_T/N/M    | Pathologic TNM staging                                                     |
| TUMOR_SITE             | Antrum / Body / Cardia / Lesser Curve / Greater Curve / Pylorus / Incisura |
| SEX                    | Male / Female                                                              |


### Sample-level variables


| Variable          | Description                                                                            |
| ----------------- | -------------------------------------------------------------------------------------- |
| MUTATION_COUNT    | Somatic mutation count                                                                 |
| TMB_NONSYNONYMOUS | Nonsynonymous TMB (mut/Mb)                                                             |
| ONCOTREE_CODE     | STAD                                                                                   |
| SOMATIC_STATUS    | Matched                                                                                |
| XCELL_*           | xCell immune cell deconvolution scores (partial; available for SG and some HK samples) |


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


| Variable           | Description                                   |
| ------------------ | --------------------------------------------- |
| AGE                | Age at diagnosis (range: ~32–89)              |
| EBV_STATUS         | 0 (negative) / 1 (positive)                   |
| H_PYLORI_INFECTION | 0 / 1                                         |
| M_STAGE            | Distant metastasis (0/1)                      |
| N_STAGE            | Nodal stage (0–3)                             |
| NEOADJUVANT_CHEMO  | 0 / 1                                         |
| PATHOLOGY_LAUREN   | intestinal / diffuse / mixed                  |
| TUMOR_SITE         | antrum / body / cardia / diffuse/site unknown |
| UICC_TUMOR_STAGE   | IA / IB / II / IIIA / IIIB / IV               |
| SEX                | Male / Female                                 |


### Sample-level variables


| Variable              | Description                                             |
| --------------------- | ------------------------------------------------------- |
| MSI_STATUS            | microsatellite stable or low-level MSI / high-level MSI |
| MUTATION_COUNT        | Somatic mutation count                                  |
| TMB_NONSYNONYMOUS     | Nonsynonymous TMB (mut/Mb)                              |
| TUMOR_CONTENT         | Estimated tumor cellularity (%)                         |
| TUMOR_DIFFERENTIATION | well / moderate / poor                                  |
| T_STAGE               | T1–T4                                                   |
| ONCOTREE_CODE         | STAD                                                    |
| SOMATIC_STATUS        | Matched                                                 |


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


| Variable           | Description                                     |
| ------------------ | ----------------------------------------------- |
| AGE                | Age at diagnosis (range: ~37–87)                |
| DFS_MONTHS         | Disease-free survival in months                 |
| DFS_STATUS         | 0:DiseaseFree / 1:Recurred/Progressed           |
| EBV_STATUS         | Negative / Positive                             |
| H_PYLORI_INFECTION | 0 / 1                                           |
| M_STAGE            | Distant metastasis (0/1)                        |
| N_STAGE            | Nodal stage (0–3)                               |
| OS_MONTHS          | Overall survival in months                      |
| OS_STATUS          | 0:LIVING / 1:DECEASED                           |
| PATHOLOGY_LAUREN   | intestinal / diffuse (signet ring cell) / mixed |
| TUMOR_SITE         | antrum / body / cardia                          |
| UICC_TUMOR_STAGE   | 0 / IA / IB / II / IIIA / IIIB / IV             |
| VITAL_STATUS       | Alive or censored / Died of disease             |
| SEX                | Male / Female                                   |


### Sample-level variables


| Variable              | Description                     |
| --------------------- | ------------------------------- |
| MSI_STATUS            | MSS / MSI                       |
| MUTATION_COUNT        | Somatic mutation count          |
| TMB_NONSYNONYMOUS     | Nonsynonymous TMB (mut/Mb)      |
| TUMOR_CONTENT         | Estimated tumor cellularity (%) |
| TUMOR_DIFFERENTIATION | well / moderate / poor          |
| TYPE_OF_SURGERY       | curative / palliative           |
| T_STAGE               | Tis / T1–T4                     |
| ONCOTREE_CODE         | STAD                            |
| SOMATIC_STATUS        | Matched                         |


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


| Variable         | Description                                                 |
| ---------------- | ----------------------------------------------------------- |
| AGE_GROUP        | Age range (e.g., 40–80, 60–70, 70–80, 80–90); ~half missing |
| BONE_SCAN_RESULT | Integer score (0–5); ~half missing                          |
| LAUREN_CLASS     | All **Diffuse**                                             |
| M_STAGE          | M0 / M1                                                     |
| N_STAGE          | N0 / N1 / N2 / N3a / N3b                                    |
| TUMOR_SITE       | U (upper) / M (middle) / L (lower)                          |
| UICC_TUMOR_STAGE | IB / IIA / IIB / IIIA / IIIB / IIIC / IV; ~half missing     |
| SEX              | Male / Female                                               |


### Sample-level variables


| Variable          | Description                                                 |
| ----------------- | ----------------------------------------------------------- |
| GRADE             | Poor / Poor (signet); partial                               |
| HER2_IHC_SCORE    | 0 / 1; partial                                              |
| MUTATION_COUNT    | Somatic mutation count                                      |
| TMB_NONSYNONYMOUS | Nonsynonymous TMB (mut/Mb)                                  |
| RHOA_MUTATION     | Specific RHOA variant (Y42C, L22R, R5W, Y74D, R5Q); partial |
| ONCOTREE_CODE     | STAD                                                        |
| SOMATIC_STATUS    | Matched                                                     |
| T_STAGE           | T1b / T2 / T3 / T4a; partial                                |


### Notes

- **All diffuse-type gastric cancer** — dataset specifically designed to study diffuse/signet-ring histology
- **RHOA mutations** explicitly annotated; known driver in diffuse gastric cancer
- **No survival data** (no OS or DFS)
- Age given as grouped ranges, not exact values; ~half of patients have missing clinical annotations
- Bone scan result included (unusual variable; possibly for bone metastasis assessment)
- HER2 IHC score available for partially annotated cases
- Tumor site uses anatomical thirds (U/M/L) rather than named regions

---

---

## 6. egc_msk_2017

**N = 305 patients / 341 samples** | Esophagogastric Cancer (STAD + GEJ + ESCA) | MSK-IMPACT341 (targeted panel)  
MSK; all metastatic disease; PMID: 29122777

### Patient-level variables


| Variable                                            | Description                            |
| --------------------------------------------------- | -------------------------------------- |
| AGE_AT_DIAGNOSIS                                    | Age at diagnosis (decimal years)       |
| LAUREN_CLASS                                        | Intestinal / Diffuse / Mixed           |
| ECOG_AT_STAGE_4_DX                                  | ECOG performance status at stage IV dx |
| OS_MONTHS                                           | Overall survival in months             |
| OS_MONTHS_STAGE_IV                                  | OS from stage IV diagnosis             |
| VITAL_STATUS                                        | DOD / NED / AWD / LTF                  |
| RACE                                                | White / Asian / etc.                   |
| SEX                                                 | Male / Female                          |
| TUMOR_GRADE                                         | Mod / Poor                             |
| TREATMENT                                           | First-line regimen (free text)         |
| CLINICAL_TRIAL_1L/2L/3L                             | Clinical trial participation by line   |
| NUMBER_OF_LINES_OF_TREATMENT_FOR_MET_DISEASE        | Total lines of therapy                 |
| RECEIVED_IO_AFTER_1ST_LINE                          | Yes / No                               |
| IO_OS/PFS_MONTHS_START_OF_IO                        | IO-specific OS and PFS                 |
| IO_PFS_STATUS                                       | Progressed / etc.                      |
| PLATINUM_OS/PFS_MONTHS                              | Platinum-specific outcomes             |
| PLATINUM_BEST_RESPONSE                              | SD / PR / CR / etc.                    |
| PLATINUM_GERMLINE_TESTED / HRD_MUTATION / LST_SCORE | HRD biomarkers                         |
| EBV_TESTED                                          | YES / NO / INSUFFICIENT TISSUE         |


### Sample-level variables


| Variable                                     | Description                                                                |
| -------------------------------------------- | -------------------------------------------------------------------------- |
| CANCER_TYPE_DETAILED                         | Stomach Adenocarcinoma / Adenocarcinoma of GEJ / Esophageal Adenocarcinoma |
| ONCOTREE_CODE                                | STAD / GEJ / ESCA                                                          |
| GENE_PANEL                                   | IMPACT341                                                                  |
| MOLECULAR_SUBTYPE                            | CIN / MSI / EBV / GS                                                       |
| MSI_SENSOR_SCORE                             | Continuous MSI score                                                       |
| MUTATION_COUNT / TMB_NONSYNONYMOUS           | Mutation burden                                                            |
| FRACTION_GENOME_ALTERED (FGA)                | CNA burden                                                                 |
| HER2_FISH / HER2_IHC_OR_FISH / HER2_IMPACT   | HER2 status and amplification                                              |
| HER2_RESPONSE / HER2_TTP_ON_1ST_LINE_MONTHS  | HER2-targeted therapy outcomes                                             |
| H_PYLORI                                     | Positive / Negative                                                        |
| IO_PDL1_TESTED / IO_RESPONSE / IO_MUT_GROUPS | IO biomarkers and response                                                 |
| MET_LIVER / MET_LUNG / MET_PERITONEUM        | Metastatic site flags                                                      |
| SITE_OF_TISSUE_BIOPSY                        | Biopsy location                                                            |
| STAGE_AT_DIAGNOSIS                           | Stage at initial dx                                                        |
| SAMPLE_TYPE                                  | Primary / Metastasis                                                       |
| SOMATIC_STATUS                               | **Unmatched** (no paired normal)                                           |
| IMPACT_PRE_POST                              | Pre- or post-treatment biopsy                                              |


### Notes

- **All patients have metastatic disease**; treatment history richly annotated
- **Unmatched** tumor-only sequencing (no paired normal) — important for mutation calling interpretation
- HER2 therapy outcomes and platinum outcomes tracked as separate endpoints
- IO (immunotherapy) response and PDL1 testing status available
- Sample IDs shared with egc_mskcc_2020 and egc_msk_2023 (same P-XXXXXXX namespace)

---

## 7. egc_mskcc_2020

**N = 487 patients** | Esophagogastric Cancer (STAD + GEJ + ESCA) | MSK-IMPACT (targeted panel)  
MSK; PMID: 33795256

### Patient-level variables


| Variable                          | Description                   |
| --------------------------------- | ----------------------------- |
| COHORT                            | Cohort label                  |
| CSTAGE_CATEGORY / PSTAGE_CATEGORY | Clinical and pathologic stage |
| OS_MONTHS                         | Overall survival in months    |
| OS_STATUS                         | 0:LIVING / 1:DECEASED         |
| SEX                               | Male / Female                 |
| TUMOR_GRADE                       | Tumor grade                   |


### Sample-level variables


| Variable                           | Description                                              |
| ---------------------------------- | -------------------------------------------------------- |
| AGE_AT_SAMPLING                    | Age at time of biopsy                                    |
| CANCER_TYPE_DETAILED               | Stomach Adenocarcinoma / GEJ / Esophageal Adenocarcinoma |
| ONCOTREE_CODE                      | STAD / GEJ / ESCA                                        |
| MSI_SCORE                          | Continuous MSI sensor score                              |
| MUTATION_COUNT / TMB_NONSYNONYMOUS | Mutation burden                                          |
| FRACTION_GENOME_ALTERED (FGA)      | CNA burden                                               |
| TUMOR_PURITY                       | Estimated tumor purity                                   |
| WGD                                | Whole-genome doubling (0/1)                              |
| SAMPLE_TYPE                        | Primary / Metastasis                                     |
| SAMPLE_TYPE_DETAIL                 | e.g., Primary Esophagus / Lymph Node                     |
| SOMATIC_STATUS                     | Matched                                                  |


### Notes

- Minimal patient-level clinical annotation (stage, OS, grade only)
- Age captured at sample level (AGE_AT_SAMPLING), not patient level
- WGD flag and tumor purity available — useful for CNA analyses
- OS data available; no DFS or PFS

---

## 8. egc_msk_tp53_ccr_2022

**N = 237 patients** | Esophagogastric Cancer (STAD + GEJ + ESCA) | MSK-IMPACT (targeted panel)  
MSK; perioperative/neoadjuvant chemotherapy cohort; TP53-focused; PMID: 35377946

### Patient-level variables


| Variable                | Description                         |
| ----------------------- | ----------------------------------- |
| IMPACT_TMB_SCORE        | TMB score                           |
| OS_MONTHS / OS_STATUS   | Overall survival                    |
| PFS_MONTHS / PFS_STATUS | Progression-free survival           |
| SEX                     | Male / Female                       |
| TREATMENT_EFFECT        | Pathologic response score (numeric) |
| TUMOR_GRADE             | Tumor grade                         |


### Sample-level variables


| Variable                           | Description                                              |
| ---------------------------------- | -------------------------------------------------------- |
| CANCER_TYPE_DETAILED               | Stomach Adenocarcinoma / GEJ / Esophageal Adenocarcinoma |
| ONCOTREE_CODE                      | STAD / GEJ / ESCA                                        |
| CT_CATEGORY / PT_CATEGORY          | Clinical and pathologic T stage                          |
| CN_CATEGORY / PN_CATEGORY          | Clinical and pathologic N stage                          |
| CSTAGE_CATEGORY / PSTAGE_CATEGORY  | Clinical and pathologic overall stage                    |
| MSI_TYPE                           | Stable / Hypermutant                                     |
| MUTATION_COUNT / TMB_NONSYNONYMOUS | Mutation burden                                          |
| FRACTION_GENOME_ALTERED (FGA)      | CNA burden                                               |
| TP53_BIALLELIC_GROUPS              | TP53 biallelic status category                           |
| TRG_SCORE                          | Tumor regression grade (pathologic response)             |
| NT_RESPONDER                       | Near-total pathologic response (Yes/No)                  |
| PET_RESPONDER                      | PET imaging response (Yes/No)                            |
| PATHOLOGIC_DOWNSTAGING             | Yes / No                                                 |
| RESIDUAL_NODAL_DISEASE             | Yes / No                                                 |
| SURGERY_TYPE                       | Esophagectomy / etc.                                     |
| SAMPLE_TYPE                        | Primary / Metastasis                                     |
| TUMOR_PURITY / WGD                 | Purity and whole-genome doubling                         |
| SOMATIC_STATUS                     | Matched                                                  |


### Notes

- **Both OS and PFS** available
- Study focused on TP53 mutation status and biallelic inactivation in neoadjuvant-treated EGC
- Pathologic response richly annotated: TRG score, NT_responder, PET_responder, downstaging
- Pre-treatment biopsies from patients who received perioperative chemotherapy
- TP53_BIALLELIC_GROUPS distinguishes truncating biallelic, missense, etc.

---

## 9. egc_msk_2023

**N = 902 patients** | Esophagogastric Cancer (STAD + GEJ + ESCA) | MSK-IMPACT341 (targeted panel)  
MSK; early-onset vs average-onset EGC study; PMID: 37699004

### Patient-level variables


| Variable                       | Description                                                |
| ------------------------------ | ---------------------------------------------------------- |
| AGE_AT_DIAGNOSIS               | Age at diagnosis (decimal years)                           |
| AGE_CATEGORY                   | **Early Onset** (≤40) / **Average Onset** (>40)            |
| BMI_CATEGORIES                 | Underweight / Normal_Weight / Overweight / Obese / Unknown |
| ECOG_PS                        | ECOG performance status                                    |
| ETHNICITY                      | Hispanic or Latino / Not Hispanic or Latino                |
| HISTOLOGY                      | Adenocarcinoma / Signet_Diffuse                            |
| MET_NONMET_STATUS              | Metastatic / Non-Metastatic                                |
| OS_MONTHS / OS_MONTHS_STAGE_IV | Overall survival (total and from stage IV)                 |
| OS_STATUS                      | 0:LIVING / 1:DECEASED                                      |
| PRIMARY_SITE_TRI               | GEJ (Siewert I-II) / Gastric                               |
| RACE                           | White / Asian / Black / Other                              |
| SEX                            | Male / Female                                              |
| STAGE                          | I–IV                                                       |
| TIME_SX_DX_MONTHS              | Time from symptoms to diagnosis                            |


### Sample-level variables (107 columns total)


| Variable                                                                           | Description                                              |
| ---------------------------------------------------------------------------------- | -------------------------------------------------------- |
| CANCER_TYPE_DETAILED                                                               | Stomach Adenocarcinoma / GEJ / Esophageal Adenocarcinoma |
| ONCOTREE_CODE                                                                      | STAD / GEJ / ESCA                                        |
| UPDATED_TUMOR_SUBTYPE                                                              | CIN / MSI / EBV / GS / Other                             |
| MSI_SCORE / MSI_TYPE                                                               | MSI sensor score and classification                      |
| MSI_HYPERMUTANT                                                                    | Yes / No                                                 |
| MUTATION_COUNT / TMB_NONSYNONYMOUS                                                 | Mutation burden                                          |
| FRACTION_GENOME_ALTERED                                                            | CNA burden                                               |
| TUMOR_PURITY                                                                       | Estimated tumor purity                                   |
| GENE_PANEL                                                                         | IMPACT341                                                |
| SAMPLE_COVERAGE                                                                    | Sequencing depth                                         |
| SAMPLE_TYPE                                                                        | Primary / Metastasis                                     |
| SOMATIC_STATUS                                                                     | Matched                                                  |
| DENOVO_PRIMARY_METS                                                                | De novo metastatic (Yes/No)                              |
| LIVER_METS / LUNG_METS / PERITONEUM_METS / LYMPHNODE_METS / BONE_METS / BRAIN_METS | Metastatic site flags                                    |
| SYMPTOMATIC_METS                                                                   | Yes / No                                                 |
| WEIGHT_LOSS / DYSPHAGIA / DECR_APPETITE / NAUSEA_VOMITING / ABDONIMAL_PAIN         | Symptom flags                                            |
| GERD / BARRETTS_ESOPHAGUS / H_PYLORI (via GERD field)                              | GI risk factors                                          |
| MULTIPLE_CANCERS / ALTERNATIVE_DX_BINARY                                           | Second malignancy flags                                  |
| USED_FOR_GENOMICS / USED_FOR_SURVIVAL / USE_FOR_SYMPTOM_ANALYSIS                   | Analysis-eligibility flags                               |
| HTN / DM / CAD / CKD / COPD / HLD / etc.                                           | Comorbidity flags                                        |


### Notes

- **Largest MSK dataset** (N=902); primary focus is early-onset EGC biology
- Most extensively annotated sample file (107 columns) — symptoms, comorbidities, metastatic sites
- AGE_CATEGORY is the key stratification variable (Early Onset ≤40 vs Average Onset >40)
- TIME_SX_DX_MONTHS captures delay to diagnosis (relevant for early-onset biology)
- Molecular subtypes (UPDATED_TUMOR_SUBTYPE) available for all samples
- OS available; no DFS/PFS

---

## Cross-Dataset Comparison


| Feature                 | egc_tmucih_2015   | stad_oncosg_2018     | stad_pfizer_uhongkong | stad_uhongkong | stad_utokyo       | egc_msk_2017         | egc_mskcc_2020 | egc_msk_tp53_ccr_2022 | egc_msk_2023         |
| ----------------------- | ----------------- | -------------------- | --------------------- | -------------- | ----------------- | -------------------- | -------------- | --------------------- | -------------------- |
| N patients              | 78                | ~146                 | 100                   | 22             | 30                | 305                  | 487            | 237                   | 902                  |
| Sequencing              | WXS               | WGS/WXS              | WGS                   | WXS            | WXS               | IMPACT panel         | IMPACT panel   | IMPACT panel          | IMPACT panel         |
| Pure gastric only       | Yes               | Yes                  | Yes                   | Yes            | Yes               | No (EGC)             | No (EGC)       | No (EGC)              | No (EGC)             |
| OS data                 | Yes               | No                   | No                    | Yes            | No                | Yes                  | Yes            | Yes                   | Yes                  |
| DFS/PFS data            | No                | No                   | No                    | DFS            | No                | IO/Platinum PFS      | No             | PFS                   | No                   |
| Lauren's classification | Indirect          | Yes                  | Yes                   | Yes            | Yes (all Diffuse) | Yes                  | No             | No                    | Histology field      |
| Molecular subtype       | Partial           | Yes (EBV/MSI/CIN/GS) | MSI only              | MSI only       | No                | Yes (CIN/MSI/EBV/GS) | No             | MSI type              | Yes (CIN/MSI/EBV/GS) |
| H. pylori               | No                | No                   | Yes                   | Yes            | No                | Yes (sample)         | No             | No                    | No                   |
| EBV status              | Yes               | Partial              | Yes                   | Yes            | No                | Yes (sample)         | No             | No                    | No                   |
| TMB                     | Yes               | Yes                  | Yes                   | Yes            | Yes               | Yes                  | Yes            | Yes                   | Yes                  |
| HER2                    | No                | No                   | No                    | No             | Partial           | Yes (detailed)       | No             | No                    | No                   |
| Treatment history       | No                | No                   | No                    | No             | No                | Yes (rich)           | No             | Partial               | No                   |
| IO response data        | No                | No                   | No                    | No             | No                | Yes                  | No             | No                    | No                   |
| Pathologic response     | No                | No                   | No                    | No             | No                | No                   | No             | Yes (TRG/NT/PET)      | No                   |
| Comorbidities/symptoms  | No                | No                   | No                    | No             | No                | No                   | No             | No                    | Yes (extensive)      |
| Age at onset            | Exact             | Partial              | Exact                 | Exact          | Range only        | Exact                | Sample-level   | No                    | Yes + category       |
| Race / Ethnicity        | No                | No                   | No                    | No             | No                | Race                 | No             | No                    | Race + Ethnicity     |
| Tumor purity            | No                | No                   | Yes                   | Yes            | No                | No                   | Yes            | Yes                   | Yes                  |
| Matched normal          | Yes               | Yes                  | Yes                   | Yes            | Yes               | **No**               | Yes            | Yes                   | Yes                  |
| RHOA mutation           | No                | No                   | No                    | No             | Yes               | No                   | No             | No                    | No                   |
| Immune deconvolution    | No                | Partial (xCell)      | No                    | No             | No                | No                   | No             | No                    | No                   |
| Population              | Chinese (Beijing) | Multi (SG/HK)        | Chinese (HK)          | Chinese (HK)   | Japanese          | Mixed (US)           | Mixed (US)     | Mixed (US)            | Mixed (US)           |


