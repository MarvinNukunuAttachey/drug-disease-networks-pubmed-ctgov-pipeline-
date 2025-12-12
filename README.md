# Drug–Disease Networks: PubMed + ClinicalTrials.gov Pipeline

This repository contains an end-to-end R pipeline to extract, harmonize, and analyze **drug–disease treatment relationships** from:

- **PubMed** (via `rentrez` and MeSH treatment annotations)
- **ClinicalTrials.gov v2 API** (via `httr`/`jsonlite`)

The focus is on:
- Mapping **which diseases each drug treats**
- Quantifying **how strong** the treatment evidence is (e.g., number of PMIDs / trials)
- Building **drug–disease** and **drug–drug co-treatment** **networks**
- Summarizing **trial phases, status, and temporal trends**

---

## Pipeline Overview

### 1. Shared Drug List

The analysis starts from a curated list of commonly prescribed medications  
(e.g., atorvastatin, metformin, lisinopril, furosemide, tirzepatide, etc.).  
This same list is used across **both** data sources to keep everything aligned.

---

## Part I – PubMed Pipeline

**File:** `pubmed_pipeline.R` (or the corresponding section in your main script)

### Data Source
- PubMed, queried via the **`rentrez`** package.
- Uses **MeSH terms** + subheadings to focus on **treatment-related** articles:
  - `"Therapeutic Use"[Subheading]`
  - `"Drug Therapy"[Subheading]`
  - `"Treatment Outcome"[MeSH Terms]`
  - `"clinical trial"[Publication Type]`

### Key Steps

1. **Query PubMed** for each drug using `entrez_search` and `entrez_fetch`.
2. **Extract MeSH metadata**:
   - PMID  
   - Article title  
   - MeSH descriptor terms  
   - Publication year  
   - Publication types  
3. **Identify candidate diseases** from MeSH terms:
   - Keep terms containing patterns like `Disease`, `Syndrome`, `Disorder`, `Cancer`, `Infection`, `Failure`, `Diabetes`, `Pain`, `Stroke`, `Hypertension`.
   - Drop non-disease concepts such as `Animals`, `Healthy Volunteers`, `Method`, `Therapy`, etc.
4. Build a **drug–disease edge list**:
   - One row per `(drug, disease)` pair with **weight = number of supporting PMIDs**.

### PubMed Outputs & Visuals

- **Top diseases** by total PMID evidence  
- **Drug breadth**:
  - Number of unique diseases per drug  
  - Total PMIDs per drug  
- **Scatter plot** of:
  - X = number of diseases treated  
  - Y = total PMIDs  
- **Weighted network metrics**:
  - Degree / strength for each drug  
- **Disease co-annotation density**:
  - How many drugs are mapped to each disease  
- **Publication year distribution**:
  - Temporal trend of treatment evidence  
- **Top publication types**:
  - Analogue to trial status categories  
- **Networks**:
  - Drug–disease bipartite network (using `igraph` + `ggraph`)  
  - Drug–drug co-treatment network:
    - Two drugs are connected if they treat at least one **shared disease**
    - Edge weight = shared evidence (min of weights, aggregated across diseases)

---

## Part II – ClinicalTrials.gov Pipeline

**File:** `ctgov_pipeline.R` (or the corresponding section in your main script)

### Data Source
- **ClinicalTrials.gov v2 API**, accessed via `httr` + `jsonlite`.

### Key Steps

1. For each drug, query `https://clinicaltrials.gov/api/v2/studies` using:
   - `query.intr = <drug>`
   - Pagination with `pageToken`.
2. For each study, extract:
   - NCT ID  
   - Brief title  
   - Overall status (e.g., Completed, Recruiting)  
   - Phases (e.g., Phase 2, Phase 3)  
   - Conditions (diseases)  
   - First submit date  
3. Clean and **unnest conditions**:
   - Split multiple conditions per study  
   - Drop generic non-disease phrases:
     - e.g., `Healthy Volunteer`, `Progression`, `Participant`, `Study`
4. Build a **drug–disease edge list**:
   - One row per `(drug, disease)` with **weight = number of trials**.

### ClinicalTrials.gov Outputs & Visuals

- **Top diseases** by total trial count  
- **Top drugs** by:
  - Total number of trials  
  - Breadth: number of unique diseases per drug  
- **Breadth vs total trials** scatter for top drugs  
- **Distribution of trial phases** (Phase 1–4, etc.)  
- **Trial status distribution** (Completed, Recruiting, etc.)  
- **Disease co-annotation density**:
  - Number of drugs per disease  
- **Yearly trend**:
  - Number of trials by first submit year  
- **Networks**:
  - Drug–disease bipartite network  
  - Drug–drug co-treatment network for top drugs  

---

## Network Analysis

For both PubMed and ClinicalTrials.gov, the pipeline computes network-level stats using **`igraph`**:

- Total nodes and edges  
- Number of drugs vs diseases  
- Density  
- Number of connected components  
- Giant component size  
- Average degree and median degree  
- Average path length  
- Diameter  
- For drug–drug networks:
  - Weighted degree (co-treatment strength)

These metrics let us compare how **literature-based** evidence (PubMed) vs **trial-based** evidence (CT.gov) organize the same set of medications in the disease space.

---

## PubTator RAG Integration (Submodule)

This project integrates the **[`pubtator-rag`](https://github.com/Adi-M02/pubtator-rag)** repository as a **git submodule** to support PubTator-based retrieval, annotation, and graph construction.

The submodule is located in the folder: `pubtator-rag/`

### Updating the Submodule

When the original repository is updated, you can sync the latest changes with:

```bash
cd pubtator-rag
git pull origin main   # or 'master' if that is the default branch
cd ..
git add pubtator-rag
git commit -m "Update pubtator-rag to latest upstream version"
git push

---

## Credits

Special thanks to **Adi** (https://github.com/Adi-M02) for developing and maintaining  
the original PubTator RAG framework used in this project.

---

## Requirements

R (≥ 4.0) with the following packages:

```r
install.packages(c(
  "dplyr", "tidyr", "purrr", "stringr",
  "rentrez", "xml2", "httr", "jsonlite",
  "igraph", "ggraph", "ggplot2", "forcats",
  "ggrepel"
))

```
Author

Marvin Nukunu-Attachey
Graduate Researcher – Health Informatics, University of Iowa

Email: marvindee98@gmail.com
 · mnukunuattachey@uiowa.edu

GitHub: https://github.com/MarvinNukunuAttachey

Research Interests:
Health Informatics, Interpretable Machine Learning
