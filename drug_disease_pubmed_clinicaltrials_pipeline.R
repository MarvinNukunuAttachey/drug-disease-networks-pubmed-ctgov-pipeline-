###############################################################
# Unified PubMed + ClinicalTrials.gov Drug–Disease Pipeline
# Author: Marvin Nukunu-Attachey
# Focus: Treatment relationships + co-annotation + networks
###############################################################
file.remove("pubmed_results.csv")
file.remove("ctgov_results.csv")
# ================== Libraries ==================
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(rentrez)
library(xml2)
library(igraph)
library(ggplot2)
library(forcats)
library(httr)
library(jsonlite)
library(ggraph)

USER_EMAIL <- "enter-email here"

# Drug list
drugs <- c(
  "Atorvastatin", "Metformin", "Levothyroxine", "Lisinopril", "Amlodipine",
  "Metoprolol", "Albuterol", "Losartan", "Gabapentin", "Omeprazole",
  "Sertraline", "Rosuvastatin", "Pantoprazole", "Escitalopram",
  "Dextroamphetamine; Dextroamphetamine Saccharate; Amphetamine; Amphetamine Aspartate",
  "Hydrochlorothiazide", "Bupropion", "Fluoxetine", "Semaglutide", "Montelukast",
  "Trazodone", "Simvastatin", "Amoxicillin", "Tamsulosin",
  "Acetaminophen; Hydrocodone", "Fluticasone", "Meloxicam", "Apixaban",
  "Furosemide", "Insulin Glargine", "Duloxetine", "Ibuprofen", "Famotidine",
  "Empagliflozin", "Carvedilol", "Tramadol", "Alprazolam", "Prednisone",
  "Hydroxyzine", "Buspirone", "Clopidogrel", "Glipizide", "Citalopram",
  "Potassium Chloride", "Allopurinol", "Aspirin", "Cyclobenzaprine",
  "Ergocalciferol", "Oxycodone", "Methylphenidate", "Venlafaxine",
  "Spironolactone", "Ondansetron", "Zolpidem", "Cetirizine", "Estradiol",
  "Pravastatin", "Hydrochlorothiazide; Lisinopril", "Lamotrigine",
  "Quetiapine", "Fluticasone; Salmeterol", "Clonazepam", "Dulaglutide",
  "Azithromycin", "Hydrochlorothiazide; Losartan", "Amoxicillin; Clavulanate",
  "Latanoprost", "Cholecalciferol", "Propranolol", "Ezetimibe",
  "Topiramate", "Paroxetine", "Diclofenac", "Budesonide; Formoterol",
  "Atenolol", "Lisdexamfetamine", "Doxycycline", "Pregabalin",
  "Ethinyl Estradiol; Norethindrone", "Glimepiride", "Tizanidine",
  "Clonidine", "Fenofibrate", "Insulin Lispro", "Valsartan", "Cephalexin",
  "Baclofen", "Rivaroxaban", "Ferrous Sulfate", "Amitriptyline",
  "Finasteride", "Dapagliflozin", "Acetaminophen; Oxycodone", "Folic Acid",
  "Aripiprazole", "Olmesartan", "Ethinyl Estradiol; Norgestimate",
  "Valacyclovir", "Mirtazapine", "Lorazepam", "Levetiracetam",
  "Insulin Aspart", "Naproxen", "Cyanocobalamin", "Loratadine",
  "Diltiazem", "Sumatriptan", "Triamcinolone", "Hydralazine",
  "Tirzepatide", "Celecoxib", "Acetaminophen", "Alendronate",
  "Oxybutynin", "Hydrochlorothiazide; Triamterene", "Warfarin",
  "Progesterone", "Fluticasone; Umeclidinium; Vilanterol", "Testosterone",
  "Nifedipine", "Methocarbamol", "Benzonatate", "Sitagliptin",
  "Chlorthalidone", "Isosorbide", "Donepezil", "Dexmethylphenidate",
  "Sulfamethoxazole; Trimethoprim", "Clobetasol", "Methotrexate",
  "Hydroxychloroquine", "Lovastatin", "Pioglitazone", "Irbesartan",
  "Methylprednisolone", "Norethindrone", "Meclizine",
  "Ethinyl Estradiol; Levonorgestrel", "Fluticasone; Vilanterol",
  "Ketoconazole", "Thyroid", "Azelastine", "Nitrofurantoin",
  "Adalimumab", "Memantine", "Prednisolone", "Esomeprazole", "Docusate",
  "Clindamycin", "Acyclovir"
)

###############################################################
# I. PubMed via rentrez
###############################################################

# ---- Helper: extract metadata including year & publication type ----
extract_mesh_metadata <- function(xml_txt) {
  doc <- read_xml(xml_txt)
  articles <- xml_find_all(doc, ".//PubmedArticle")
  
  map_dfr(articles, function(a) {
    pmid  <- xml_text(xml_find_first(a, ".//PMID"))
    title <- xml_text(xml_find_first(a, ".//ArticleTitle"))
    mesh  <- xml_text(xml_find_all(a, ".//MeshHeading/DescriptorName"))
    
    year_node <- xml_find_first(a, ".//Article/Journal/JournalIssue/PubDate/Year")
    year_val  <- xml_text(year_node)
    if (is.na(year_val) || year_val == "") {
      year_val <- NA
    }
    
    
    ptypes <- xml_text(xml_find_all(a, ".//PublicationType"))
    
    tibble(
      pmid       = pmid,
      title      = title,
      mesh_terms = if (length(mesh)) paste(mesh, collapse = "; ") else NA_character_,
      pub_year   = year_val,
      pub_types  = if (length(ptypes)) paste(ptypes, collapse = "; ") else NA_character_
    )
  })
}

# ---- Fetch PubMed data for one drug ----
fetch_pubmed <- function(drug, retmax = 120, sleep_sec = 0.6) {
  query <- sprintf('("%s"[MeSH Terms]) AND 
                   ("Therapeutic Use"[Subheading] OR "Drug Therapy"[Subheading] 
                   OR "Treatment Outcome"[MeSH Terms])
                   AND "clinical trial"[Publication Type]', drug)
  message("🔍 PubMed: ", drug)
  
  srch <- tryCatch(
    entrez_search(db = "pubmed", term = query, retmax = retmax, email = USER_EMAIL),
    error = function(e) NULL
  )
  if (is.null(srch) || length(srch$ids) == 0)
    return(tibble(drug, disease = NA, pmid = NA, positive = FALSE,
                  pub_year = NA, pub_types = NA))
  
  batches <- split(srch$ids, ceiling(seq_along(srch$ids) / 10))
  records <- map_dfr(batches, function(ids) {
    txt <- tryCatch(
      entrez_fetch(db = "pubmed", id = paste(ids, collapse = ","), rettype = "xml"),
      error = function(e) NULL
    )
    Sys.sleep(sleep_sec)
    if (is.null(txt)) return(tibble())
    extract_mesh_metadata(txt)
  })
  
  # Extract candidate diseases from MeSH
  disease_df <- records %>%
    filter(!is.na(mesh_terms)) %>%
    separate_rows(mesh_terms, sep = ";") %>%
    mutate(mesh_terms = str_trim(mesh_terms)) %>%
    filter(
      str_detect(mesh_terms,
                 "(Disease|Syndrome|Disorder|Cancer|Infection|Failure|Diabetes|Pain|Stroke|Hypertension)"),
      !str_detect(mesh_terms,
                  "Animals|Healthy|Volunteers|Method|Study|Therapy")
    ) %>%
    distinct(pmid, mesh_terms, .keep_all = TRUE) %>%
    rename(disease = mesh_terms)
  
  tibble(
    drug      = drug,
    disease   = disease_df$disease,
    pmid      = disease_df$pmid,
    positive  = TRUE,
    pub_year  = disease_df$pub_year,
    pub_types = disease_df$pub_types
  )
}

# ---- Get or load PubMed data ----
if (file.exists("pubmed_results.csv")) {
  pub_df <- read.csv("pubmed_results.csv", stringsAsFactors = FALSE)
} else {
  pub_df <- map_df(drugs, fetch_pubmed)
  write.csv(pub_df, "pubmed_results.csv", row.names = FALSE)
}

pub_df <- pub_df %>% filter(!is.na(disease))

# ================== PubMed: summaries ==================

# Edges: drug–disease with weight = # PMIDs
edges_pub <- pub_df %>% count(drug, disease, name = "weight")

# Top diseases by PMID support
pub_top_diseases <- edges_pub %>%
  group_by(disease) %>%
  summarise(pmids = sum(weight), .groups = "drop") %>%
  slice_max(pmids, n = 20)

ggplot(pub_top_diseases,
       aes(x = fct_reorder(disease, pmids), y = pmids)) +
  geom_col(fill = "#d62728", alpha = .9) +
  coord_flip() +
  labs(title = "Top 20 Diseases Treated (PubMed)",
       x = "Disease", y = "Number of Supporting PMIDs") +
  theme_minimal(base_size = 13)

# Drug breadth: how many unique diseases
pub_drug_breadth <- edges_pub %>%
  group_by(drug) %>%
  summarise(unique_diseases = n_distinct(disease),
            total_pmids     = sum(weight), .groups = "drop")

pub_top_breadth <- pub_drug_breadth %>%
  slice_max(unique_diseases, n = 20)

ggplot(pub_top_breadth,
       aes(x = fct_reorder(drug, unique_diseases), y = unique_diseases)) +
  geom_col(fill = "#1f77b4", alpha = .9) +
  coord_flip() +
  labs(title = "Top 20 Drugs by Number of Diseases Treated (PubMed)",
       x = "Drug", y = "Unique Diseases") +
  theme_minimal(base_size = 13)

ggplot(pub_top_breadth,
       aes(x = unique_diseases,
           y = total_pmids,
           color = unique_diseases,
           size = total_pmids,
           label = drug)) +
  geom_point(alpha = 0.9) +
  scale_color_viridis_c() +
  geom_text(vjust = -0.6, size = 2.8, alpha = 0.9) +
  scale_size(range = c(2, 10), guide = "legend") +
  labs(
    title = "Drug Breadth vs Total Evidence (PubMed)",
    x = "Unique Diseases Treated",
    y = "Total PMIDs (Treatment Evidence)",
    color = "Unique Diseases",
    size = "Total PMIDs"
  ) +
  theme_minimal(base_size = 13)

# Weighted evidence: sum of PMID-weighted edges
g_pub <- graph_from_data_frame(edges_pub, directed = FALSE)
V(g_pub)$group <- ifelse(V(g_pub)$name %in% drugs, "Drug", "Disease")

valid_drugs_pub <- intersect(drugs, V(g_pub)$name)

pub_drug_metrics <- tibble(
  drug           = valid_drugs_pub,
  degree         = degree(g_pub, v = valid_drugs_pub),
  weighted_degree = strength(g_pub, vids = valid_drugs_pub,
                             weights = E(g_pub)$weight)
)

pub_top_weighted <- pub_drug_metrics %>%
  slice_max(weighted_degree, n = 20)

ggplot(pub_top_weighted,
       aes(x = fct_reorder(drug, weighted_degree), y = weighted_degree)) +
  geom_col(fill = "#1f77b4", alpha = .9) +
  coord_flip() +
  labs(title = "Top 20 Drugs by Weighted Evidence (PubMed)",
       x = "Drug",
       y = "Weighted Degree (Sum of PMID-Weighted Links)") +
  theme_minimal(base_size = 13)

# Disease co-annotation density (how many drugs treat each disease)
pub_disease_density <- edges_pub %>%
  group_by(disease) %>%
  summarise(n_drugs = n_distinct(drug), .groups = "drop") %>%
  mutate(category = ifelse(n_drugs > median(n_drugs),
                           "high-density", "low-density"))

ggplot(pub_disease_density,
       aes(x = n_drugs, fill = category)) +
  geom_histogram(bins = 25, alpha = .8, position = "identity") +
  scale_fill_manual(values = c("high-density" = "#d62728",
                               "low-density"  = "#1f77b4")) +
  labs(title = "Distribution of Disease Co-annotation Density (PubMed)",
       x = "Number of Drugs Treating the Disease",
       y = "Count") +
  theme_minimal(base_size = 13)

# Publication year distribution (phase-style analogue)
pub_year_counts <- pub_df %>%
  filter(!is.na(pub_year), pub_year != "") %>%
  distinct(pmid, pub_year) %>%
  count(pub_year) %>%
  mutate(pub_year = as.integer(pub_year)) %>%
  arrange(pub_year)

ggplot(pub_year_counts,
       aes(x = pub_year, y = n)) +
  geom_col(fill = "steelblue", alpha = .8) +
  labs(title = "Distribution of Publication Years (PubMed)",
       x = "Publication Year",
       y = "Number of Articles") +
  theme_minimal(base_size = 13)

# Publication type distribution (status-style analogue)
pub_type_counts <- pub_df %>%
  filter(!is.na(pub_types)) %>%
  distinct(pmid, pub_types) %>%
  separate_rows(pub_types, sep = ";") %>%
  mutate(pub_types = str_trim(pub_types)) %>%
  count(pub_types, sort = TRUE) %>%
  slice_max(n, n = 10)

ggplot(pub_type_counts,
       aes(x = fct_reorder(pub_types, n), y = n)) +
  geom_col(fill = "steelblue", alpha = .8) +
  coord_flip() +
  labs(title = "Top Publication Types (PubMed)",
       x = "Publication Type",
       y = "Number of Articles") +
  theme_minimal(base_size = 13)

# ================== PubMed networks ==================

# ---- Drug–disease bipartite network (PubMed) ----
set.seed(123)
layout_pub <- create_layout(g_pub, layout = "fr")

drug_nodes_pub <- as.data.frame(layout_pub) %>%
  filter(name %in% drugs)

ggraph(layout_pub) +
  geom_edge_link(aes(width = weight), colour = "grey80", alpha = 0.6) +
  geom_node_point(aes(colour = group), size = 3) +
  geom_text(
    data = drug_nodes_pub,
    aes(x = x, y = y, label = name),
    size = 3,
    vjust = -0.4
  ) +
  scale_colour_manual(values = c("Drug" = "#1f77b4", "Disease" = "#d62728")) +
  scale_edge_width(range = c(0.2, 2)) +
  labs(title = "Drug–Disease Treatment Network (PubMed)",
       colour = "Node Type",
       width  = "PMID Support") +
  theme_void()

# ---- Drug–drug co-treatment network (PubMed) ----
drug_disease_list_pub <- split(edges_pub, edges_pub$disease)

drug_drug_edges_pub <- map_dfr(drug_disease_list_pub, function(df) {
  if (nrow(df) < 2) return(NULL)
  combos <- t(combn(df$drug, 2))
  tibble(
    drug1 = combos[, 1],
    drug2 = combos[, 2],
    w1    = df$weight[match(combos[, 1], df$drug)],
    w2    = df$weight[match(combos[, 2], df$drug)],
    edge_weight = pmin(w1, w2)
  )
}) %>%
  group_by(drug1, drug2) %>%
  summarise(weight = sum(edge_weight), .groups = "drop")

g_dd_pub <- graph_from_data_frame(drug_drug_edges_pub, directed = FALSE)

dd_metrics_pub <- tibble(
  drug           = V(g_dd_pub)$name,
  weighted_degree = strength(g_dd_pub, weights = E(g_dd_pub)$weight)
)

top20_dd_pub <- dd_metrics_pub %>%
  slice_max(weighted_degree, n = 20)

# Plot: top 20 drugs in co-treatment network
ggplot(top20_dd_pub,
       aes(x = fct_reorder(drug, weighted_degree), y = weighted_degree)) +
  geom_col(fill = "#ff7f0e", alpha = .9) +
  coord_flip() +
  labs(title = "Top 20 Drugs in Drug–Drug Co-treatment Network (PubMed)",
       x = "Drug",
       y = "Weighted Degree (Co-treatment Strength)") +
  theme_minimal(base_size = 13)

# Network plot restricted to top 20 drugs
g_dd_pub_top <- induced_subgraph(g_dd_pub, vids = top20_dd_pub$drug)
set.seed(123)
layout_dd_pub <- create_layout(g_dd_pub_top, layout = "fr")

ggraph(layout_dd_pub) +
  geom_edge_link(aes(width = weight), colour = "grey80", alpha = 0.6) +
  geom_node_point(colour = "#ff7f0e", size = 3) +
  geom_text(aes(x = x, y = y, label = name), vjust = -0.4, size = 3) +
  scale_edge_width(range = c(0.2, 2)) +
  labs(title = "Drug–Drug Co-treatment Network (PubMed, Top 20 Drugs)",
       width = "Co-treatment Weight") +
  theme_void()


pub_stats <- list(
  total_nodes         = gorder(g_pub),
  total_edges         = gsize(g_pub),
  drugs               = sum(V(g_pub)$name %in% drugs),
  diseases            = sum(!(V(g_pub)$name %in% drugs)),
  density             = edge_density(g_pub),
  components          = count_components(g_pub),
  giant_component_size = max(components(g_pub)$csize),
  avg_deg             = mean(degree(g_pub)),
  median_deg          = median(degree(g_pub)),
  avg_path            = mean_distance(g_pub, directed = FALSE),
  diameter            = diameter(g_pub, directed = FALSE)
)


pub_stats

pub_dd_stats <- list(
  total_nodes = gorder(g_dd_pub),
  total_edges = gsize(g_dd_pub),
  density = edge_density(g_dd_pub),
  components = count_components(g_dd_pub),
  giant_component = max(components(g_dd_pub)$csize),
  avg_deg = mean(degree(g_dd_pub)),
  diameter = diameter(g_dd_pub)
)

pub_dd_stats

###############################################################
# II. ClinicalTrials.gov v2
###############################################################

fetch_ctgov_v2_studies <- function(drug, max_pages = 3, page_size = 100) {
  message("🔍 ClinicalTrials.gov: ", drug)
  
  base_url <- "https://clinicaltrials.gov/api/v2/studies"
  all_records <- list()
  next_page <- NULL
  
  for (i in seq_len(max_pages)) {
    params <- list(
      format = "json",
      `query.intr` = drug,
      fields = paste(
        "protocolSection.identificationModule.nctId",
        "protocolSection.identificationModule.briefTitle",
        "protocolSection.statusModule.overallStatus",
        "protocolSection.conditionsModule.conditions",
        "protocolSection.designModule.phases",
        sep = ","
      ),
      pageSize = page_size
    )
    if (!is.null(next_page)) params$pageToken <- next_page
    
    res <- tryCatch(GET(base_url, query = params), error = function(e) NULL)
    if (is.null(res) || status_code(res) != 200) break
    
    js <- fromJSON(content(res, as = "text", encoding = "UTF-8"), flatten = TRUE)
    if (is.null(js$studies) || length(js$studies) == 0) break
    
    df <- tibble(
      drug   = drug,
      nct_id = sapply(js$studies$protocolSection.identificationModule.nctId, `[`, 1),
      title  = sapply(js$studies$protocolSection.identificationModule.briefTitle, `[`, 1),
      status = sapply(js$studies$protocolSection.statusModule.overallStatus, `[`, 1),
      phase  = sapply(js$studies$protocolSection.designModule.phases,
                      function(x) ifelse(length(x) == 0, NA, paste(x, collapse = "; "))),
      disease = sapply(js$studies$protocolSection.conditionsModule.conditions,
                       function(x) ifelse(length(x) == 0, NA, paste(x, collapse = "; ")))
    )
    
    all_records[[i]] <- df
    next_page <- js$nextPageToken
    if (is.null(next_page)) break
  }
  
  bind_rows(all_records)
}

# ---- Get or load ClinicalTrials.gov data ----
if (file.exists("ctgov_results.csv")) {
  ct_df <- read.csv("ctgov_results.csv", stringsAsFactors = FALSE)
} else {
  ct_df <- map_df(drugs, fetch_ctgov_v2_studies)
  write.csv(ct_df, "ctgov_results.csv", row.names = FALSE)
}

ct_df <- ct_df %>%
  filter(!is.na(disease) & disease != "") %>%
  filter(!str_detect(disease,
                     regex("Healthy|Volunteer|Progression|Study|Participant",
                           ignore_case = TRUE)))

# For disease-level analysis, split multiple conditions
ct_long <- ct_df %>%
  separate_rows(disease, sep = ";") %>%
  mutate(disease = str_trim(disease)) %>%
  filter(disease != "")

# ================== ClinicalTrials.gov: summaries ==================

edges_ct <- ct_long %>% count(drug, disease, name = "weight")

# Top diseases
ct_top_diseases <- edges_ct %>%
  group_by(disease) %>%
  summarise(trials = sum(weight), .groups = "drop") %>%
  slice_max(trials, n = 20)

ggplot(ct_top_diseases,
       aes(x = fct_reorder(disease, trials), y = trials)) +
  geom_col(fill = "#d62728", alpha = .9) +
  coord_flip() +
  labs(title = "Top 20 Diseases (ClinicalTrials.gov)",
       x = "Disease", y = "Number of Trials") +
  theme_minimal(base_size = 13)

# Top drugs
ct_drug_breadth <- edges_ct %>%
  group_by(drug) %>%
  summarise(unique_diseases = n_distinct(disease),
            total_trials    = sum(weight), .groups = "drop")

ct_top_drugs <- ct_drug_breadth %>%
  slice_max(total_trials, n = 20)

ggplot(ct_top_drugs,
       aes(x = fct_reorder(drug, total_trials), y = total_trials)) +
  geom_col(fill = "#1f77b4", alpha = .9) +
  coord_flip() +
  labs(title = "Top 20 Drugs by Trial Count (ClinicalTrials.gov)",
       x = "Drug", y = "Total Trials") +
  theme_minimal(base_size = 13)

# Breadth vs total trials (scatter)
ggplot(ct_drug_breadth,
       aes(x = unique_diseases, y = total_trials, label = drug)) +
  geom_point(aes(size = total_trials, colour = unique_diseases), alpha = .8) +
  geom_text(vjust = -0.6, size = 3) +
  scale_color_viridis_c() +
  labs(title = "Drug Breadth vs Total Trials (ClinicalTrials.gov)",
       x = "Unique Diseases",
       y = "Total Trials") +
  theme_minimal(base_size = 13)


# --- Top 20 drugs by total number of trials ---
ct_top20 <- ct_drug_breadth %>%
  arrange(desc(total_trials)) %>%
  slice_head(n = 20)

# --- Breadth vs total trials (scatter, top 20 only) ---
ggplot(ct_top20,
       aes(x = unique_diseases, y = total_trials, label = drug)) +
  geom_point(aes(size = total_trials, colour = unique_diseases), alpha = .8) +
  geom_text(vjust = -0.6, size = 3) +
  scale_color_viridis_c() +
  labs(title = "Drug Breadth vs Total Trials (ClinicalTrials.gov, Top 20)",
       x = "Unique Diseases",
       y = "Total Trials") +
  theme_minimal(base_size = 13)

# Distribution of trial phases
ct_phase_counts <- ct_df %>%
  filter(!is.na(phase)) %>%
  separate_rows(phase, sep = ";") %>%
  mutate(phase = str_trim(phase)) %>%
  count(phase, sort = TRUE)

ggplot(ct_phase_counts,
       aes(x = fct_reorder(phase, n), y = n)) +
  geom_col(fill = "lightblue", alpha = .9) +
  coord_flip() +
  labs(title = "Distribution of Clinical Trial Phases (All Drugs)",
       x = "Phase", y = "Number of Studies") +
  theme_minimal(base_size = 13)

# Distribution of trial status
ct_status_counts <- ct_df %>%
  filter(!is.na(status)) %>%
  count(status, sort = TRUE)

ggplot(ct_status_counts,
       aes(x = fct_reorder(status, n), y = n)) +
  geom_col(fill = "lightblue", alpha = .9) +
  coord_flip() +
  labs(title = "Trial Status Distribution (Across All Drugs)",
       x = "Trial Status", y = "Count") +
  theme_minimal(base_size = 13)

# Disease co-annotation density (trials)
ct_disease_density <- edges_ct %>%
  group_by(disease) %>%
  summarise(n_drugs = n_distinct(drug), .groups = "drop") %>%
  mutate(category = ifelse(n_drugs > median(n_drugs),
                           "high-density", "low-density"))

ggplot(ct_disease_density,
       aes(x = n_drugs, fill = category)) +
  geom_histogram(bins = 25, alpha = .8, position = "identity") +
  scale_fill_manual(values = c("high-density" = "#d62728",
                               "low-density"  = "#1f77b4")) +
  labs(title = "Distribution of Disease Co-annotation Density (ClinicalTrials.gov)",
       x = "Number of Drugs Treating the Disease",
       y = "Count") +
  theme_minimal(base_size = 13)

# ================== ClinicalTrials.gov networks ==================

# ---- Drug–disease network (CT.gov) ----
g_ct <- graph_from_data_frame(edges_ct, directed = FALSE)
V(g_ct)$group <- ifelse(V(g_ct)$name %in% drugs, "Drug", "Disease")

set.seed(123)
layout_ct <- create_layout(g_ct, layout = "fr")

drug_nodes_ct <- as.data.frame(layout_ct) %>%
  filter(name %in% drugs)

ggraph(layout_ct) +
  geom_edge_link(aes(width = weight), colour = "grey80", alpha = 0.6) +
  geom_node_point(aes(colour = group), size = 3) +
  geom_text(
    data = drug_nodes_ct,
    aes(x = x, y = y, label = name),
    size = 3,
    vjust = -0.4
  ) +
  scale_colour_manual(values = c("Drug" = "#1f77b4", "Disease" = "#d62728")) +
  scale_edge_width(range = c(0.2, 2)) +
  labs(title = "Drug–Disease Treatment Network (ClinicalTrials.gov)",
       colour = "Node Type",
       width  = "Number of Trials") +
  theme_void()

# ---- Drug–drug co-treatment network (CT.gov) ----
drug_disease_list_ct <- split(edges_ct, edges_ct$disease)

drug_drug_edges_ct <- map_dfr(drug_disease_list_ct, function(df) {
  if (nrow(df) < 2) return(NULL)
  combos <- t(combn(df$drug, 2))
  tibble(
    drug1 = combos[, 1],
    drug2 = combos[, 2],
    w1    = df$weight[match(combos[, 1], df$drug)],
    w2    = df$weight[match(combos[, 2], df$drug)],
    edge_weight = pmin(w1, w2)
  )
}) %>%
  group_by(drug1, drug2) %>%
  summarise(weight = sum(edge_weight), .groups = "drop")

g_dd_ct <- graph_from_data_frame(drug_drug_edges_ct, directed = FALSE)

dd_metrics_ct <- tibble(
  drug            = V(g_dd_ct)$name,
  weighted_degree = strength(g_dd_ct, weights = E(g_dd_ct)$weight)
)

top20_dd_ct <- dd_metrics_ct %>%
  slice_max(weighted_degree, n = 20)

ggplot(top20_dd_ct,
       aes(x = fct_reorder(drug, weighted_degree), y = weighted_degree)) +
  geom_col(fill = "#ff7f0e", alpha = .9) +
  coord_flip() +
  labs(title = "Top 20 Drugs in Drug–Drug Co-treatment Network (ClinicalTrials.gov)",
       x = "Drug",
       y = "Weighted Degree (Co-treatment Strength)") +
  theme_minimal(base_size = 13)

g_dd_ct_top <- induced_subgraph(g_dd_ct, vids = top20_dd_ct$drug)

set.seed(123)
layout_dd_ct <- create_layout(g_dd_ct_top, layout = "fr")

ggraph(layout_dd_ct) +
  geom_edge_link(aes(width = weight), colour = "grey80", alpha = 0.6) +
  geom_node_point(colour = "#ff7f0e", size = 3) +
  geom_text(aes(x = x, y = y, label = name), vjust = -0.4, size = 3) +
  scale_edge_width(range = c(0.2, 2)) +
  labs(
    title = "Drug–Drug Co-treatment Network (ClinicalTrials.gov, Top 20 Drugs)",
    edge_width = "Co-treatment Weight"
  ) +
  theme_void()



########################################################################
###############################################################
# 📘 ClinicalTrials.gov Drug–Disease + Drug–Drug Pipeline
# Option B – aligned with PubMed rentrez pipeline
# Author: Marvin Nukunu-Attachey
###############################################################

# ================== Libraries ==================
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(httr)
library(jsonlite)
library(igraph)
library(ggplot2)
library(forcats)
library(ggraph)
library(ggrepel)

# ================== Helper: Safe Extract ==================
get_or_na <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  as.character(x)
}

# ================== Fetch from ClinicalTrials.gov v2 ==================
fetch_ctgov <- function(drug, max_pages = 3, page_size = 100) {
  message("🔍 Fetching ClinicalTrials.gov data for: ", drug)
  
  base_url <- "https://clinicaltrials.gov/api/v2/studies"
  all_records <- list()
  next_page <- NULL
  
  for (i in seq_len(max_pages)) {
    query_params <- list(
      format    = "json",
      "query.intr" = drug,
      pageSize  = page_size
      # (no fields=… to avoid 400 errors; we parse from full payload)
    )
    if (!is.null(next_page)) {
      query_params$pageToken <- next_page
    }
    
    res <- tryCatch(GET(base_url, query = query_params), error = function(e) NULL)
    if (is.null(res) || httr::status_code(res) != 200) {
      message("⚠️ API error for ", drug,
              " (status ", ifelse(is.null(res), NA, httr::status_code(res)), ")")
      break
    }
    
    js <- fromJSON(content(res, as = "text", encoding = "UTF-8"),
                   simplifyDataFrame = FALSE)
    
    if (is.null(js$studies) || length(js$studies) == 0) {
      message("⚠️ No studies found for ", drug)
      break
    }
    
    df <- map_dfr(js$studies, function(s) {
      ps  <- s$protocolSection
      idm <- ps$identificationModule
      stm <- ps$statusModule
      cdm <- ps$conditionsModule
      dsm <- ps$designModule
      
      nct_id       <- get_or_na(idm$nctId)
      title        <- get_or_na(idm$briefTitle)
      status       <- get_or_na(stm$overallStatus)
      first_submit <- get_or_na(stm$studyFirstSubmitDate)
      
      phase <- if (!is.null(dsm$phases) && length(dsm$phases) > 0) {
        paste(dsm$phases, collapse = "; ")
      } else NA_character_
      
      disease <- if (!is.null(cdm$conditions) && length(cdm$conditions) > 0) {
        paste(cdm$conditions, collapse = "; ")
      } else NA_character_
      
      tibble(
        drug         = drug,
        nct_id       = nct_id,
        title        = title,
        status       = status,
        phase        = phase,
        disease      = disease,
        first_submit = first_submit
      )
    })
    
    all_records[[length(all_records) + 1]] <- df
    next_page <- js$nextPageToken
    if (is.null(next_page)) break
  }
  
  bind_rows(all_records)
}

# ================== Drug List (same as PubMed script) ==================
drugs <- c(
  "Atorvastatin", "Metformin", "Levothyroxine", "Lisinopril", "Amlodipine",
  "Metoprolol", "Albuterol", "Losartan", "Gabapentin", "Omeprazole",
  "Sertraline", "Rosuvastatin", "Pantoprazole", "Escitalopram",
  "Dextroamphetamine; Dextroamphetamine Saccharate; Amphetamine; Amphetamine Aspartate",
  "Hydrochlorothiazide", "Bupropion", "Fluoxetine", "Semaglutide", "Montelukast",
  "Trazodone", "Simvastatin", "Amoxicillin", "Tamsulosin",
  "Acetaminophen; Hydrocodone", "Fluticasone", "Meloxicam", "Apixaban",
  "Furosemide", "Insulin Glargine", "Duloxetine", "Ibuprofen", "Famotidine",
  "Empagliflozin", "Carvedilol", "Tramadol", "Alprazolam", "Prednisone",
  "Hydroxyzine", "Buspirone", "Clopidogrel", "Glipizide", "Citalopram",
  "Potassium Chloride", "Allopurinol", "Aspirin", "Cyclobenzaprine",
  "Ergocalciferol", "Oxycodone", "Methylphenidate", "Venlafaxine",
  "Spironolactone", "Ondansetron", "Zolpidem", "Cetirizine", "Estradiol",
  "Pravastatin", "Hydrochlorothiazide; Lisinopril", "Lamotrigine",
  "Quetiapine", "Fluticasone; Salmeterol", "Clonazepam", "Dulaglutide",
  "Azithromycin", "Hydrochlorothiazide; Losartan", "Amoxicillin; Clavulanate",
  "Latanoprost", "Cholecalciferol", "Propranolol", "Ezetimibe",
  "Topiramate", "Paroxetine", "Diclofenac", "Budesonide; Formoterol",
  "Atenolol", "Lisdexamfetamine", "Doxycycline", "Pregabalin",
  "Ethinyl Estradiol; Norethindrone", "Glimepiride", "Tizanidine",
  "Clonidine", "Fenofibrate", "Insulin Lispro", "Valsartan", "Cephalexin",
  "Baclofen", "Rivaroxaban", "Ferrous Sulfate", "Amitriptyline",
  "Finasteride", "Dapagliflozin", "Acetaminophen; Oxycodone", "Folic Acid",
  "Aripiprazole", "Olmesartan", "Ethinyl Estradiol; Norgestimate",
  "Valacyclovir", "Mirtazapine", "Lorazepam", "Levetiracetam",
  "Insulin Aspart", "Naproxen", "Cyanocobalamin", "Loratadine",
  "Diltiazem", "Sumatriptan", "Triamcinolone", "Hydralazine",
  "Tirzepatide", "Celecoxib", "Acetaminophen", "Alendronate",
  "Oxybutynin", "Hydrochlorothiazide; Triamterene", "Warfarin",
  "Progesterone", "Fluticasone; Umeclidinium; Vilanterol", "Testosterone",
  "Nifedipine", "Methocarbamol", "Benzonatate", "Sitagliptin",
  "Chlorthalidone", "Isosorbide", "Donepezil", "Dexmethylphenidate",
  "Sulfamethoxazole; Trimethoprim", "Clobetasol", "Methotrexate",
  "Hydroxychloroquine", "Lovastatin", "Pioglitazone", "Irbesartan",
  "Methylprednisolone", "Norethindrone", "Meclizine",
  "Ethinyl Estradiol; Levonorgestrel", "Fluticasone; Vilanterol",
  "Ketoconazole", "Thyroid", "Azelastine", "Nitrofurantoin",
  "Adalimumab", "Memantine", "Prednisolone", "Esomeprazole", "Docusate",
  "Clindamycin", "Acyclovir"
)

# ================== Fetch or Load ClinicalTrials Data ==================
if (file.exists("ctgov_results.csv")) {
  ct_df <- read.csv("ctgov_results.csv", stringsAsFactors = FALSE)
} else {
  ct_df <- map_df(drugs, fetch_ctgov)
  write.csv(ct_df, "ctgov_results.csv", row.names = FALSE)
}

# ================== Basic Cleaning & Long Format ==================
ct_long <- ct_df %>%
  filter(!is.na(disease), disease != "") %>%
  separate_rows(disease, sep = ";") %>%
  mutate(disease = str_trim(disease)) %>%
  filter(disease != "") %>%
  # drop generic / non-disease phrases
  filter(!str_detect(disease, regex("Healthy|Volunteer|Progression|Participant|Study", ignore_case = TRUE)))

# ================== Edge List: Drug–Disease (weight = #trials) ==================
edges_ct <- ct_long %>%
  count(drug, disease, name = "weight")

# ================== Summaries ==================
# Top diseases by total trial links
ct_top_diseases <- edges_ct %>%
  count(disease, wt = weight, name = "total_trials") %>%
  slice_max(total_trials, n = 20)

# Top drugs by total trial links
ct_top_drugs <- edges_ct %>%
  count(drug, wt = weight, name = "total_trials") %>%
  slice_max(total_trials, n = 20)

# Breadth vs total trials (per drug)
ct_drug_breadth <- edges_ct %>%
  group_by(drug) %>%
  summarise(
    unique_diseases = n_distinct(disease),
    total_trials    = sum(weight),
    .groups = "drop"
  )

# Disease co-annotation density: how many drugs per disease
ct_disease_density <- edges_ct %>%
  count(disease, name = "n_drugs") %>%
  mutate(
    category = ifelse(n_drugs > median(n_drugs), "high-density", "low-density")
  )

# ================== Yearly Trend (First Submit Date) ==================
ct_year_counts <- ct_df %>%
  filter(!is.na(first_submit), first_submit != "") %>%
  mutate(year = suppressWarnings(as.integer(substr(first_submit, 1, 4)))) %>%
  filter(!is.na(year)) %>%
  count(year)

# ================== Build Bipartite Graph ==================
g_ct <- graph_from_data_frame(edges_ct, directed = FALSE)
V(g_ct)$type   <- V(g_ct)$name %in% drugs
V(g_ct)$degree <- degree(g_ct, mode = "all")

# ================== Helper: Build Drug–Drug Co-Treatment Edges ==================
build_drug_drug_edges <- function(edges_df) {
  by_disease <- split(edges_df, edges_df$disease)
  
  dd_raw <- map_dfr(by_disease, function(df) {
    if (nrow(df) < 2) return(NULL)
    
    combos <- t(combn(df$drug, 2))
    tibble(
      drug1 = combos[, 1],
      drug2 = combos[, 2],
      w1    = df$weight[match(combos[, 1], df$drug)],
      w2    = df$weight[match(combos[, 2], df$drug)],
      edge_weight = pmin(w1, w2)
    )
  })
  
  dd_raw %>%
    group_by(drug1, drug2) %>%
    summarise(weight = sum(edge_weight), .groups = "drop")
}

dd_edges_ct <- build_drug_drug_edges(edges_ct)
g_dd_ct <- graph_from_data_frame(dd_edges_ct, directed = FALSE)

# ================== Drug-Level Metrics (Bipartite Graph) ==================
valid_ct_drugs <- intersect(drugs, V(g_ct)$name)

ct_drug_metrics <- tibble(
  drug            = valid_ct_drugs,
  degree          = degree(g_ct, v = valid_ct_drugs),
  betweenness     = betweenness(g_ct, v = valid_ct_drugs, directed = FALSE),
  weighted_degree = strength(g_ct, vids = valid_ct_drugs, weights = E(g_ct)$weight)
)

# ================== Visuals ==================

## 1. Top 15 diseases
ggplot(ct_top_diseases,
       aes(x = fct_reorder(disease, total_trials), y = total_trials)) +
  geom_col(fill = "#d62728", alpha = 0.8) +
  coord_flip() +
  labs(
    title = "Top 20 Diseases by Trial Count (ClinicalTrials.gov)",
    x     = "Disease",
    y     = "Total Trials"
  ) +
  theme_minimal(base_size = 13)

## 2. Top 15 drugs
ggplot(ct_top_drugs,
       aes(x = fct_reorder(drug, total_trials), y = total_trials)) +
  geom_col(fill = "#1f77b4", alpha = 0.8) +
  coord_flip() +
  labs(
    title = "Top 20 Drugs by Trial Count (ClinicalTrials.gov)",
    x     = "Drug",
    y     = "Total Trials"
  ) +
  theme_minimal(base_size = 13)

## 3. Breadth vs total trials – top 20 drugs only
ct_breadth_top20 <- ct_drug_breadth %>%
  slice_max(total_trials, n = 20)

ggplot(ct_breadth_top20,
       aes(x = unique_diseases, y = total_trials, label = drug)) +
  geom_point(aes(size = total_trials, colour = unique_diseases), alpha = 0.8) +
  geom_text(vjust = -0.6, size = 3) +
  scale_color_viridis_c() +
  labs(
    title = "Drug Breadth vs Total Trials (Top 20 Drugs, ClinicalTrials.gov)",
    x     = "Unique Diseases",
    y     = "Total Trials"
  ) +
  theme_minimal(base_size = 13)

## 4. Disease co-annotation density
ggplot(ct_disease_density,
       aes(x = n_drugs, fill = category)) +
  geom_histogram(bins = 20, alpha = 0.8, position = "identity") +
  scale_fill_manual(values = c("high-density" = "#d62728",
                               "low-density"  = "#1f77b4")) +
  labs(
    title = "Distribution of Number of Drugs per Disease (ClinicalTrials.gov)",
    x     = "Number of Drugs Treating Disease",
    y     = "Count of Diseases"
  ) +
  theme_minimal(base_size = 13)

## 5. Top 20 drugs by weighted evidence (bipartite)
ct_top20_weighted <- ct_drug_metrics %>%
  arrange(desc(weighted_degree)) %>%
  slice_head(n = 20)

ggplot(ct_top20_weighted,
       aes(x = fct_reorder(drug, weighted_degree), y = weighted_degree)) +
  geom_col(fill = "#1f77b4") +
  coord_flip() +
  labs(
    title = "Top 20 Drugs by Weighted Treatment Evidence (ClinicalTrials.gov)",
    x     = "Drug",
    y     = "Weighted Degree (Trials × Diseases)"
  ) +
  theme_minimal(base_size = 13)

## 6. Yearly trend: number of trials by first submit year
ggplot(ct_year_counts,
       aes(x = year, y = n)) +
  geom_col(fill = "#6baed6") +
  labs(
    title = "Yearly Trend of Registered Trials (ClinicalTrials.gov)",
    x     = "First Submit Year",
    y     = "Number of Trials"
  ) +
  theme_minimal(base_size = 13)

# ================== Networks ==================

## 7. Drug–Disease bipartite network (static, drug labels only)
set.seed(123)
layout_ct <- create_layout(g_ct, layout = "fr")

drug_nodes_ct <- as.data.frame(layout_ct) %>%
  filter(name %in% drugs)

ggraph(layout_ct) +
  geom_edge_link(aes(width = weight), colour = "grey80", alpha = 0.6) +
  geom_node_point(aes(colour = ifelse(name %in% drugs, "Drug", "Disease")),
                  size = 3.5) +
  geom_text_repel(
    data = drug_nodes_ct,
    aes(x = x, y = y, label = name, colour = "Drug"),
    size = 3, max.overlaps = Inf, box.padding = 0.25
  ) +
  scale_color_manual(values = c("Drug" = "#1f77b4", "Disease" = "#d62728")) +
  scale_edge_width(range = c(0.2, 2)) +
  labs(title = "Drug–Disease Treatment Network (ClinicalTrials.gov)",
       colour = "Node Type",
       width  = "Trials per Drug–Disease Pair") +
  theme_void()

## 8. Drug–Drug co-treatment network (top 20 drugs)
# Use the same top 20 by weighted_degree
top20_drugs_names <- ct_top20_weighted$drug

dd_edges_ct_top <- dd_edges_ct %>%
  filter(drug1 %in% top20_drugs_names,
         drug2 %in% top20_drugs_names)

g_dd_ct_top <- graph_from_data_frame(dd_edges_ct_top, directed = FALSE)

set.seed(123)
layout_dd_ct <- create_layout(g_dd_ct_top, layout = "fr")

ggraph(layout_dd_ct) +
  geom_edge_link(aes(width = weight), colour = "grey80", alpha = 0.6) +
  geom_node_point(colour = "#ff7f0e", size = 3) +
  geom_text(
    aes(x = x, y = y, label = name),
    vjust = -0.4,
    size = 3
  ) +
  scale_edge_width(range = c(0.2, 2)) +
  labs(
    title = "Drug–Drug Co-treatment Network (Top 20 Drugs, ClinicalTrials.gov)",
    subtitle = "Edge width reflects shared trial evidence"
  ) +
  theme_void()

ct_stats <- list(
  total_nodes = gorder(g_ct),
  total_edges = gsize(g_ct),
  drugs = sum(V(g_ct)$name %in% drugs),
  diseases = sum(!(V(g_ct)$name %in% drugs)),
  density = edge_density(g_ct),
  components = count_components(g_ct),
  giant_component = max(components(g_ct)$csize),
  avg_deg = mean(degree(g_ct)),
  median_deg = median(degree(g_ct)),
  avg_path = mean_distance(g_ct, directed = FALSE),
  diameter = diameter(g_ct, directed = FALSE)
)
ct_stats

ct_dd_stats <- list(
  total_nodes = gorder(g_dd_ct),
  total_edges = gsize(g_dd_ct),
  density = edge_density(g_dd_ct),
  components = count_components(g_dd_ct),
  giant_component = max(components(g_dd_ct)$csize),
  avg_deg = mean(degree(g_dd_ct)),
  diameter = diameter(g_dd_ct)
)
ct_dd_stats


