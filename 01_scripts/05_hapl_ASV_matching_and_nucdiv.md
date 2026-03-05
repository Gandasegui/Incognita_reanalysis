# Comparing ASVs vs haplotype counts and estimaitng Pi (nucleotide diversity)

## This first part links the ASVs with haplotypes with 100% identity and compare their frequencies

### The following files are needed:
- run_haplotype_counts.with_species.sorted.tsv; the haplotyope count and species generated before
- hap_vs_asv.with_header.b6; ASVs and haplotypes with 100% identity
- ASVs_clean_samples.csv; ASV count per sample adapted form the Supplementary Data 3 in https://www.nature.com/articles/s41467-025-64516-6
- ASVs_clean_species.xlsx; ASV species assgiment adapted form the Supplementary Data 3 in https://www.nature.com/articles/s41467-025-64516-6

```R
library(tidyverse)
library(scales)
library(readxl)
library(ggh4x)

# -----------------------------
# Input files
# -----------------------------
setwd("~/R/ITS_trichuris")
hap_counts_fp <- "data/run_haplotype_counts.with_species.polished.sorted.tsv"
hap_vs_asv_fp <- "data/hap_vs_asv.with_header.b6"
asv_tab_fp    <- "data/ASVs_clean/ASVs_clean_samples.csv"
asv_sp_fp     <- "data/ASVs_clean/ASVs_clean_species.xlsx"

hap_counts <- read_tsv(hap_counts_fp, show_col_types = FALSE)
hap_vs_asv <- read_tsv(hap_vs_asv_fp, show_col_types = FALSE)
asv_tab <- read_csv(asv_tab_fp, show_col_types = FALSE)
asv_sp <- read_xlsx(asv_sp_fp, col_names = FALSE)

# 1) Haplotype total reads
hap_totals <- hap_counts %>%
  group_by(Haplotype) %>%
  summarise(hap_reads_total = sum(Reads), .groups = "drop")

# 2) Map haplotypes to ASVs (100% matches)
hap_asv_map <- hap_vs_asv %>%
  transmute(Haplotype = qseqid, ASV_raw = sseqid) %>%
  mutate(ASV_ID = sub("_[^_]*$", "", ASV_raw))

hap_asv_map2 <- hap_asv_map %>%
  left_join(hap_totals, by = "Haplotype")

# Max signal per ASV
asv_hap_signal <- hap_asv_map2 %>%
  group_by(ASV_ID) %>%
  summarise(hap_reads_total = max(hap_reads_total, na.rm = TRUE), .groups = "drop")

# 3) ASV total reads (Rahman)
asv_cols <- names(asv_tab)[str_detect(names(asv_tab), "^ASV_Clus_")]
asv_totals <- asv_tab %>%
  select(all_of(asv_cols)) %>%
  summarise(across(everything(), ~ sum(as.numeric(.x), na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "ASV_ID", values_to = "asv_reads_total")

# ASV -> species map
asv_map <- asv_sp %>%
  rename(key = 1) %>% mutate(key = str_replace_all(key, "\\s+", "_")) %>%
  pivot_longer(-key, names_to = "col", values_to = "value") %>%
  pivot_wider(names_from = key, values_from = value) %>%
  rename(ASV_ID = Country, Species_raw = SpeciesID) %>%
  mutate(Species = case_when(
    str_detect(Species_raw, regex("incognita", ignore_case = TRUE)) ~ "incognita",
    str_detect(Species_raw, regex("trichiura", ignore_case = TRUE)) ~ "trichiura",
    TRUE ~ "other"
  )) %>%
  select(ASV_ID, Species)

# 4) Merge for plotting
plot_df <- asv_totals %>%
  left_join(asv_map, by = "ASV_ID") %>%
  left_join(asv_hap_signal, by = "ASV_ID") %>%
  mutate(Species = factor(Species, levels = c("incognita", "trichiura", "other"))) %>%
  arrange(desc(asv_reads_total)) %>%
  mutate(ASV_ID = factor(ASV_ID, levels = ASV_ID)) %>%
  filter(Species != 'other') %>%
  mutate(scale_group = ifelse(asv_reads_total > 500000, "High (>500k)", "Low (<500k)")) %>%
  # Ensure High is on top and Low is on bottom
  mutate(scale_group = factor(scale_group, levels = c("High (>500k)", "Low (<500k)")))

# 1. Define your custom species labels
species_names <- c("incognita" = "T. incognita", 
                   "trichiura" = "T. trichiura")

# 2. Plotting
p_asv_haplink_split <- ggplot(plot_df, 
                              aes(x = ASV_ID, y = asv_reads_total, color = hap_reads_total, shape = Species)) +
  geom_point(size = 3, alpha = 0.9) +
  facet_grid(scale_group ~ Species, 
             scales = "free", 
             space = "free_x",
             labeller = labeller(Species = species_names)) + 
  force_panelsizes(rows = c(1, 0.6)) +
  scale_y_continuous("ASV count", labels = scales::comma) +
  scale_x_discrete("ASV (sorted by ASV count)") +
  scale_color_viridis_c(name = "Haplotype count", option = "H", na.value = "grey75") + theme_classic() +
  scale_shape(guide = "none") + 
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    panel.border = element_rect(color = "black", fill = NA)
  )

p_asv_haplink_split
```

p_asv_haplink_split corresponds to the figure 1c

## This second part describes how to estimate nucleotide diversity accountung for ASV count

```R

```
