# Analising the ASVs from the supplementary material provided in the NatCom paper

### For this analysis the next files are needed:
- ASVs_clean_samples.csv; adapted form the Supplementary Data 3 in https://www.nature.com/articles/s41467-025-64516-6
- metadata_fixed.wide.nospace.tsv; metadata genereated from SRA archive in the preivous script "02_getting_NatCom_paper_fastqs_and_metadata.md"

```R
library(tidyverse)
library(readxl)
library(scales)
library(ggpubr)
library(ggsci)

# -----------------------------
# Input files
# -----------------------------
setwd("~/R/ITS_trichuris")
asv_fp <- "data/ASVs_clean/ASVs_clean_samples.csv"
meta_fp <- "data/metadata/metadata_fixed.wide.nospace.tsv"

# Read ASV data and depth
asv_tab <- read_csv(asv_fp, show_col_types = FALSE)
meta <- read_tsv(meta_fp, 
                 col_names = c("sra_id", "biosample_id", "sample_id", "fastq_path", "Country"),
                 show_col_types = FALSE)

# Get ASV columns
asv_cols <- grep("^ASV_Clus_", colnames(asv_tab), value = TRUE)

# ASV -> species map from species Excel
asv_sp <- read_xlsx("data/ASVs_clean/ASVs_clean_species.xlsx", col_names = FALSE)
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

# Compute % incognita
asv_long <- asv_tab %>%
  pivot_longer(all_of(asv_cols), names_to = "ASV_ID", values_to = "Reads") %>%
  mutate(Reads = as.numeric(Reads)) %>% filter(!is.na(Reads), Reads > 0) %>%
  left_join(asv_map, by = "ASV_ID")

inc_df <- asv_long %>%
  filter(Species %in% c("incognita", "trichiura")) %>%
  group_by(sample) %>%
  summarise(
    incognita_reads = sum(Reads[Species == "incognita"], na.rm = TRUE),
    trichiura_reads = sum(Reads[Species == "trichiura"], na.rm = TRUE),
    total_trichuris_reads = incognita_reads + trichiura_reads,
    incognita_pct = if_else(total_trichuris_reads > 0,
                            100 * incognita_reads / total_trichuris_reads,
                            NA_real_),
    .groups = "drop"
  )

# Merge with metadata and ONT depth
plot_df <- inc_df %>%
  left_join(select(meta, sample_id, sra_id, Country), by = c("sample" = "sample_id"))

# Add depth group
plot_df <- plot_df %>%
  mutate(depth_group = if_else(total_trichuris_reads >= 10000, "ASV count > 10k", "ASV count < 10k"))

# Boxplot + stat comparison
# Adding statistics
# Wilcoxon test per country:
# Compare incognita % between depth groups within each country
stats_df <- plot_df %>%
  na.omit() %>%
  group_by(Country) %>%
  summarise(
    p_value = wilcox.test(incognita_pct ~ depth_group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_label = paste0("p = ", signif(p_value, 3)))


# Choose a y-position per facet (country) for placing p-value text:
# Start from 95th percentile, then (in your final workflow) overwrite with manual values
p_text_df <- plot_df %>%
  na.omit() %>%
  group_by(Country) %>%
  summarise(
    y_pos = quantile(incognita_pct, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(stats_df, by = "Country")

# Manual overrides (as per the final figure tuning)
p_text_df$y_pos <- c(15, 75, 75, 75)
p_text_df$p_label <- str_remove_all(p_text_df$p_label, ' ')

p_box <- ggplot(na.omit(plot_df), aes(x = depth_group, y = incognita_pct, color = Country)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  facet_wrap(~Country, nrow = 1) +
  scale_y_continuous("% T. incognita") +
  scale_color_lancet() +
  theme(legend.position = 'none') +
  theme_classic() + theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)
  ) +
  geom_text(
    data = p_text_df,
    aes(
      x = 1.75,     # centered between the two depth groups
      y = y_pos,
      label = p_label
    ),
    inherit.aes = FALSE,
    size = 3,
    fontface = "italic"
  )

p_box

# Log scatter
p_log <- ggplot(na.omit(plot_df),
                aes(x = total_trichuris_reads, y = incognita_pct, color = Country)) +
  geom_point(alpha = 0.8) +
  geom_vline(
    xintercept = 10000, linewidth = 1,
    color = "grey60", linetype = "dashed") +
  scale_x_log10("ASV count per sample (log10)", labels = comma) +
  scale_y_continuous("% T. incognita") +
  facet_wrap(~Country, scales = "free_x", ncol = 1) +
  scale_color_lancet() +
  theme_classic() +
  theme(legend.position = 'none')

p_log
```
p_log corresponds to the Figure 1a

p_box corresponds to the Figure 1b
