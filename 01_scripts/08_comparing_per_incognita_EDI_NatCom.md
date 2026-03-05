# Reading ASV and haplotype counts and generating the plot

### For this code we need:
- ASV_FWD_vs_ITS1refs.min10.b6 and ASV_FWD_vs_ITS2refs.min10.b6; which species assgimend of EDI ASVs
- The ASV tables for DADA two obtiane in the previous script
- The metadta of the SRA archive of EDI paper (PRJNA1131306.RunInfo.csv)
- run_haplotype_counts.with_species.sorted.tsv; the haplotype count obtained in "04_generating_haplotype_counts.md"
- Nautre_paper_metadata.csv; extracted from the Supplementary Data 6 from https://www.nature.com/articles/s41467-025-64516-6
- metadata_fixed.wide.nospace.tsv; genertae in "/02_getting_NatCom_paper_fastqs_and_metadata.md"


```R
library(tidyverse)
library(readxl)
library(ggpubr)
library(ggsci)
library(dada2)

# Load EDI ASV classification results
its1_fp <- "data/ASV_FWD_vs_ITS1refs.min10.b6"
its2_fp <- "data/ASV_FWD_vs_ITS2refs.min10.b6"
asv1_tab <- "dada2_FWD_only_shared40/marker_ITS-1/ASV_table.tsv"
asv2_tab <- "dada2_FWD_only_shared40/marker_ITS-2/ASV_table.tsv"
edi_meta_fp <- "data/PRJNA1131306.RunInfo.csv"

read_besthit <- function(fp, cov_min = 0.7, pid_min = 90) {
  read_tsv(fp, col_names = FALSE, show_col_types = FALSE) %>%
    setNames(c("qseqid","sseqid","pident","length","qlen","slen",
               "qstart","qend","sstart","send","evalue","bitscore")) %>%
    mutate(cov = length / qlen) %>%
    filter(pident >= pid_min, cov >= cov_min) %>%
    group_by(qseqid) %>%
    slice_max(bitscore, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(ASV_ID = qseqid, Species = case_when(
      str_detect(sseqid, "incognita") ~ "incognita",
      str_detect(sseqid, "trichiura") ~ "trichiura",
      TRUE ~ "other"
    ))
}

asv_to_pct <- function(tab_fp, blast_map) {
  read_tsv(tab_fp, show_col_types = FALSE) %>%
    pivot_longer(-Run, names_to = "ASV_ID", values_to = "Reads") %>%
    filter(Reads > 0) %>%
    left_join(blast_map, by = "ASV_ID") %>%
    filter(Species %in% c("incognita", "trichiura")) %>%
    group_by(Run, Species) %>%
    summarise(Reads = sum(Reads), .groups = "drop") %>%
    pivot_wider(names_from = Species, values_from = Reads, values_fill = 0) %>%
    mutate(total = incognita + trichiura,
           incognita_pct_edi = if_else(total > 0, 100 * incognita / total, NA_real_))
}

# Read and compute % T. incognita for EDI
b1 <- read_besthit(its1_fp)
b2 <- read_besthit(its2_fp)
e1 <- asv_to_pct(asv1_tab, b1) %>% mutate(Marker = "ITS-1")
e2 <- asv_to_pct(asv2_tab, b2) %>% mutate(Marker = "ITS-2")
edi_asv <- bind_rows(e1, e2)

edi_meta <- read_csv(edi_meta_fp, show_col_types = FALSE) %>%
  select(Run, host_subject_id, geo_loc_name)

edi_ref <- edi_asv %>%
  left_join(edi_meta, by = "Run") %>%
  filter(!is.na(incognita_pct_edi)) %>%
  rename(Country = geo_loc_name)

# NAT haplotypes %
nat_hap_fp <- "data/run_haplotype_counts.with_species.sorted.tsv"
nat_meta_fp <- "data/nat_metadata.csv"
nat_fixed_fp <- "data/metadata/metadata_fixed.wide.nospace.tsv"

nat_hap <- read_tsv(nat_hap_fp, show_col_types = FALSE)
nat_meta <- read_csv(nat_meta_fp, show_col_types = FALSE)
nat_fixed <- read_tsv(nat_fixed_fp, col_names = c("sra_id", "biosample_id", "sample_id", "fastq_path", "country"), show_col_types = FALSE)

nat_pct <- nat_hap %>%
  group_by(Run, Species) %>%
  summarise(Reads = sum(Reads), .groups = "drop") %>%
  pivot_wider(names_from = Species, values_from = Reads, values_fill = 0) %>%
  transmute(Run, incognita = Trichuris_incognita_ITS2, trichiura = Trichuris_trichiura_ITS2) %>%
  mutate(total = incognita + trichiura,
         incognita_pct_nat = if_else(total > 0, 100 * incognita / total, NA_real_))

nat_ref <- nat_meta %>%
  select(Participant_ID, Sample_No) %>%
  left_join(nat_fixed, by = c("Sample_No" = "sample_id")) %>%
  left_join(nat_pct, by = c("sra_id" = "Run"))

# NAT ASV %
asv_clean_fp <- "data/ASVs_clean/ASVs_clean_samples.csv"
asv_map_fp <- "data/ASVs_clean/ASVs_clean_species.xlsx"
asv_samples <- read_csv(asv_clean_fp, show_col_types = FALSE)
asv_map_raw <- read_xlsx(asv_map_fp, col_names = FALSE)

asv_map <- asv_map_raw %>%
  rename(key = 1) %>%
  mutate(key = str_replace_all(key, "\\s+", "_")) %>%
  pivot_longer(-key, names_to = "col", values_to = "value") %>%
  pivot_wider(names_from = key, values_from = value) %>%
  rename(ASV_ID = Country, Species_raw = SpeciesID) %>%
  mutate(Species = case_when(
    str_detect(Species_raw, "incognita") ~ "incognita",
    str_detect(Species_raw, "trichiura") ~ "trichiura",
    TRUE ~ "other"
  )) %>%
  select(ASV_ID, Species)

asv_cols <- names(asv_samples)[str_detect(names(asv_samples), "^ASV_Clus_")]
asv_long <- asv_samples %>%
  pivot_longer(cols = all_of(asv_cols), names_to = "ASV_ID", values_to = "Reads") %>%
  left_join(asv_map, by = "ASV_ID") %>%
  filter(Species %in% c("incognita", "trichiura")) %>%
  group_by(sample, Species) %>%
  summarise(Reads = sum(Reads), .groups = "drop") %>%
  pivot_wider(names_from = Species, values_from = Reads, values_fill = 0) %>%
  mutate(total = incognita + trichiura,
         incognita_pct_asv = if_else(total > 0, 100 * incognita / total, NA_real_)) %>%
  left_join(select(nat_meta, Sample_No, Participant_ID), by = c("sample" = "Sample_No"))

nat_asv <- asv_long %>%
  group_by(Participant_ID) %>%
  summarise(incognita_pct_asv = mean(incognita_pct_asv, na.rm = TRUE), .groups = "drop")

# Join NAT + EDI
nat_edi <- nat_ref %>%
  left_join(edi_ref, by = c("Participant_ID" = "host_subject_id")) %>%
  left_join(nat_asv, by = "Participant_ID") %>%
  drop_na(incognita_pct_edi, incognita_pct_asv, incognita_pct_nat)

# p1: ITS-1 vs ITS-2 from EDI
p1 <- nat_edi %>%
  filter(Marker %in% c("ITS-1", "ITS-2")) %>%
  distinct(Participant_ID, Marker, .keep_all = TRUE) %>%
  group_by(Participant_ID) %>%
  filter(n() == 2) %>%
  ggplot(aes(x = Marker, y = incognita_pct_edi, color = Country, group = Participant_ID)) +
  geom_line(alpha = 0.4) + geom_point(size = 2) +
  scale_y_continuous("% T. incognita", limits = c(0, 100)) +
  scale_x_discrete(labels = c("ITS-1" = "Venkatesan et al., ITS-1", "ITS-2" = "Venkatesan et al., ITS-2")) +
  scale_color_lancet() +
  theme_classic() + theme(axis.title.x = element_blank(),
                            axis.text.x = element_text(angle = 25, hjust = 1))

# p2: ITS-1 vs Rahman (ASV/Haplotype)
p2 <- nat_edi %>%
  filter(Marker == "ITS-1") %>%
  select(Participant_ID, Country, incognita_pct_edi, incognita_pct_asv, incognita_pct_nat) %>%
  pivot_longer(cols = starts_with("incognita_pct"), names_to = "source", values_to = "value") %>%
  mutate(source = factor(source, levels = c("incognita_pct_edi", "incognita_pct_asv", "incognita_pct_nat"),
                         labels = c("Venkatesan et al., ITS-1", "Rahman et al., ITS-2 (ASVs)", "Rahman et al., ITS-2 (Haplotypes)"))) %>%
  ggplot(aes(x = source, y = value, color = Country, group = Participant_ID)) +
  geom_line(alpha = 0.4) + geom_point(size = 2) +
  scale_y_continuous("% T. incognita", limits = c(0, 100)) +
  scale_color_lancet() +
  theme_classic() + theme(axis.title.x = element_blank(),
                          axis.text.x = element_text(angle = 25, hjust = 1))

# p3: ITS-2 vs Rahman (ASV/Haplotype)
p3 <- nat_edi %>%
  filter(Marker == "ITS-2") %>%
  select(Participant_ID, Country, incognita_pct_edi, incognita_pct_asv, incognita_pct_nat) %>%
  pivot_longer(cols = starts_with("incognita_pct"), names_to = "source", values_to = "value") %>%
  mutate(source = factor(source, levels = c("incognita_pct_edi", "incognita_pct_asv", "incognita_pct_nat"),
                         labels = c("Venkatesan et al., ITS-2", "Rahman et al., ITS-2 (ASVs)", "Rahman et al., ITS-2 (Haplotypes)"))) %>%
  ggplot(aes(x = source, y = value, color = Country, group = Participant_ID)) +
  geom_line(alpha = 0.4) + geom_point(size = 2) +
  scale_y_continuous("% T. incognita", limits = c(0, 100)) +
  scale_color_lancet() +
  theme_classic() + theme(axis.title.x = element_blank(),
                          axis.text.x = element_text(angle = 25, hjust = 1))

# Final panel
per_inc <- ggarrange(p1, p2, p3, nrow = 1, widths = c(1, 1.2, 1.2), common.legend = TRUE, legend = "right")
```

per_inc corresponds to the Figure 1e
