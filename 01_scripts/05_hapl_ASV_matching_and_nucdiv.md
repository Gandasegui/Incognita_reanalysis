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
This script needs the following files, adapted from the Supplementary Data 3 in https://www.nature.com/articles/s41467-025-64516-6:
- cote_sequence.txt
- laos_sequence.txt
- anzania_sequence.txt
- uganda_sequence.txt

```R
library(tidyverse)
library(Biostrings)
library(pwalign)


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

bootstrap_pi <- function(dist_mat, freqs, n_boot = 1000) {
  n <- length(freqs)
  boot_pis <- numeric(n_boot)
  for (b in 1:n_boot) {
    idx    <- sample(1:n, n, replace = TRUE)
    d_boot <- dist_mat[idx, idx]
    f_boot <- freqs[idx] / sum(freqs[idx])
    fm     <- outer(f_boot, f_boot)
    boot_pis[b] <- sum(fm[upper.tri(fm)] * d_boot[upper.tri(d_boot)]) * 2
  }
  boot_pis
}

calc_pi_from_dist <- function(dist_mat, freqs) {
  fm <- outer(freqs, freqs)
  sum(fm[upper.tri(fm)] * dist_mat[upper.tri(dist_mat)]) * 2
}


# -----------------------------------------------------------------------------
# Core function
# -----------------------------------------------------------------------------

calculate_nucleotide_diversity <- function(input_file, n_boot = 1000) {
  
  library(Biostrings)
  library(pwalign)
  library(tidyverse)
  
  # -------------------------
  # 1. Load sequences
  # -------------------------
  df <- read_tsv(input_file, col_types = cols())
  df <- df %>% mutate(uid = paste0("seq_", row_number()))
  
  # -------------------------
  # 2. Compute FULL pairwise distance matrix once
  # -------------------------
  message("Computing full pairwise distance matrix...")
  all_seqs  <- DNAStringSet(df$Sequence)
  names(all_seqs) <- df$uid
  full_dist <- pwalign::stringDist(all_seqs, method = "levenshtein") %>% as.matrix()
  full_dist <- full_dist / mean(width(all_seqs))
  
  # -------------------------
  # 3. Calculate π per species
  # -------------------------
  species_list <- unique(df$Species_ID)
  
  results <- map_dfr(species_list, function(sp) {
    
    sp_df  <- df %>% filter(Species_ID == sp)
    sp_idx <- which(df$Species_ID == sp)
    
    if (length(sp_idx) < 2) {
      message("Skipping ", sp, ": only 1 sequence")
      return(tibble(Species_ID                  = sp,
                    pi_unweighted               = NA, pi_unweighted_ci_low  = NA, pi_unweighted_ci_high = NA,
                    pi_coverage_weighted        = NA, pi_weighted_ci_low    = NA, pi_weighted_ci_high   = NA,
                    pi_weighted_vs_unweighted_p = NA,
                    n_sequences                 = 1,
                    total_coverage              = sum(sp_df$Coverage)))
    }
    
    message("Processing: ", sp, " (", length(sp_idx), " sequences)")
    
    dist_mat <- full_dist[sp_idx, sp_idx]
    
    freq_uw <- rep(1 / length(sp_idx), length(sp_idx))
    freq_w  <- sp_df$Coverage / sum(sp_df$Coverage)
    
    pi_unweighted <- calc_pi_from_dist(dist_mat, freq_uw)
    pi_weighted   <- calc_pi_from_dist(dist_mat, freq_w)
    
    # Bootstrap CIs
    boot_uw <- bootstrap_pi(dist_mat, freq_uw, n_boot)
    boot_w  <- bootstrap_pi(dist_mat, freq_w,  n_boot)
    ci_uw   <- quantile(boot_uw, c(0.025, 0.975))
    ci_w    <- quantile(boot_w,  c(0.025, 0.975))
    
    # Within-species weighted vs unweighted permutation test
    obs_diff   <- pi_weighted - pi_unweighted
    perm_diffs <- replicate(n_boot, {
      perm_freq <- sample(freq_w)
      calc_pi_from_dist(dist_mat, perm_freq / sum(perm_freq)) - pi_unweighted
    })
    p_val <- mean(abs(perm_diffs) >= abs(obs_diff))
    
    tibble(Species_ID                  = sp,
           pi_unweighted               = round(pi_unweighted, 6),
           pi_unweighted_ci_low        = round(ci_uw[1], 6),
           pi_unweighted_ci_high       = round(ci_uw[2], 6),
           pi_coverage_weighted        = round(pi_weighted, 6),
           pi_weighted_ci_low          = round(ci_w[1], 6),
           pi_weighted_ci_high         = round(ci_w[2], 6),
           pi_weighted_vs_unweighted_p = round(p_val, 4),
           n_sequences                 = length(sp_idx),
           total_coverage              = sum(sp_df$Coverage))
  })
  
  # -------------------------
  # 4. Between-species permutation tests
  # -------------------------
  if (length(species_list) > 1) {
    message("\nRunning pairwise between-species permutation tests...")
    sp_pairs <- combn(species_list, 2, simplify = FALSE)
    
    between_tests <- map_dfr(sp_pairs, function(pair) {
      
      idx1 <- which(df$Species_ID == pair[1])
      idx2 <- which(df$Species_ID == pair[2])
      
      if (length(idx1) < 2 | length(idx2) < 2) {
        return(tibble(species_1 = pair[1], species_2 = pair[2],
                      p_value_unweighted = NA, p_value_weighted = NA))
      }
      
      pool_idx <- c(idx1, idx2)
      n1       <- length(idx1)
      n2       <- length(idx2)
      
      # Submatrix for pooled sequences - computed once
      pool_dist  <- full_dist[pool_idx, pool_idx]
      pool_cov   <- df$Coverage[pool_idx]
      
      # Observed π difference
      obs_uw <- abs(calc_pi_from_dist(pool_dist[1:n1, 1:n1],         rep(1/n1, n1)) -
                      calc_pi_from_dist(pool_dist[(n1+1):(n1+n2), (n1+1):(n1+n2)], rep(1/n2, n2)))
      
      obs_w  <- abs(calc_pi_from_dist(pool_dist[1:n1, 1:n1],
                                      pool_cov[1:n1] / sum(pool_cov[1:n1])) -
                      calc_pi_from_dist(pool_dist[(n1+1):(n1+n2), (n1+1):(n1+n2)],
                                        pool_cov[(n1+1):(n1+n2)] / sum(pool_cov[(n1+1):(n1+n2)])))
      
      # Permutation: shuffle labels, reindex precomputed distance matrix
      perm_results <- replicate(n_boot, {
        idx      <- sample(n1 + n2)
        pi1_idx  <- idx[1:n1]
        pi2_idx  <- idx[(n1+1):(n1+n2)]
        
        d1 <- pool_dist[pi1_idx, pi1_idx]
        d2 <- pool_dist[pi2_idx, pi2_idx]
        c1 <- pool_cov[pi1_idx]
        c2 <- pool_cov[pi2_idx]
        
        perm_uw <- abs(calc_pi_from_dist(d1, rep(1/n1, n1)) -
                         calc_pi_from_dist(d2, rep(1/n2, n2)))
        perm_w  <- abs(calc_pi_from_dist(d1, c1 / sum(c1)) -
                         calc_pi_from_dist(d2, c2 / sum(c2)))
        c(perm_uw, perm_w)
      })
      
      tibble(species_1          = pair[1],
             species_2          = pair[2],
             p_value_unweighted = round(mean(perm_results[1, ] >= obs_uw), 4),
             p_value_weighted   = round(mean(perm_results[2, ] >= obs_w),  4))
    })
    
    print(between_tests)
    base <- tools::file_path_sans_ext(basename(input_file))
    write_tsv(between_tests, paste0(base, "_between_species_tests.tsv"))
    message("Between-species tests saved to: ", paste0(base, "_between_species_tests.tsv"))
  }
  
  # -------------------------
  # 5. Save and return results
  # -------------------------
  base     <- tools::file_path_sans_ext(basename(input_file))
  out_file <- paste0(base, "_nucleotide_diversity.tsv")
  write_tsv(results, out_file)
  message("Results saved to: ", out_file)
  
  print(results)
  return(invisible(list(diversity = results)))
}





run_calculate_nucleotide_diversity <- function(input_files, n_boot = 1000) {
  
  all_diversity     <- list()
  all_between_tests <- list()
  
  for (f in input_files) {
    message("\n========== Processing: ", f, " ==========")
    result <- calculate_nucleotide_diversity(f, n_boot = n_boot)
    
    # Tag results with source file
    base <- tools::file_path_sans_ext(basename(f))
    all_diversity[[base]] <- result$diversity %>% mutate(file = base)
    
    # Load between-species tests if they were saved
    between_file <- paste0(base, "_between_species_tests.tsv")
    if (file.exists(between_file)) {
      all_between_tests[[base]] <- read_tsv(between_file, col_types = cols()) %>%
        mutate(file = base)
    }
  }
  
  # Combine and save summary tables
  combined_diversity     <- bind_rows(all_diversity)
  combined_between_tests <- bind_rows(all_between_tests)
  
  write_tsv(combined_diversity,     "combined_nucleotide_diversity.tsv")
  write_tsv(combined_between_tests, "combined_between_species_tests.tsv")
  
  message("\nCombined diversity saved to: combined_nucleotide_diversity.tsv")
  message("Combined between-species tests saved to: combined_between_species_tests.tsv")
  
  print(combined_diversity)
  print(combined_between_tests)
  
  return(invisible(list(diversity = combined_diversity,
                        between_tests = combined_between_tests)))
}

# Run all files
results <- run_calculate_nucleotide_diversity(
  input_files = c("data/cote_sequence.txt",
                  "data/laos_sequence.txt",
                  "data/tanzania_sequence.txt",
                  "data/uganda_sequence.txt"),
  n_boot = 1000
)




plot_nucleotide_diversity <- function(combined_file = "combined_nucleotide_diversity.tsv",
                                      output_file   = "nucleotide_diversity_plot.pdf",
                                      plot_width    = 12,
                                      plot_height   = 6) {
  library(tidyverse)
  
  df <- read_tsv(combined_file, col_types = cols())
  
  # Reshape to long format for faceting
  long_df <- bind_rows(
    df %>% transmute(Species_ID,
                     file,
                     pi       = pi_unweighted,
                     ci_low   = pi_unweighted_ci_low,
                     ci_high  = pi_unweighted_ci_high,
                     type     = "Unweighted"),
    df %>% transmute(Species_ID,
                     file,
                     pi       = pi_coverage_weighted,
                     ci_low   = pi_weighted_ci_low,
                     ci_high  = pi_weighted_ci_high,
                     type     = "ASV count weighted")
  ) %>%
    mutate(type = factor(type, levels = c("Unweighted", "ASV count weighted")))
  long_df$Species_ID <- str_replace(long_df$Species_ID, 'Trichuris_incognita', 'T. incognita')
  long_df$Species_ID <- str_replace(long_df$Species_ID, 'Trichuris_trichiura', 'T. trchiura')
  p <- ggplot(long_df, aes(x = Species_ID, y = pi, color = file)) +
    geom_pointrange(aes(ymin = ci_low, ymax = ci_high),
                    position = position_dodge(width = 0.5),
                    size = 0.6) +
    facet_wrap(~ type, ncol = 2) +
    scale_color_lancet() +
    labs(y = expression("Nucleotide diversity (" ~ pi ~ ")")) +
    theme_classic() +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          strip.text    = element_text(face = "bold"),
          legend.position = "none")
  
  p
}

plot_nucleotide_diversity()
```
plot_nucleotide_diversity() generates the Fig 1d
