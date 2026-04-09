# =========================================
# R Script: Run permutation testing on ASV vs depth
# =========================================



library(tidyverse)

# -------------------------
# 1. Load data
# -------------------------
df <- read_csv("ASV_per_inc.csv", col_types = cols())




# -------------------------
# 2. Permutation test
# -------------------------
run_permutation_test <- function(data, n_perm = 10000) {
  
  # Observed correlation per country
  observed <- data %>%
    group_by(Country) %>%
    summarise(
      obs_cor   = cor(.data[["ASV_count"]], .data[["incognita_pct"]], method = "spearman"),
      n_samples = n(),
      .groups   = "drop"
    )
  
  # Permutation: shuffle ASV_count within each country
  perm_results <- map_dfr(1:n_perm, function(i) {
    data %>%
      group_by(Country) %>%
      mutate(ASV_count_perm = sample(.data[["ASV_count"]])) %>%
      summarise(
        perm_cor = cor(.data[["ASV_count_perm"]], .data[["incognita_pct"]], method = "spearman"),
        .groups  = "drop"
      ) %>%
      mutate(perm = i)
  })
  
  # P-value: proportion of permuted correlations >= observed
  pval_df <- observed %>%
    left_join(
      perm_results %>%
        group_by(Country) %>%
        summarise(
          mean_perm_cor = mean(perm_cor),
          sd_perm_cor   = sd(perm_cor),
          .groups       = "drop"
        ),
      by = "Country"
    ) %>%
    left_join(
      map_dfr(unique(data$Country), function(ctry) {
        obs  <- observed %>% filter(Country == ctry) %>% pull(obs_cor)
        null <- perm_results %>% filter(Country == ctry) %>% pull(perm_cor)
        tibble(
          Country = ctry,
          p_value = mean(abs(null) >= abs(obs))
        )
      }),
      by = "Country"
    ) %>%
    mutate(
      p_adj  = p.adjust(p_value, method = "BH"),
      signif = case_when(p_adj < 0.001 ~ "***",
                         p_adj < 0.01  ~ "**",
                         p_adj < 0.05  ~ "*",
                         TRUE          ~ "ns")
    )
  
  list(observed = observed, perm_results = perm_results, stats = pval_df)
}

results <- run_permutation_test(df, n_perm = 10000)

# -------------------------
# 3. Print results table
# -------------------------
message("\n--- Permutation test results ---")
print(
  results$stats %>%
    select(Country, n_samples, obs_cor, mean_perm_cor, p_value, p_adj, signif) %>%
    mutate(
      obs_cor       = round(obs_cor, 4),
      mean_perm_cor = round(mean_perm_cor, 4)
    )
)

# -------------------------
# 4. Plot: observed correlation with raw data
# -------------------------
p1 <- df %>%
  ggplot(aes(x = ASV_count, y = incognita_pct, color = Country)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  geom_text(data = results$stats,
            aes(x     = Inf,
                y     = Inf,
                label = paste0("r=",       round(obs_cor, 2),
                               "\np=",     p_value,
                               "\np.adj=", p_adj)),
            hjust = 1.1, vjust = 1.5, size = 3.5,
            inherit.aes = FALSE, color = "black") +
  facet_wrap(~ Country, scales = "free") +
  scale_color_brewer(palette = "Set1") +
  labs(x        = "ASV count",
       y        = "% Trichuris incognita",
       title    = "Correlation: ASV Count vs % Incognita by Country",
       subtitle = "Spearman correlation with permutation-derived significance") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text      = element_text(face = "bold"))

# -------------------------
# 5. Plot: null distribution per country
# -------------------------
p2 <- results$perm_results %>%
  ggplot(aes(x = perm_cor)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(data = results$stats,
             aes(xintercept = obs_cor),
             color = "red", linewidth = 1, linetype = "dashed") +
  geom_text(data = results$stats,
            aes(x     = obs_cor,
                y     = Inf,
                label = paste0("obs r=",   round(obs_cor, 2),
                               "\np=",     p_value,
                               "\np.adj=", p_adj)),
            hjust = -0.1, vjust = 1.5, color = "red", size = 2.5) +
  facet_wrap(~ Country) +
  labs(x        = "Permuted Spearman correlation",
       y        = "Count",
       title    = "Permutation Test: ASV Count vs % Incognita",
       subtitle = "Red line = observed correlation | n = 10,000 permutations per country") +
  xlim(-1,1) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"))

# -------------------------
# 6. Save
# -------------------------
ggsave("permutation_test.pdf",
       plot   = gridExtra::arrangeGrob(p1, p2, ncol = 1),
       width  = 12, height = 14)
message("Plot saved to: permutation_test.pdf")

write_csv(results$stats, "permutation_test_results.csv")
message("Results saved to: permutation_test_results.csv")