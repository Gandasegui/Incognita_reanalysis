# =========================================
# R Script: Comparison of ONT / Illumina sequencing depth from matched samples
# =========================================

library(tidyverse)
library(scales)

df <- read.delim("depth_EDI_Nat.txt", header = TRUE, sep = "\t")

max_depth <- max(c(df$illumina_depth, df$ONT_depth), na.rm = TRUE)

ggplot(df, aes(x = illumina_depth, y = ONT_depth, color = ONT_incognita_pct, shape = Country)) +
  geom_point(size = 3, alpha = 0.9) +
  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  
  scale_x_continuous(
    limits = c(0, max_depth),
    labels = label_comma()
  ) +
  scale_y_continuous(
    limits = c(0, max_depth),
    labels = label_comma()
  ) +
  
  # 🎨 nice continuous colour scale
  scale_color_viridis_c(
    name = "% Trichuris incognita",
    option = "viridis"
  ) +
  
  labs(
    title = "Sequencing Depth Comparison",
    x = "Illumina Depth",
    y = "ONT Depth"
  ) +
  theme_minimal()