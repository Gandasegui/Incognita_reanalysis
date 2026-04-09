# =========================================
# R Script: Minimum 5' Prefix Length for Discrimination
# =========================================

library(proxy)
library(ggplot2)

# ----------------------------
# 1. Load sequences from tab-delimited file
# ----------------------------
# Replace with your file path
file_path <- "all_asv_sequence.txt"  # your file
df <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)

# Extract sequences and labels
seqs <- df$Sequence
labels <- df$Species_ID  # use Species_ID as group label

# ----------------------------
# 2. Function to compute Hamming distance
# ----------------------------
hamming_dist <- function(s1, s2) {
  sum(strsplit(s1, "")[[1]] != strsplit(s2, "")[[1]]) / nchar(s1)
}

# ----------------------------
# 3. Function to compute distance matrix
# ----------------------------
compute_dist_matrix <- function(seqs_trimmed) {
  n <- length(seqs_trimmed)
  D <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n) {
    for(j in i:n) {
      d <- hamming_dist(seqs_trimmed[i], seqs_trimmed[j])
      D[i,j] <- d
      D[j,i] <- d
    }
  }
  return(D)
}

# ----------------------------
# 4. 1-NN leave-one-out classification
# ----------------------------
knn_loo_accuracy <- function(D, labels) {
  n <- nrow(D)
  correct <- 0
  for(i in 1:n) {
    dists <- D[i,]
    dists[i] <- Inf  # exclude self
    nn <- which.min(dists)
    if(labels[nn] == labels[i]) correct <- correct + 1
  }
  return(correct / n)
}

# ----------------------------
# 5. Loop over prefix lengths
# ----------------------------
#max_len <- min(nchar(seqs))  # maximum safe prefix
max_len <- 250
prefix_lengths <- seq(1, max_len, by = 1)
accuracy_results <- numeric(length(prefix_lengths))

for(k in seq_along(prefix_lengths)) {
  L <- prefix_lengths[k]
  # trim sequences to prefix
  seqs_trimmed <- substr(seqs, 1, L)
  
  # distance matrix
  D <- compute_dist_matrix(seqs_trimmed)
  
  # 1-NN accuracy
  acc <- knn_loo_accuracy(D, labels)
  accuracy_results[k] <- acc
  cat(sprintf("Prefix length %d bp -> Accuracy: %.3f\n", L, acc))
}

# ----------------------------
# 6. Find minimum prefix length with >=95% accuracy
# ----------------------------
threshold <- 0.95
min_idx <- which(accuracy_results >= threshold)[1]
if(!is.na(min_idx)) {
  cat(sprintf("\nMinimum prefix length for >= %.0f%% accuracy: %d bp\n", threshold*100, prefix_lengths[min_idx]))
} else {
  cat("\nNo prefix length reached the desired accuracy.\n")
}

# ----------------------------
# 7. Plot accuracy vs prefix length
# ----------------------------
df_plot <- data.frame(Prefix = prefix_lengths, Accuracy = accuracy_results)
ggplot(df_plot, aes(x = Prefix, y = Accuracy)) +
  geom_line(color = "blue", size = 1.2) +
  geom_point(color = "red", size = 2) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "darkgreen") +
  labs(title = "5' Prefix Length vs 1-NN Classification Accuracy",
       x = "Prefix length (bp)",
       y = "Accuracy") +
  theme_minimal()