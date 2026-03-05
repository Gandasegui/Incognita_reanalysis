# Generatign ASVs from fastq files from Emerging Infectious Dieases (EDI) paper

To do so, we will use only forward reads. 
In the paper original paper by Venkatesan, they used fwd and rev reads, but the lenght of the PE reads is not enought to generate ASVS
This is not critinal to species assgiment, consideirn that the ITS sequences differ enought to classifiy the ASV by species.
However, the ASVs obtained should not be considered as comparable to the original paper.

```R
library(tidyverse)
library(readxl)
library(ggpubr)
library(ggsci)
library(dada2)


# STEP 1: Run DADA2 on EDI Illumina samples

setwd("~/R/ITS_trichuris")
fastq_dir <- "data/fastq_EDI"
meta_fp   <- "data/metadata/metadata_EDI.tsv"
# metadata_EDI.tsv is a three clum tsv with the fastq sra names, marker (ITS-1/2), and country
out_dir   <- "dada2_FWD_only_shared40"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
meta <- read_tsv(meta_fp, show_col_types = FALSE)

fnF_all <- sort(list.files(fastq_dir, pattern = "_1\\.fastq(\\.gz)?$", full.names = TRUE))
get_run <- function(x) str_remove(basename(x), "_1\\.fastq(\\.gz)?$")
run_ids <- vapply(fnF_all, get_run, character(1))
names(fnF_all) <- run_ids
fnF <- fnF_all[run_ids %in% meta$Run]
stopifnot(length(fnF) > 0)

filt_dir <- file.path(out_dir, "01_filtered")
dir.create(filt_dir, showWarnings = FALSE)
filtF <- file.path(filt_dir, paste0(names(fnF), "_F_filt.fastq.gz"))

out_filt <- filterAndTrim(fnF, filtF, truncLen = 240, maxEE = 3, truncQ = 2,
                          maxN = 0, minLen = 50, rm.phix = TRUE,
                          compress = TRUE, multithread = TRUE)
write_tsv(as_tibble(out_filt), file.path(out_dir, "01_filterAndTrim_summary.tsv"))

errF <- learnErrors(filtF, multithread = TRUE)
saveRDS(errF, file.path(out_dir, "02_errF.rds"))
derepF <- derepFastq(filtF)
names(derepF) <- names(fnF)
dadaF <- dada(derepF, err = errF, multithread = TRUE)

seqtab <- makeSequenceTable(dadaF)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
saveRDS(seqtab.nochim, file.path(out_dir, "03_seqtab_nochim_FWD.rds"))

# Helper to export marker-specific results
export_marker_simple <- function(marker_name, meta, seqtab, out_dir) {
  out_m <- file.path(out_dir, paste0("marker_", marker_name))
  dir.create(out_m, showWarnings = FALSE, recursive = TRUE)
  
  runs_m <- meta %>%
    filter(Marker == marker_name) %>%
    pull(Run) %>%
    intersect(rownames(seqtab))
  if (length(runs_m) == 0) return(NULL)
  
  st <- seqtab[runs_m, , drop = FALSE]
  st <- st[, colSums(st) > 0, drop = FALSE]
  asv_seqs <- colnames(st)
  asv_ids  <- paste0(marker_name, "_ASV_", sprintf("%05d", seq_along(asv_seqs)))
  
  write_tsv(as.data.frame(st) %>% setNames(asv_ids) %>% rownames_to_column("Run"),
            file.path(out_m, "ASV_table.tsv"))
  writeLines(c(rbind(paste0(">", asv_ids), asv_seqs)),
             file.path(out_m, "ASVs.fasta"))
  write_tsv(tibble(ASV_ID = asv_ids, sequence = asv_seqs),
            file.path(out_m, "ASV_id_to_sequence.tsv"))
}

meta2 <- meta %>% filter(Run %in% rownames(seqtab.nochim))
export_marker_simple("ITS-1", meta2, seqtab.nochim, out_dir)
export_marker_simple("ITS-2", meta2, seqtab.nochim, out_dir)
```
### From this scritp, the following elements should be obtained:
- The output directory "dada2_FWD_only_shared40", with marker specific subdirectories
- In each subdirectory, the marker-specific ASV table (ASV_table.tsv) and ASV sequences (ASVs.fasta)
