# ============================================================
# Jeju ARMS Day 1 — DADA2 ASV Inference (Core Practice)
# ============================================================

library(dada2)  # fastq QC and ASV inference

# 1. Setup paths and load data -------------------------------

if (!exists("WORKSHOP_ROOT")) WORKSHOP_ROOT <- getwd()
raw_path <- file.path(WORKSHOP_ROOT, "data", "fastq")
out_path <- file.path(WORKSHOP_ROOT, "outputs", "day1_ASV_outputs")
filt_path <- file.path(out_path, "filtered_fastq")

dir.create(filt_path, recursive = TRUE, showWarnings = FALSE)

r1_pattern <- "(_R1(_001)?|_1)\\.(fastq)(\\.gz)?$"
r2_pattern <- "(_R2(_001)?|_2)\\.(fastq)(\\.gz)?$"
end_pattern <- "(_R[12](_001)?|_[12])\\.(fastq)(\\.gz)?$"

fwd <- sort(list.files(raw_path, r1_pattern, full.names = TRUE, ignore.case = TRUE))
rev <- sort(list.files(raw_path, r2_pattern, full.names = TRUE, ignore.case = TRUE))

if (!length(fwd) || length(fwd) != length(rev)) {
  stop("Error: Please check the number of R1/R2 FASTQ files.")
}

sample_id <- sub(end_pattern, "", basename(fwd), ignore.case = TRUE)


# 2. Inspect Read Quality (Raw) ------------------------------

plotQualityProfile(fwd[1:2])
plotQualityProfile(rev[1:2])


# 3. Filter and Trim -----------------------------------------

primer_len <- c(26, 26)  # COI FWD and REV primer lengths
trunc_len <- c(210, 180) 

filt_f <- file.path(filt_path, paste0(sample_id, "_F.fastq.gz"))
filt_r <- file.path(filt_path, paste0(sample_id, "_R.fastq.gz"))

filtered <- filterAndTrim(
  fwd, filt_f, rev, filt_r,
  truncLen = trunc_len, trimLeft = primer_len,
  maxN = 0, maxEE = c(2, 4), truncQ = 2,
  rm.phix = TRUE, matchIDs = TRUE, compress = TRUE)


# 4. Compare Before & After Trimming -------------------------

# Plot 4 panels in one window for the first sample
# (Forward Raw, Forward Filtered, Reverse Raw, Reverse Filtered)
plotQualityProfile(c(fwd[1], filt_f[1], rev[1], filt_r[1]))

# 5. Learn error rates ---------------------------------------

err_f <- learnErrors(filt_f)
err_r <- learnErrors(filt_r)

plotErrors(err_f, nominalQ = TRUE)
plotErrors(err_r, nominalQ = TRUE)


# 6. Dereplication and ASV inference (Denoise) -------------------------
# 1) Dereplication
derepFs <- derepFastq(filt_f)
derepRs <- derepFastq(filt_r)

# 2) ASV inference
asv_f <- dada(derepFs, err = err_f)
asv_r <- dada(derepRs, err = err_r)


# 7. Merge paired reads --------------------------------------

min_overlap <- 20

merged <- mergePairs(
  asv_f, derepFs, asv_r, derepRs,
  minOverlap = min_overlap, maxMismatch = 0)

names(merged) <- sample_id
sequence_table <- makeSequenceTable(merged)


# 8. Remove chimeras -----------------------------------------

asv_matrix <- removeBimeraDenovo(
  sequence_table, method = "consensus")


# 9. Export ASV table and representative sequences -----------

abundance_order <- order(colSums(asv_matrix), decreasing = TRUE)
asv_matrix <- asv_matrix[, abundance_order, drop = FALSE]

sequences <- colnames(asv_matrix)
asv_names <- paste0("ASV", seq_along(sequences))
colnames(asv_matrix) <- asv_names

asv_table <- data.frame(
  ASV_ID = asv_names, t(asv_matrix),
  check.names = FALSE, row.names = NULL)

write.table(
  asv_table, file.path(out_path, "05_ASV_table_inferred.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(
  as.vector(rbind(paste0(">", asv_names), sequences)),
  file.path(out_path, "05_ASV_representative_sequences.fasta"))


# 10. Read tracking ------------------------------------------

count_reads <- function(x) as.numeric(sum(getUniques(x)))

tracking <- data.frame(
  SampleID = sample_id,
  Input = filtered[, "reads.in"],
  Filtered = filtered[, "reads.out"],
  Denoised_forward = vapply(asv_f, count_reads, numeric(1)),
  Denoised_reverse = vapply(asv_r, count_reads, numeric(1)),
  Merged = vapply(merged, count_reads, numeric(1)),
  Non_chimeric = as.numeric(rowSums(asv_matrix)[sample_id]),
  row.names = NULL)

write.table(
  tracking, file.path(out_path, "06_read_tracking.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n[SUCCESS] Completed:", nrow(asv_matrix), "samples /", ncol(asv_matrix), "ASVs\n")
cat("Note: For Day 2, we will use the pre-prepared ASV table.\n")