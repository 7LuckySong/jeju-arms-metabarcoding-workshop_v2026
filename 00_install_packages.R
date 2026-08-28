# ============================================================
# Workshop Package Installation Script (Day 1 & Day 2)
# ============================================================

# 1. Install CRAN packages (For Day 2 and data manipulation)
# - tidyverse: Data manipulation and visualization
# - vegan: Ecological community analysis and diversity metrics
# - indicspecies: Indicator species analysis
install.packages(c(
  "tidyverse",
  "vegan",
  "indicspecies"
))

# 2. Prepare Bioconductor installation
# Install BiocManager only if it is not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 3. Install Bioconductor packages (For Day 1)
# - dada2: FASTQ QC and ASV (Amplicon Sequence Variant) inference
BiocManager::install("dada2", ask = FALSE, update = FALSE)

# ============================================================
# 4. Installation Check
# ============================================================
# If the message below prints without any errors, 
# the installation was successfully completed.
library(tidyverse)
library(vegan)
library(indicspecies)
library(dada2)

packageVersion("dada2")

cat("\n[SUCCESS] All required packages for the workshop have been successfully installed and loaded!\n")
