# Marine eDNA Metabarcoding Analysis Using R

## Jeju ARMS Case Study

This repository contains the lecture slides, R scripts, and practice datasets for a two-day workshop on marine eDNA metabarcoding analysis.

The workshop follows a complete analytical workflow, from raw sequencing reads to ecological interpretation, using a teaching dataset reconstructed from historical Jeju Autonomous Reef Monitoring Structures (ARMS) data.

The emphasis is on understanding the purpose of each processing step and interpreting the outputs. Detailed programming syntax, mathematical formulas, and software development are outside the scope of the workshop.

## Workshop overview

| Day | Topic | Main workflow |
| --- | --- | --- |
| Day 1 | NGS raw data to ASVs | FASTQ inspection → filtering → error learning → denoising → paired-read merging → chimera removal → ASV table and representative sequences |
| Day 2 | Taxonomy and community ecology | BLAST interpretation → taxonomic filtering → ecological filtering → community matrix → rarefaction → alpha diversity → community composition → beta diversity → indicator taxa |

By the end of the workshop, participants will be able to:

- explain the main steps in a marine eDNA metabarcoding workflow
- inspect FASTQ quality and infer ASVs with DADA2
- interpret precomputed BLAST results and apply taxonomic cutoffs
- distinguish target organisms from contaminants and non-target records
- construct a taxon-by-sample community matrix
- calculate and interpret alpha- and beta-diversity results
- connect R outputs to ecological questions

## Workshop format

- The workshop focuses on **data processing principles and ecological interpretation** rather than formulas or detailed programming theory.
- Participants run prepared R scripts and inspect the results step by step.
- No previous R programming experience is required.
- Bringing a laptop with R and RStudio installed is strongly recommended.
- Participants who cannot run R during the session can still follow the lecture and interpretation workflow.
- BLAST searches have been computed in advance. Installing BLAST+ or downloading a reference database is not required.

## Biological case study

The practice dataset is based on a previous Jeju ARMS study comparing sessile benthic communities from two contrasting habitats:

- **Gangjeong (GJ):** macroalgae-dominated habitat
- **Bomok (BM):** coral-dominated habitat

The Day 2 teaching dataset contains eight reconstructed samples: four Gangjeong samples and four Bomok samples. It is designed for education and does not represent a new independent biological survey.

## Software requirements

Install the following software before the workshop:

- [R](https://cran.r-project.org/)
- [RStudio Desktop](https://posit.co/download/rstudio-desktop/)

The required R packages are installed by `00_install_packages.R`. Key packages include:

- `dada2`
- `tidyverse`
- `vegan`
- `indicspecies`

Package installation is required only once. Windows users should keep the workshop's default `multithread = FALSE` settings.

## Download and start

### Recommended: download the workshop package

1. Open the **Releases** section of this repository.
2. Download the latest `Jeju_ARMS_Metabarcoding_Workshop.zip` file.
3. Extract the ZIP file completely. Do not run the project from inside the compressed folder.
4. Open `ARMS_Metabarcoding_Workshop.Rproj` in RStudio.
5. Open and run `00_install_packages.R` if the required packages are not already installed.
6. Open the relevant practice script for Day 1 or Day 2.

In RStudio, run the current line or selected lines with `Ctrl + Enter` on Windows or `Command + Enter` on macOS.

### Alternative: download from the repository

Select **Code → Download ZIP**, extract the downloaded archive, and open the `.Rproj` file in RStudio.

## Repository structure

```text
jeju-arms-metabarcoding-workshop/
├── README.md
├── ARMS_Metabarcoding_Workshop.Rproj
├── 00_install_packages.R
├── scripts/
│   ├── Day01_practice_core.R
│   └── Day02_practice_core.R
├── data/
│   ├── fastq/
│   └── day2/
│       ├── ASV_table.tsv
│       ├── ASV_representatives.fasta
│       ├── taxonomy_final_ASV.tsv
│       └── metadata.tsv
├── slides/
│   ├── ARMS_metabarcoding_Day1.pdf
│   └── ARMS_metabarcoding_Day2.pdf
└── outputs/
    └── day1_ASV_outputs/
```

Do not rename or move files after downloading the package. The practice scripts use this folder structure to locate their inputs and save their outputs.

## Important note about the two datasets

Day 1 and Day 2 use intentionally different teaching datasets.

- **Day 1** uses a reduced pseudo-FASTQ dataset so that participants can run the DADA2 workflow during the workshop.
- **Day 2** uses a prepared dataset reconstructed from historical Jeju ARMS ASV and BLAST results so that the workshop can focus on taxonomy and community ecology.

The ASV identifiers generated during Day 1 do not correspond directly to the prepared ASV identifiers used during Day 2. Do not join the Day 1 output to the Day 2 BLAST or taxonomy files by ASV ID.

## Day 1: NGS raw data to ASVs

Day 1 introduces the structure of paired-end FASTQ files and the DADA2 workflow:

1. inspect forward- and reverse-read quality;
2. filter and trim low-quality reads;
3. learn sequencing error rates;
4. dereplicate and denoise reads;
5. merge paired-end reads;
6. remove chimeric sequences; and
7. export an ASV abundance table and representative ASV sequences.

The main outputs are saved under `outputs/day1_ASV_outputs/`.

## Day 2: Taxonomy and community ecology

Day 2 begins with precomputed BLAST results and follows four ecological questions:

1. What organisms match our ASVs?
2. Which records belong in the ecological analysis?
3. How different are the Gangjeong and Bomok communities?
4. Which taxa characterize each site?

The workshop-specific taxonomic acceptance thresholds are:

| Accepted rank | Minimum sequence identity |
| --- | ---: |
| Species | 98% |
| Genus | 95% |
| Family | 90% |
| Order | 85% |

These values are teaching rules for this dataset, not universal metabarcoding thresholds. Appropriate cutoffs should be evaluated for each marker, reference database, and study objective.

Ecological analyses include rarefaction, observed richness, Shannon and Simpson indices, community composition, log-transformed Bray–Curtis dissimilarity, nMDS, ANOSIM, and indicator-taxon analysis.

## Troubleshooting

- **R cannot find a file:** confirm that the ZIP file was fully extracted and that the `.Rproj` file was opened before running the script.
- **A package is missing:** run `00_install_packages.R`, restart RStudio, and try again.
- **The working directory is incorrect:** reopen `ARMS_Metabarcoding_Workshop.Rproj`; do not use `setwd()` to point to individual folders.
- **Your result differs slightly:** small differences can occur between package versions. Focus on the overall pattern and interpretation.
- **A step takes too long:** wait for the current command to finish before running the next block.

For unresolved problems, open an issue in this repository and include the error message, the script name, and the step where the error occurred.

## Reference

The biological context and source study are described in:

Lee, K.-T., Kim, T., Park, G.-H., Oh, C., Park, H.-S., Kang, D.-H., Kang, H.-S., & Yang, H.-S. (2024). Assessment of Sessile Benthic Communities in Jeju Island, Republic of Korea, Using Autonomous Reef Monitoring Structures (ARMS). *Diversity, 16*, 83. <https://doi.org/10.3390/d16020083>

## Data-use note

The files in this repository are provided for workshop and educational use. The practice dataset has been reduced and reconstructed to support classroom execution and interpretation. Results obtained from this teaching dataset should not be treated as independent estimates from the complete original study.

