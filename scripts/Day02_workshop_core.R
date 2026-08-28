# ============================================================
# Jeju ARMS Day 2 — Taxonomy to community ecology
# ============================================================

library(vegan)
library(indicspecies)

# 1. Load data ------------------------------------------------

asv <- read.delim("data/Day2/ASV_table.tsv", check.names = FALSE)
metadata <- read.delim("data/Day2/metadata.tsv", stringsAsFactors = FALSE)
taxonomy <- read.delim("data/Day2/taxonomy_final_ASV.tsv",
                       stringsAsFactors = FALSE)

sample_names <- metadata$SampleID
metadata$Site <- factor(metadata$Site,
                        levels = c("Gangjeong", "Bomok"))

# Check BLAST results
nrow(taxonomy)
table(taxonomy$BLAST_status)
table(taxonomy$Status)


# 2. Ecological filtering ------------------------------------

keep <- taxonomy$BLAST_status == "Accepted" & 
  taxonomy$Status == "Include"
summary(keep)

target_taxonomy <- taxonomy[keep, ]

target_counts <- asv[
  match(target_taxonomy$ASV_ID, asv$ASV_ID),
  sample_names
]

nrow(target_taxonomy)    # 190 ASVs


# 3. Collapse ASVs into taxa ---------------------------------

taxon_counts <- rowsum(
  as.matrix(target_counts),
  target_taxonomy$Taxon_label
)

community <- t(taxon_counts)

dim(community)           # 8 samples x 88 taxa


# 4. Rarefaction ---------------------------------------------

site_colors <- c(Gangjeong = "#2C7FB8",Bomok = "#E67E22")

rarecurve(community, step = 100, label = FALSE,
          col = site_colors[metadata$Site],
          lwd = 2)
legend("bottomright", names(site_colors), 
       col = site_colors, lwd = 2, bty = "n")


# 5. Rarefy to equal sequencing depth ------------------------

sample_depth <- rowSums(community)
min_depth <- min(sample_depth)

sample_depth
min_depth

set.seed(20260821)

community_rare <- rrarefy(community, sample = min_depth)
rowSums(community_rare)


# 6. Alpha diversity -----------------------------------------

alpha <- data.frame(
  SampleID = sample_names,
  Site = metadata$Site,
  Observed = specnumber(community_rare),
  Shannon = diversity(community_rare, "shannon"),
  Simpson = diversity(community_rare, "simpson")
)

alpha

par(mfrow = c(1, 3))
for (metric in c("Observed", "Shannon", "Simpson")) {
    test <- wilcox.test(alpha[[metric]] ~ alpha$Site, exact = TRUE)
    boxplot(alpha[[metric]] ~ alpha$Site,
            col = adjustcolor(site_colors, alpha.f = 0.35),
            main = paste0(metric, "\np = ", signif(test$p.value, 3)), xlab = "",
            ylab = metric)
}

par(mfrow = c(1, 1))


# 7. Community composition -----------------------------------

taxon_group <- target_taxonomy$Group[
  match(rownames(taxon_counts), target_taxonomy$Taxon_label)]

composition_counts <- rowsum(taxon_counts, taxon_group, reorder = FALSE)

composition_relative <- sweep(composition_counts, 2,
                              colSums(composition_counts), "/")

# Plotting
composition_colors <- hcl.colors(nrow(composition_relative), "Set 3")
par(mar = c(9, 5, 4, 10), xpd = TRUE)
barplot(composition_relative, col = composition_colors, border = NA, las = 2,
        ylab = "Relative read proportion",
        main = "Community composition")
legend("topright",inset = c(-0.28, 0),
       legend = rev(rownames(composition_relative)),
       fill = rev(composition_colors),
       bty = "n",cex = 0.85)

# Reset graphics parameters
par(mar = c(5, 4, 4, 2), xpd = FALSE)

# 8. Bray-Curtis and nMDS ------------------------------------

community_log <- log1p(community)

bray <- vegdist(community_log, method = "bray")
print(bray)


set.seed(20260821)


nmds <- metaMDS(community_log, distance = "bray", k = 2, trymax = 100,
                autotransform = FALSE,  trace = FALSE)

nmds$stress

plot(nmds, type = "n")
points(nmds, display = "sites", pch = 19, cex = 1.5,
       col = site_colors[metadata$Site])

text(nmds, display = "sites", labels = sample_names, pos = 3, cex = 0.7)
legend("topright", names(site_colors), col = site_colors, pch = 19, bty = "n")    # Site legend
legend("topleft", legend = paste0("Stress = ", sprintf("%.4f", nmds$stress)),     # Stress value
       bty = "n", cex = 0.8)

# 9. ANOSIM --------------------------------------------------

set.seed(20260821)

anosim_result <- anosim(
  bray,
  metadata$Site,
  permutations = 999
)

anosim_result


# 10. Indicator taxa -----------------------------------------

set.seed(20260821)

indval <- multipatt(
  community,
  metadata$Site,
  func = "IndVal.g",
  control = permute::how(nperm = 999)
)

summary(indval, alpha = 0.05, indvalcomp = TRUE)

# 10-1. Select top 5 indicator taxa --------------------------

# Top 5 taxa from the IndVal summary
top_gj <- c(
  "Botryllus sp.",
  "Megabalanus volcano",
  "Sphacelariaceae",
  "Rhodomelaceae",
  "Papenfussiella kuromo"
)

top_bm <- c(
  "Scytosiphonaceae",
  "Tubuliporida",
  "Halopsis sp.",
  "Dendostrea sandvichensis",
  "Cirratulidae"
)

top_taxa <- c(top_gj, top_bm)

# 10-2. Heatmap ----------------------------------------------

heat_mat <- t(community_log[, top_taxa])

rownames(heat_mat) <- c(
  paste0("GJ | ", top_gj),
  paste0("BM | ", top_bm)
)

heat_colors <- colorRampPalette(c("white", "orange", "red"))(50)

heatmap(heat_mat, Rowv = NA, Colv = NA, scale = "none", col = heat_colors,
        margins = c(8, 18), main = "Top indicator taxa")

legend("topright", legend = c("Low", "Medium", "High"),
       fill = heat_colors[c(1, 25, 50)], title = "Abundance", bty = "n",
       cex = 0.8)

# 이쁘게 그리기
heat_mat <- t(community_log[, top_taxa])

rownames(heat_mat) <- c(
  paste0("GJ | ", top_gj),
  paste0("BM | ", top_bm)
)

# BM at top, GJ at bottom
heat_mat <- heat_mat[nrow(heat_mat):1, ]

heat_colors <- colorRampPalette(
  c("white", "orange", "red")
)(50)

par(mar = c(6, 18, 6, 3))

image(
  1:ncol(heat_mat),
  1:nrow(heat_mat),
  t(heat_mat),
  col = heat_colors,
  axes = FALSE,
  xlab = "",
  ylab = ""
)

# Sample names
axis(
  3,
  at = 1:ncol(heat_mat),
  labels = colnames(heat_mat),
  tick = FALSE,
  cex.axis = 1.1
)

# Taxon names
axis(
  2,
  at = 1:nrow(heat_mat),
  labels = rownames(heat_mat),
  las = 1,
  tick = FALSE,
  cex.axis = 1.1
)

box()

legend(
  "bottom",
  inset = -0.24,
  xpd = TRUE,
  horiz = TRUE,
  legend = c("Low", "Medium", "High"),
  fill = heat_colors[c(1, 25, 50)],
  title = "Abundance",
  bty = "n",
  cex = 1.1
)

