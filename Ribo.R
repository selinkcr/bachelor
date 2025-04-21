# Load libraries
library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)

# Input BAM file
bam_file <- "WT_Rep1_rpf_UMI_less_rRNA_less_tRNA.bam"
bam <- readGAlignments(bam_file)

# A-site offset
offset <- 11

# Gene list to zoom into
genes_of_interest <- c(
  "BW25113_RS13540_WP_001235102.1_clpB", "BW25113_RS13545_WP_000040169.1_pgeF",
  "BW25113_RS13550_WP_000079100.1_rluD", "BW25113_RS13555_WP_000197686.1_bamD",
  "BW25113_RS13565_WP_000178456.1_raiA", "BW25113_RS00005_WP_001386572.1_thrL",
  "BW25113_RS00010_WP_001264707.1_thrA", "BW25113_RS00015_WP_000241662.1_thrB",
  "BW25113_RS00020_WP_000781074.1_thrC", "BW25113_RS00030_WP_000738719.1_yaaX",
  "BW25113_RS00035_WP_000906197.1_yaaA", "BW25113_RS10525_WP_000019197.1_plaP",
  "BW25113_RS12650_WP_000290230.1_cysP", "BW25113_RS13510_WP_000841103.1_kgtP"
)

# Frame colors for visualization
frame_colors <- c("darkgreen", "dodgerblue", "firebrick")

# Strand, 3' end, A-site
strand_bam <- strand(bam)
three_prime <- ifelse(strand_bam == "+", end(bam), start(bam))
a_site <- ifelse(strand_bam == "+", three_prime - offset, three_prime + offset)

# Create GRanges
a_site_gr <- GRanges(seqnames = seqnames(bam),
                     ranges = IRanges(start = a_site, end = a_site),
                     strand = strand_bam)

# Coverage
cov_reads <- coverage(bam)
cov_asite <- coverage(a_site_gr)

# ========== 1. ZOOMED COVERAGE PLOTS ==========
dir.create("gene_plots", showWarnings = FALSE)

for (gene in genes_of_interest) {
  if (gene %in% names(cov_reads)) {
    cov_r <- as.numeric(cov_reads[[gene]])
    cov_a <- as.numeric(cov_asite[[gene]])
    
    len <- min(length(cov_r), length(cov_a))
    len_to_plot <- min(150, len)
    
    frames <- (1:len_to_plot) %% 3
    bar_colors <- frame_colors[frames + 1]
    
    png(file = paste0("gene_plots/", gene, "_zoom_frames.png"), width = 1000, height = 600)
    par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
    
    barplot(cov_r[1:len_to_plot],
            col = bar_colors,
            border = NA,
            main = paste("Read Start Coverage (Frames):", gene),
            xlab = "Position", ylab = "Read Count")
    
    barplot(cov_a[1:len_to_plot],
            col = bar_colors,
            border = NA,
            main = paste("A-site Coverage (Frames):", gene),
            xlab = "Position", ylab = "A-site Count")
    
    dev.off()
  }
}

# ========== 2. TRIPLET PERIODICITY FRAME PLOT ==========
# (a) Read start frame
read_start <- ifelse(strand_bam == "+", start(bam), end(bam))
frame_reads <- table(read_start %% 3)

# (b) A-site frame
frame_asite <- table(a_site %% 3)

# Plot both
png("triplet_periodicity.png", width = 1000, height = 500)
par(mfrow = c(1, 2))
barplot(frame_reads,
        names.arg = c("Frame 0", "Frame 1", "Frame 2"),
        col = frame_colors,
        main = "Triplet Periodicity (Read Starts)",
        ylab = "Read Count")
barplot(frame_asite,
        names.arg = c("Frame 0", "Frame 1", "Frame 2"),
        col = frame_colors,
        main = "Triplet Periodicity (A-site Offset 11)",
        ylab = "Read Count")
dev.off()

# ========== 3. METAGENE PROFILE (First 25 nt of all genes) ==========

# Get matrix of first 25 positions for all genes (whole reads)
meta_matrix_reads <- sapply(names(cov_reads), function(g) {
  v <- as.numeric(cov_reads[[g]])
  if (length(v) >= 25) return(v[1:25]) else return(rep(NA, 25))
})
meta_matrix_reads <- meta_matrix_reads[, colSums(is.na(meta_matrix_reads)) == 0]
meta_profile_reads <- rowMeans(meta_matrix_reads)

# Get matrix for A-site
meta_matrix_asite <- sapply(names(cov_asite), function(g) {
  v <- as.numeric(cov_asite[[g]])
  if (length(v) >= 25) return(v[1:25]) else return(rep(NA, 25))
})
meta_matrix_asite <- meta_matrix_asite[, colSums(is.na(meta_matrix_asite)) == 0]
meta_profile_asite <- rowMeans(meta_matrix_asite)

# Plot both
png("metagene_profiles.png", width = 1000, height = 500)
par(mfrow = c(1, 2))
barplot(meta_profile_reads,
        col = frame_colors[((1:25 - 1) %% 3) + 1],
        main = "Metagene Profile (Read Start)",
        xlab = "Position from TSS", ylab = "Mean Coverage")
barplot(meta_profile_asite,
        col = frame_colors[((1:25 - 1) %% 3) + 1],
        main = "Metagene Profile (A-site)",
        xlab = "Position from TSS", ylab = "Mean A-site Coverage")
dev.off()
