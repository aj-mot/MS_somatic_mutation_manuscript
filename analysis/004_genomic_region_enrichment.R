#!/usr/bin/env Rscript

# Perform sSNV enrichment analysis vs gene expression and DNA replication timing

# sSNVs in each bin (gene expression and DNA replication timing deciles) are 
# fitted to the de novo mutational signatures by SCAN2

library("ggplot2")

################################################################################
# Enrichment of sSNVs in highly expressed genes using 54 tissues from GTEx

# Load results of SCAN2 permtool
load("../data/SCAN2_permtool_GTEx_mutsigs_results.rda")

snvs <- do.call(rbind, lapply(seq_along(es), function(i) {
  e <- es[[i]]
  e$bin <- sapply(strsplit(rownames(e), split="\\|\\|\\|"), function(z) z[1])
  e$sig <- sapply(strsplit(rownames(e), split="\\|\\|\\|"), function(z) z[2])
  e
  e <- e[e$bin %in% 1:10, ]
  e$bin <- as.integer(e$bin)
  cbind(e, Tissue=strsplit(basename(names(es)[i]), split="\\.")[[1]][1])
}))

snvs$brain.tissue <- sapply(snvs$Tissue, function(y) {
  substr(y, 1, 5)
})
snvs$brain.tissue[snvs$brain.tissue != "Brain"] <- "Nonbrain tissue"
snvs$brain.tissue[snvs$brain.tissue == "Brain"] <- "Brain tissue"

# Rename signatures
snvs$sig[snvs$sig == "MS-rank-3 sig 1"] <- "MS sig 2"
snvs$sig[snvs$sig == "MS-rank-3 sig 2"] <- "MS sig 3"
snvs$sig[snvs$sig == "MS-rank-3 sig 3"] <- "MS sig 1"

# Reorder signatures
snvs$sig <- factor(snvs$sig, levels = c("MS sig 1", "MS sig 2", "MS sig 3", "norm", "pct.resid", "resid.norm"), 
                   ordered = T)

# Enrichment in top decile
mean(snvs$enr[(snvs$bin == 10) & (snvs$sig == "MS sig 1") & (snvs$brain.tissue == "Brain tissue")]) # 1.434385  # 43% enrichment
mean(snvs$enr[(snvs$bin == 10) & (snvs$sig == "MS sig 2") & (snvs$brain.tissue == "Brain tissue")]) # 0.5586923 # 44% depletion
mean(snvs$enr[(snvs$bin == 10) & (snvs$sig == "MS sig 3") & (snvs$brain.tissue == "Brain tissue")]) # 0.6872308 # 31% depletion

mean(snvs$enr[(snvs$bin == 10) & (snvs$sig == "MS sig 1") & (snvs$brain.tissue != "Brain tissue")]) # 1.330439. # 33% enrichment
mean(snvs$enr[(snvs$bin == 10) & (snvs$sig == "MS sig 2") & (snvs$brain.tissue != "Brain tissue")]) # 0.5862439 # 41% depletion
mean(snvs$enr[(snvs$bin == 10) & (snvs$sig == "MS sig 3") & (snvs$brain.tissue != "Brain tissue")]) # 0.6725122 # 33% depletion

# Figure 4A - sSNV enrichment vs Gene expression
p <- ggplot(snvs[(snvs$bin %in% 1:10) & (snvs$sig %in% paste0("MS sig ", 1:3)), ], 
            aes(as.integer(bin), enr, group=Tissue, colour=brain.tissue)) +
  geom_line() +
  geom_line(data=snvs[(snvs$bin %in% 1:10) & (snvs$sig %in% paste0("MS sig ", 1:3)) & (snvs$brain.tissue == "Brain tissue"), ]) +
  geom_hline(yintercept = 1, lty=2, colour = "grey") +
  facet_grid(. ~ sig) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = c(0.92, 0.84)) +
  scale_x_continuous(breaks=seq(2, 10, by=2)) +
  scale_y_continuous(breaks=seq(0.6, 1.6, by=0.2)) +
  ylab("sSNVs (observed/expected)") + 
  xlab("GTEx gene expression decile")

print(p)
ggsave(plot = p, filename=file.path("..", "figures", "Fig_4A.pdf"), 
       width=9, height=3)
rm(emeta, es, snvs)

################################################################################
# Enrichment of sSNVs in genomic regions binned by DNA replication timing deciles (15 cell lines)

# Load results of SCAN2 permtool
load("../data/SCAN2_permtool_Repliseq_mutsigs_results.rda")

snvs <- do.call(rbind, lapply(seq_along(es), function(i) {
  e <- es[[i]]
  e$bin <- sapply(strsplit(rownames(e), split="\\|\\|\\|"), function(z) z[1])
  e$sig <- sapply(strsplit(rownames(e), split="\\|\\|\\|"), function(z) z[2])
  e <- e[e$bin %in% 1:10, ]
  e$bin <- as.numeric(e$bin)
  cbind(e, Cell_line=strsplit(basename(names(es)[i]), split="\\.")[[1]][1])
}))

# Reverse order of bins
snvs$bin <- 11 - snvs$bin

# Rename signatures
snvs$sig[snvs$sig == "MS-rank-3 sig 1"] <- "MS sig 2"
snvs$sig[snvs$sig == "MS-rank-3 sig 2"] <- "MS sig 3"
snvs$sig[snvs$sig == "MS-rank-3 sig 3"] <- "MS sig 1"

# Reorder signatures
snvs$sig <- factor(snvs$sig, levels = c("MS sig 1", "MS sig 2", "MS sig 3", "norm", "pct.resid", "resid.norm"), 
                   ordered = T)

# Enrichment in first decile - early replicating
mean(snvs$enr[(snvs$bin == 1) & (snvs$sig == "MS sig 1")]) # 1.1884    # 19% enrichment
mean(snvs$enr[(snvs$bin == 1) & (snvs$sig == "MS sig 2")]) # 0.6306    # 37% depletion
mean(snvs$enr[(snvs$bin == 1) & (snvs$sig == "MS sig 3")]) # 0.7764667 # 22% depletion  
# Enrichment in last decile - late replicating
mean(snvs$enr[(snvs$bin == 10) & (snvs$sig == "MS sig 1")]) # 0.9626   # 4% depletion
mean(snvs$enr[(snvs$bin == 10) & (snvs$sig == "MS sig 2")]) # 1.632933 # 63% enrichment
mean(snvs$enr[(snvs$bin == 10) & (snvs$sig == "MS sig 3")]) # 1.359867 # 36% enrichment

# Figure 4B - sSNV enrichment vs DNA replication timing
p <- ggplot(snvs[(snvs$bin %in% 1:10) & (snvs$sig %in% paste0("MS sig ", 1:3)), ], 
            aes(bin, enr, group=Cell_line)) +
  geom_line() +
  geom_smooth(aes(group=1), colour="red", se = F) +
  geom_hline(yintercept = 1, lty=2, colour = "grey") +
  facet_grid(. ~ sig) +
  theme(legend.position = "none") +
  theme_classic() +
  scale_x_continuous(breaks=seq(2, 10, by=2)) +
  scale_y_continuous(breaks=seq(0.6, 1.6, by=0.2)) +
  ylab("sSNVs (observed/expected)") + 
  xlab("Replication timing decile")

print(p)
ggsave(plot = p, filename=file.path("..", "figures", "Fig_4B.pdf"), 
       width=9, height=3)
rm(emeta, es, snvs)
