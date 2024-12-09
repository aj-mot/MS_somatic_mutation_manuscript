#!/usr/bin/env Rscript

# Perform mutational signature analysis
# Analysis of MS/CTL PTA neurons combined with 33 PTA neurotypical neurons from 
# Luquette et al.

library("MutationalPatterns")
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library("gridExtra")
library("NMF")
library("ggplot2")
library("plyr")
library("lme4")
library("lmerTest")
library("reshape2")
library("RColorBrewer")
library("pracma")

################################################################################
# Load sSNV calls

# List of SCAN2 VCF files of sSNVs
MS.PTA.vcf.files <- list.files(file.path("..", "data", "vcf"), 
                               pattern = ".scan2.vcf", full.names = TRUE)
MS.PTA.sample.names <- sapply(basename(MS.PTA.vcf.files), function(x) strsplit(x, ".scan2.vcf")[[1]][1])

Luquette.PTA.vcf.files <- list.files(file.path("..", "external_data", "vcf"), 
                                     pattern = ".scan2.vcf", full.names = TRUE)
Luquette.PTA.sample.names <- sapply(basename(Luquette.PTA.vcf.files), function(x) strsplit(x, ".scan2.vcf")[[1]][1])

PTA.vcf.files <- c(MS.PTA.vcf.files, Luquette.PTA.vcf.files)
PTA.sample.names <- c(MS.PTA.sample.names, Luquette.PTA.sample.names)

# Load VCF files into a GRangesList
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
grl <- read_vcfs_as_granges(PTA.vcf.files, PTA.sample.names, ref_genome)

# 96 mutational profile
# Matrix 96 x #samples
mut_mat_orig <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
mut_mat <- mut_mat_orig + 0.0001

# Define metadata
MS.sample.info <- read.csv(file.path("..", "data", "MS_sSNV_burden_data.csv"))

# Read external control data
luquette.data <- read.csv(file.path("..", "external_data", "luquette_sSNV_burden_data.csv"))
MS.sample.info <- rbind.fill(data.frame(MS.sample.info, Lab="Rubio"), 
                             data.frame(luquette.data, Lab="Walsh"))

MS.sample.info$is.MSnormal <- MS.sample.info$Tissue == "MS normal"
MS.sample.info$is.MSlesion <- MS.sample.info$Tissue == "MS lesion"
MS.sample.info$Tissue[MS.sample.info$Lab == "Walsh"] <- "External control"
MS.sample.info$Tissue[MS.sample.info$Lab == "Rubio" & MS.sample.info$Disease == "CTL"] <- "Control"
MS.sample.info$is.Control <- MS.sample.info$Tissue == "Control"

# Order individuals by age
ind.order <- unique(MS.sample.info$Individual)
ind.order <- ind.order[order(MS.sample.info$Age[match(ind.order, MS.sample.info$Individual)])]
MS.sample.info$Individual <- factor(MS.sample.info$Individual, 
                                    levels = ind.order, ordered = T)

# Order samples by plot category, age, name
sample.order <- MS.sample.info$Sample[order(MS.sample.info$Tissue, MS.sample.info$Age, MS.sample.info$Sample)]
MS.sample.info$ordered.Sample <- factor(MS.sample.info$Sample, 
                                        levels=sample.order, ordered=T)

# Order Tissue
tissue.order <- c("External control", "Control", "MS normal", "MS lesion")
MS.sample.info$Tissue <- factor(MS.sample.info$Tissue, 
                                levels = tissue.order, ordered = T)

ms.tissue.colours <- c(`External control`="blueviolet", `Control`="steelblue1", 
                       `MS normal`="olivedrab3", `MS lesion`="coral1")

# Order metadata rows to match mut_mat
MS.sample.info <- MS.sample.info[match(PTA.sample.names, MS.sample.info$Sample), ]

################################################################################
# Load published/known mutational signatures: 
# -COSMIC
# -Luquette et al, PTAerr (Universal PTA artifacts)
# -Petljak et al, SBS.sc_E, SBS.sc_F
# -Lodato Sig A, B, C 
# -Miller et al 4 sigs

# COSMIC signatures (60 sigs + 18 possible artifact sigs = total 78)
known_signatures <- get_known_signatures(muttype = "snv", incl_poss_artifacts = T)
# Rows will correspond to mut count matrix

# Add MDA artefact signatures, SBS.sc_E, SBS.sc_F
petljak.sigs <- read.csv(file.path("..", "external_data", "Petljak2019_signatures.csv"))
known_signatures <- cbind(known_signatures, as.matrix(petljak.sigs[, -(1:3)]))

# Add Universal PTA artifact signature (PTAerr) Luquette et al
lysis.sig <- read.csv(file.path("..", "external_data", "PTA_sSNV_artifact_signature.csv"))
known_signatures <- cbind(known_signatures, PTAerr=lysis.sig)

# Lodato sigs
lodato.sigs <- read.csv(file.path("..", "external_data", "Lodato2018_SignatureData_Aging.csv"))
known_signatures <- cbind(known_signatures, lodato.sigs[, -(1:2)])

# Miller sigs
miller.sigs.W <- read.csv(file.path("..", "external_data", "Miller_Sig_W1_W4.csv"))
miller.sigs.N <- read.table(file.path("..", "external_data", "Miller_Sig_N1_N4.tsv"), 
                            sep="\t", header=T)
known_signatures <- cbind(known_signatures, miller.sigs.N[, -1], miller.sigs.W[, -1])

known_signatures.cols <- colnames(known_signatures)  
known_signatures <- as.matrix(do.call(cbind, lapply(1:ncol(known_signatures), function (i) {
  known_signatures[, i] / sum(known_signatures[, i])
})))
rownames(known_signatures) <- paste(lodato.sigs$X, lodato.sigs$X.1, sep=":")
colnames(known_signatures) <- known_signatures.cols

################################################################################
# COSMIC mutational signature refitting

# Mutational signature analysis diagnostics
# Determine an appropriate number of COSMIC signatures to use for signature 
# refitting. Perform forward-stepwise regression, iteratively adding COSMIC 
# signatures and computing the reconstructed profile SSE.

# Fit 78 COSMIC signatures + PTA artefact signature (Luquette et al)
sigs2fit <- known_signatures[, c(1:78, 82)]
# Set all sigs to sum to 1
for (i in 1:ncol(sigs2fit)) {
  sigs2fit[, i] <- sigs2fit[, i] / sum(sigs2fit[, i])
}

# Start with 0 fitted signatures - constant profile
sig.used <- rep(F, ncol(sigs2fit))
resid.0 <- sum(mut_mat^2)
prev.resid <- resid.0

# Find the first best-fitting signature
next.resids <- sapply(which(!sig.used), function(i) {
  this.ssq <- do.call(sum, lapply(1:ncol(mut_mat), function(j) {
    lsqnonneg(sigs2fit[, i, drop=F], mut_mat[, j])$resid.norm
  }))
})
sig.used[!sig.used][which.min(next.resids)] <- T
print(colnames(sigs2fit)[which(sig.used)]) # "SBS5" is best fitting signature
next.resid <- min(next.resids)

# Iteratively add signatures to the fit
fit.resid <- rep(NA, ncol(sigs2fit))
fit.resid[1] <- next.resid
i <- 2
while ((prev.resid - next.resid) / prev.resid > 0.01) {
  next.resids <- sapply(which(!sig.used), function(i) {
    this.sig.used <- sig.used
    this.sig.used[i] <- T
    this.ssq <- do.call(sum, lapply(1:ncol(mut_mat), function(j) {
      lsqnonneg(sigs2fit[, this.sig.used, drop=F], mut_mat[, j])$resid.norm
    }))
  })
  sig.used[!sig.used][which.min(next.resids)] <- T
  print(colnames(sigs2fit)[which(sig.used)])
  prev.resid <- next.resid
  next.resid <- min(next.resids)
  fit.resid[i] <- next.resid
  i <- i + 1
}
fit.resid.improvement <- - fit.resid[2:length(fit.resid)] + fit.resid[1:(length(fit.resid)-1)]

# Extended Data Figure 6A - Reconstructed profile SSE vs Number of COSMIC signatures

pdf(file.path("..", "figures", "Ext_Data_Fig_6A.pdf"), 
    width = 6, height = 6)
plot(2:(sum(!is.na(fit.resid.improvement))+1), fit.resid.improvement[!is.na(fit.resid.improvement)], 
     type="l",
     xlab = "Number of COSMIC signatures", ylab = "Improvement in reconstructed profile SSE")
abline(h=2000, col="red")
dev.off()

# Perform signature refitting with 10 COSMIC signatures
# Signatures to fit, based on forward-stepwise regression above:
use.cosmic.sigs <- c("SBS1", "SBS5", "SBS12", "SBS16", "SBS19", "SBS30", "SBS42", "SBS54", "SBS88", "SBS89")

pta_cosmic_sigs <- known_signatures[, use.cosmic.sigs]
pta_cosmic_fit <- fit_to_signatures(mut_mat, pta_cosmic_sigs)

# Scale contributions so that sum of contributions for each neuron equals the 
# total estimated sSNV burden
pta_cosmic_fit_contribution <- pta_cosmic_fit$contribution
for (i in 1:ncol(pta_cosmic_fit_contribution)) {
  pta_cosmic_fit_contribution[, i] <- pta_cosmic_fit_contribution[, i] * 
    MS.sample.info$burden[match(colnames(pta_cosmic_fit_contribution)[i], MS.sample.info$Sample)] / 
    sum(pta_cosmic_fit$reconstructed[, i])
}
# Order signatures by contribution
sig.order <- order(rowSums(pta_cosmic_fit_contribution), decreasing = F)
pta_cosmic_fit_contribution <- pta_cosmic_fit_contribution[sig.order, ]
pta_cosmic_sigs <- pta_cosmic_sigs[, sig.order]

# Order samples by Tissue / Age / Individual / Sample
sample.order2 <- order(
  MS.sample.info$Tissue[match(colnames(pta_cosmic_fit_contribution), MS.sample.info$Sample)], 
  MS.sample.info$Age[match(colnames(pta_cosmic_fit_contribution), MS.sample.info$Sample)], 
  MS.sample.info$Individual[match(colnames(pta_cosmic_fit_contribution), MS.sample.info$Sample)], 
  MS.sample.info$Sample[match(colnames(pta_cosmic_fit_contribution), MS.sample.info$Sample)])
pta_cosmic_fit_contribution <- pta_cosmic_fit_contribution[, sample.order2]

sig.colours <- brewer.pal(10, "Set3")
names(sig.colours) <- colnames(pta_cosmic_sigs)

# Make data frame for plotting
pta_cosmic_fit_contribution.df <- as.data.frame(t(pta_cosmic_fit_contribution))
pta_cosmic_fit_contribution.df$Sample <- rownames(pta_cosmic_fit_contribution.df)

pta_cosmic_fit_contribution.df <- merge(pta_cosmic_fit_contribution.df, MS.sample.info)
pta_cosmic_fit_contribution.df.long <- melt(pta_cosmic_fit_contribution.df, id.vars = colnames(MS.sample.info))

# Figure 3A - COSMIC signature contribution to each neuron

p <- ggplot(pta_cosmic_fit_contribution.df.long, 
            aes(ordered.Sample, value, fill=variable)) +
  geom_col() +
  ylab("Number of mutations") + xlab("") +
  scale_fill_manual(labels=names(sig.colours), values=sig.colours) +
  facet_grid(~ Tissue, scales = "free_x", space = "free_x") +
  theme_classic() + 
  scale_y_continuous(limits = c(0, 4000), expand = c(0, 0)) +
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position="right") +
  theme(legend.title=element_blank())
print(p)
ggsave(p, filename = file.path("..", "figures", "Fig_3A.pdf"), 
       width = unit(15, "cm"), height = unit(6, "cm"))


# Figure 3B - COSMIC signature contribution vs Age

sbs.plots <- lapply(unique(pta_cosmic_fit_contribution.df.long$variable)[10:1], function(x) {
  
  p <- ggplot(subset(pta_cosmic_fit_contribution.df.long, variable == x), 
              aes(Age, value, colour=Tissue)) +
    geom_point() +
    geom_smooth(method="lm", se = F) +
    theme_classic() +
    scale_colour_manual(labels=names(ms.tissue.colours), values=ms.tissue.colours) +
    theme(legend.position=c(0.22, 0.9)) +
    ylab(x) +
    theme(legend.title=element_blank())
  if (x != "SBS5") p <- p + guides(colour = F)
  p
})
print(sbs.plots[[1]])

ggsave(
  filename = file.path("..", "figures", "Fig_3B.pdf"), 
  plot = marrangeGrob(sbs.plots[c(1, 2, 3, 6, 10)], nrow=1, ncol=5, top = ""), 
  width = 24, height = 5
)

################################################################################
# De novo mutational signature analysis

# NMF for all PTA data jointly (n=106)
# Extract de novo signatures with MS/CTL PTA neurons (n=73) and external 
# controls (Luquette et al, n=33).

# Mutational signature analysis diagnostics
# Extended Data Figure 6B - NMF rank survey

estimate <- nmf(mut_mat, rank = 1:10, method = "brunet", 
                nrun = 100, seed = 123456, .opt = "v-p")

pdf(file.path("..", "figures", "Ext_Data_Fig_6B.pdf"),
    height=6, width=6)
plot(estimate)
dev.off()


# Perform signature extraction
nmf_res.rank3 <- extract_signatures(mut_mat, rank = 3, nrun = 100, single_core = TRUE)
colnames(nmf_res.rank3$signatures) <- paste("MS-rank-3 sig", 1:3)

# Reorder/rename signatures
rank3.reorder.sigs <- c(3, 1, 2)
# Sig 1 is age-associated
# Sig 2 shows a disease-effect
# Sig 3 shows a lesion-effect
sigs.rank3 <- nmf_res.rank3$signatures[, rank3.reorder.sigs]
colnames(sigs.rank3) <- paste("MS sig", 1:3)

# Figure 3C - De novo mutational signatures
p1 <- plot_96_profile(sigs.rank3, condensed = TRUE, ymax=0.125)
print(p1)
ggsave(plot = p1, width=6, height=4, 
       filename=file.path("..", "figures", "Fig_3C.pdf"))

write.csv(sigs.rank3, file.path("..", "data", "MS_de_novo_mut_signatures.csv"), 
          row.names = T, quote = F)


# Compute signature contributions to sSNV burden of each neuron

# Scale signature contributions so that sum of all signature contributions 
# is equal to the total sSNV burden in each neuron. This takes into account
# that the #called sSNVs used for signature extraction is not the same as 
# the total estimated #sSNVs.

MS.sample.info$`MS sig 1` <- sapply(1:nrow(MS.sample.info), function(i) {
  sn <- MS.sample.info$Sample[i]
  bd <- MS.sample.info$burden[i]
  s3.contrib <- nmf_res.rank3$contribution[3, match(sn, colnames(nmf_res.rank3$contribution))] * 
    sum(nmf_res.rank3$signatures[, 3])
  total.contrib <- do.call(sum, lapply(1:3, function(j) nmf_res.rank3$contribution[j, match(sn, colnames(nmf_res.rank3$contribution))] * 
                                         sum(nmf_res.rank3$signatures[, j])))
  s3.contrib / total.contrib * bd
})

MS.sample.info$`MS sig 2` <- sapply(1:nrow(MS.sample.info), function(i) {
  sn <- MS.sample.info$Sample[i]
  bd <- MS.sample.info$burden[i]
  s1.contrib <- nmf_res.rank3$contribution[1, match(sn, colnames(nmf_res.rank3$contribution))] * 
    sum(nmf_res.rank3$signatures[, 1])
  total.contrib <- do.call(sum, lapply(1:3, function(j) nmf_res.rank3$contribution[j, match(sn, colnames(nmf_res.rank3$contribution))] * 
                                         sum(nmf_res.rank3$signatures[, j])))
  s1.contrib / total.contrib * bd
})

MS.sample.info$`MS sig 3` <- sapply(1:nrow(MS.sample.info), function(i) {
  sn <- MS.sample.info$Sample[i]
  bd <- MS.sample.info$burden[i]
  s2.contrib <- nmf_res.rank3$contribution[2, match(sn, colnames(nmf_res.rank3$contribution))] * 
    sum(nmf_res.rank3$signatures[, 2])
  total.contrib <- do.call(sum, lapply(1:3, function(j) nmf_res.rank3$contribution[j, match(sn, colnames(nmf_res.rank3$contribution))] * 
                                         sum(nmf_res.rank3$signatures[, j])))
  s2.contrib / total.contrib * bd
})


# Make data frame for plotting
contribs.rank3.long <- melt(MS.sample.info, measure.vars = paste("MS sig", 1:3))

# Order the signatures
contribs.rank3.long$variable <- factor(contribs.rank3.long$variable, 
                                       levels=paste("MS sig", 1:3), ordered=T)


# Figure 3D - De novo mutational signature contribution vs Age

# Plot contributions for each signature
rank3.contrib.plots <- lapply(unique(contribs.rank3.long$variable), function(x) {
  
  p <- ggplot(subset(contribs.rank3.long, variable == x), aes(Age, value, colour=Tissue)) +
    geom_point() +
    geom_smooth(method="lm", se = F) +
    theme_classic() +
    scale_colour_manual(labels=names(ms.tissue.colours), values=ms.tissue.colours) +
    theme(legend.position=c(0.22, 0.9)) +
    ylab(x) +
    theme(legend.title=element_blank())
  if (x != "MS sig 1") p <- p + guides(colour = "none")
  p
})
print(rank3.contrib.plots[[1]])

ggsave(
  filename = file.path("..", "figures", "Fig_3D.pdf"), 
  plot = marrangeGrob(rank3.contrib.plots, nrow=1, ncol=3, top = ""), 
  width = 12, height = 4
)

# Extended Data Figure 7 - Cosine similarity between de novo signatures and 
# published signatures

# Heatmap: 
# columns - de novo sigs 1-3
# rows    - published single neuron sigs (A, N4, W3, C, N2, W2, PTAerr)
#           + COSMIC sigs

sigs.rank3_pub_cosmic_sim <- cos_sim_matrix(known_signatures[, c(83, 89, 92, 85, 87, 91, 82, 1:78)], sigs.rank3)
rownames(sigs.rank3_pub_cosmic_sim)[c(1, 4)] <- paste("Signature", rownames(sigs.rank3_pub_cosmic_sim)[c(1, 4)]) 
rownames(sigs.rank3_pub_cosmic_sim) <- gsub("\\.", " ", rownames(sigs.rank3_pub_cosmic_sim))
p <- plot_cosine_heatmap(sigs.rank3_pub_cosmic_sim, 
                          col_order = ,
                          cluster_cols = F, cluster_rows = F,
                          plot_values = T)

print(p)
ggsave(p, filename = file.path("..", "figures", "Ext_Data_Fig_7.pdf"), 
       width = unit(8, "cm"), height = unit(20, "cm"))


################################################################################
# Refitting of de novo mutational signatures
# Fit the COSMIC signatures (including artefact signatures) + PTAerr (PTA artefact sig from Luquette et al)

denovo.rank3.cosmic.fit <- fit_to_signatures_strict(sigs.rank3, known_signatures[, c(1:78, 82)])
denovo.rank3.cosmic.fit.contrib <- data.frame(denovo.rank3.cosmic.fit$fit_res$contribution[rowMax(denovo.rank3.cosmic.fit$fit_res$contribution) > 0, ])
for (i in 1:ncol(denovo.rank3.cosmic.fit.contrib)) denovo.rank3.cosmic.fit.contrib[, i] <- denovo.rank3.cosmic.fit.contrib[, i] / sum(denovo.rank3.cosmic.fit.contrib[, i])
denovo.rank3.cosmic.fit.contrib$cosmic <- rownames(denovo.rank3.cosmic.fit.contrib)
denovo.rank3.cosmic.fit.contrib

# Extended Data Figure 8A - COSMIC signature contribution to de novo signatures

denovo.rank3.cosmic.fit.contrib.long <- melt(denovo.rank3.cosmic.fit.contrib, id.vars="cosmic")
denovo.rank3.cosmic.fit.contrib.long$cosmic <- factor(denovo.rank3.cosmic.fit.contrib.long$cosmic, 
                                                      levels=unique(denovo.rank3.cosmic.fit.contrib.long$cosmic), ordered=T)
denovo.rank3.cosmic.fit.contrib.long$variable <- sapply(denovo.rank3.cosmic.fit.contrib.long$variable, function(x)
  gsub("\\.", " ", x))

cosmic.colours <- c(brewer.pal(8, "Set2"), 
                    brewer.pal(12, "Set3")[c(2, 8, 9, 10, 11)], 
                    brewer.pal(12, "Set3")[4])

p <- ggplot(subset(denovo.rank3.cosmic.fit.contrib.long, value > 0), 
            aes(variable, value, fill=cosmic)) +
  geom_col() +
  geom_text(aes(label = cosmic),
            position = position_stack(vjust = .5)) +
  ylab("COSMIC contribution") + xlab("De novo signature") +
  scale_fill_manual(values=cosmic.colours) +
  theme_classic() + 
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(), axis.ticks.x=element_blank()) +
  theme(legend.position="right") +
  theme(legend.title=element_blank()) +
  guides(fill="none")
print(p)
ggsave(p, filename = file.path("..", "figures", "Ext_Data_Fig_8A.pdf"), 
       width = unit(4, "cm"), height = unit(4, "cm"))


# Plot the de novo mutational signature profiles alongside potential constituent 
# COSMIC sigs

# Extended Data Figure 8B - COSMIC signature contribution to MS sig 1
de_novo_1_plus_cosmic <- cbind(sigs.rank3[, 1, drop=F], known_signatures[, match(paste0("SBS", c(5, 12, 16, 54, 88, 89)), colnames(known_signatures))])
p1 <- plot_96_profile(de_novo_1_plus_cosmic, condensed = T)
p1 <- p1 + facet_grid(sample ~ substitution, scales = "free_y") +
  coord_cartesian(ylim = NULL) + 
  scale_y_continuous(breaks = waiver())
print(p1)
ggsave(p1, filename = file.path("..", "figures", "Ext_Data_Fig_8B.pdf"), 
       width = unit(6, "cm"), height = unit(8, "cm"))

# Extended Data Figure 8C - COSMIC signature contribution to MS sig 2
de_novo_2_plus_cosmic <- cbind(sigs.rank3[, 2, drop=F], known_signatures[, match(paste0("SBS", c(1, 30, 44, 89)), colnames(known_signatures))])
p2 <- plot_96_profile(de_novo_2_plus_cosmic, condensed = T)
p2 <- p2 + facet_grid(sample ~ substitution, scales = "free_y") +
  coord_cartesian(ylim = NULL) + 
  scale_y_continuous(breaks = waiver())
print(p2)
ggsave(p2, filename = file.path("..", "figures", "Ext_Data_Fig_8C.pdf"), 
       width = unit(6, "cm"), height = unit(5.8, "cm"))

# Extended Data Figure 8D - COSMIC signature contribution to MS sig 3
de_novo_3_plus_cosmic <- cbind(sigs.rank3[, 3, drop=F], known_signatures[, match(c(paste0("SBS", c(1, 8, 11, 12, 19, 40)), "PTAerr"), colnames(known_signatures))])
p3 <- plot_96_profile(de_novo_3_plus_cosmic, condensed = T)
p3 <- p3 + facet_grid(sample ~ substitution, scales = "free_y") +
  coord_cartesian(ylim = NULL) + 
  scale_y_continuous(breaks = waiver())
print(p3)
ggsave(p3, filename = file.path("..", "figures", "Ext_Data_Fig_8D.pdf"), 
       width = unit(6, "cm"), height = unit(9, "cm"))

# Extended Data Figure 8E - COSMIC signature SBS42
x <- cbind(sigs.rank3, known_signatures[, match("SBS42", colnames(known_signatures)), drop=F])
p4 <- plot_96_profile(x[, 4, drop=F], condensed = T)
p4 <- p4 + facet_grid(sample ~ substitution, scales = "free_y") +
  coord_cartesian(ylim = NULL) + 
  scale_y_continuous(breaks = waiver())
print(p4)
ggsave(p4, filename = file.path("..", "figures", "Ext_Data_Fig_8E.pdf"), 
       width = unit(6, "cm"), height = unit(2, "cm"))

sessionInfo()
