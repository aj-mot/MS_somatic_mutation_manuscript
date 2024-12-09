#!/usr/bin/env Rscript

# Single-cell WGS QC metrics

library(ggplot2)
library(plyr)
library(lmerTest)

# Extended Data Figure 4A - sSNVs/neuron vs Power ratio
MS.sample.info <- read.csv(file.path("..", "data", "MS_sSNV_burden_data.csv"))
MS.sample.info$Tissue <- factor(MS.sample.info$Tissue, levels=c("CTL", "MS normal", "MS lesion"), ordered=F)

p1 <- ggplot(MS.sample.info, 
            aes(eff.sens, burden, colour=Tissue)) +
  geom_point() +
  ylab("sSNVs/neuron") + xlab("Power ratio") +
  theme_classic() + 
  theme(legend.position="right") +
  theme(legend.title=element_blank())
print(p1)
ggsave(p1, filename = file.path("..", "figures", "Ext_Data_Fig_4A.pdf"), 
       width = unit(6, "cm"), height = unit(4, "cm"))

# Extended Data Figure 4B - sSNVs/neuron vs MAPD
# MAPD
mapd.stats <- read.csv(file.path("..", "data", "mapd_data.csv"))
MS.mapd <- merge(MS.sample.info, mapd.stats)

p2 <- ggplot(MS.mapd, 
            aes(MAPD, burden, colour=Tissue)) +
  geom_point() +
  ylab("sSNVs/neuron") + xlab("MAPD") +
  theme_classic() + 
  theme(legend.position="right") +
  theme(legend.title=element_blank())
print(p2)
ggsave(p2, filename = file.path("..", "figures", "Ext_Data_Fig_4B.pdf"), 
       width = unit(6, "cm"), height = unit(4, "cm"))


# Extended Data Figure 4C - sSNVs/neuron vs CoV
p3 <- ggplot(MS.mapd, 
            aes(CoV, burden, colour=Tissue)) +
  geom_point() +
  ylab("sSNVs/neuron") + xlab("CoV") +
  theme_classic() + 
  theme(legend.position="right") +
  theme(legend.title=element_blank())
print(p3)
ggsave(p3, filename = file.path("..", "figures", "Ext_Data_Fig_4C.pdf"), 
       width = unit(6, "cm"), height = unit(4, "cm"))


# Extended Data Figure 4D - sSNVs/neuron vs Proportion of sites with allele 
# balance between 0.3 and 0.7
ab.frac.results <- read.csv(file.path("..", "data", "allele_balance_data.csv"))
MS.mapd <- merge(MS.mapd, ab.frac.results)

p4 <- ggplot(MS.mapd, 
            aes(ab.frac, burden, colour=Tissue)) +
  geom_point() +
  ylab("sSNVs/neuron") + xlab("Proportion of sites with allele\nbalance between 0.3 and 0.7") +
  theme_classic() + 
  theme(legend.position="right") +
  theme(legend.title=element_blank())
print(p4)
ggsave(p4, filename = file.path("..", "figures", "Ext_Data_Fig_4D.pdf"), 
       width = unit(6, "cm"), height = unit(4, "cm"))


# Extended Data Figure 4E - Residual sSNVs/neuron vs Post-mortem interval

# Read external control data
luquette.data <- read.csv(file.path("..", "external_data", "luquette_sSNV_burden_data.csv"))
MS.sample.info <- rbind.fill(data.frame(MS.sample.info, Lab="Rubio"), 
                             data.frame(luquette.data, Lab="Walsh"))

# Fit lmer
MS.sample.info$is.MSnormal <- MS.sample.info$Tissue == "MS normal"
MS.sample.info$is.MSlesion <- MS.sample.info$Tissue == "MS lesion"
burden.lmer <- lmer(burden ~ Age * is.MSlesion + is.MSnormal + (1 | Individual), 
                     MS.sample.info)
summary(burden.lmer)

# Compute the fitted sSNV burden for each neuron
MS.sample.info$burden.fitted <- summary(burden.lmer)$coefficients[1, 1] + 
  summary(burden.lmer)$coefficients[2, 1] * MS.sample.info$Age + 
  summary(burden.lmer)$coefficients[4, 1] * MS.sample.info$is.MSnormal +
  (summary(burden.lmer)$coefficients[3, 1] + 
     summary(burden.lmer)$coefficients[5, 1] * MS.sample.info$Age) * MS.sample.info$is.MSlesion

# Compute the residuals of the fitted model
MS.sample.info$burden.residual <- MS.sample.info$burden - MS.sample.info$burden.fitted

# Fit model for residual vs PMI
pmi.residual.lmer <- lmer(burden.residual ~ PMI + (1 | Individual), 
                          data = MS.sample.info)
summary(pmi.residual.lmer)

ms.tissue.colours <- c(`CTL`="steelblue1", `MS normal`="olivedrab3", 
                       `MS lesion`="coral1")

p5 <- ggplot(subset(MS.sample.info, Tissue != "External control"), 
             aes(PMI, burden.residual, fill=Tissue)) + 
  geom_jitter(width=0.7, height=0, colour="black", size=3, shape=21) +
  xlab("Post-mortem interval (hours)") + ylab("Residual sSNVs / neuron (actual - fitted)") +
  scale_fill_manual(labels=names(ms.tissue.colours), values=ms.tissue.colours) +
  geom_hline(yintercept = 0, linetype="dotted", colour="black") +
  geom_abline(intercept = summary(pmi.residual.lmer)$coefficients[1, 1], slope = summary(pmi.residual.lmer)$coefficients[2, 1], colour="black") +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(), legend.position=c(0.85, 0.9), 
        legend.key.size = unit(0.4, 'cm'), 
        legend.text = element_text(size=8)) +
  guides(linetype="none", shape="none", colour="none", fill=guide_legend(override.aes=list(shape=21)))

print(p5)
ggsave(p5, filename = file.path("..", "figures", "Ext_Data_Fig_4E.pdf"), 
       width = unit(7, "cm"), height = unit(5, "cm"))
