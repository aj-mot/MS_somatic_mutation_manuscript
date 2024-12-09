#!/usr/bin/env Rscript

# Analysis of genome-wide sSNV burden

library(ggplot2)
library(plyr)
library(lmerTest)
library(RColorBrewer)

# sSNV burden data
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

# Order Tissue
tissue.order <- c("External control", "Control", "MS normal", "MS lesion")
MS.sample.info$Tissue <- factor(MS.sample.info$Tissue, 
                                levels = tissue.order, ordered = T)

ms.tissue.colours <- c(`External control`="blueviolet", `Control`="steelblue1", 
                       `MS normal`="olivedrab3", `MS lesion`="coral1")


# Fit lmer
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

# Compute the age-adjusted burden - subtract age-trend fitted to controls
MS.sample.info$age.adjusted.burden <- MS.sample.info$burden - summary(burden.lmer)$coefficients[1, 1] - summary(burden.lmer)$coefficients[2, 1] * MS.sample.info$Age

# Figure 2A - sSNVs/neuron vs Age

# Compute fitted lines
summary(burden.lmer)$coefficients

min.age.ext.ctl <- min(MS.sample.info$Age[MS.sample.info$Tissue == "External control"])
max.age.ext.ctl <- max(MS.sample.info$Age[MS.sample.info$Tissue == "External control"])
min.age.ctl <- min(MS.sample.info$Age[MS.sample.info$Tissue == "Control"])
max.age.ctl <- max(MS.sample.info$Age[MS.sample.info$Tissue == "Control"])
min.age.msnormal <- min(MS.sample.info$Age[MS.sample.info$Tissue == "MS normal"])
max.age.msnormal <- max(MS.sample.info$Age[MS.sample.info$Tissue == "MS normal"])
min.age.mslesion <- min(MS.sample.info$Age[MS.sample.info$Tissue == "MS lesion"])
max.age.mslesion <- max(MS.sample.info$Age[MS.sample.info$Tissue == "MS lesion"])

min.age.ctl2 <- min(min.age.ext.ctl, min.age.ctl)
max.age.ctl2 <- max(max.age.ext.ctl, max.age.ctl)

ctl.line <- data.frame(Age=c(min.age.ctl2, max.age.ctl2), 
                        burden=c(summary(burden.lmer)$coefficients[1, 1] + min.age.ctl2*summary(burden.lmer)$coefficients[2, 1], 
                                 summary(burden.lmer)$coefficients[1, 1] + max.age.ctl2*summary(burden.lmer)$coefficients[2, 1]),
                        Tissue="Control", 
                        lty="1"
)
ext.ctl.line <- ctl.line
ext.ctl.line$lty <- "2"
ext.ctl.line$Tissue <- "External control"
msnormal.line <- data.frame(Age=c(min.age.msnormal, max.age.msnormal), 
                             burden=c(summary(burden.lmer)$coefficients[1, 1] + summary(burden.lmer)$coefficients[4, 1] + min.age.msnormal*summary(burden.lmer)$coefficients[2, 1], 
                                      summary(burden.lmer)$coefficients[1, 1] + summary(burden.lmer)$coefficients[4, 1] + max.age.msnormal*summary(burden.lmer)$coefficients[2, 1]),
                             Tissue="MS normal", 
                             lty="1"
)
mslesion.line <- data.frame(Age=c(min.age.mslesion, max.age.mslesion), 
                             burden=c(summary(burden.lmer)$coefficients[1, 1] + summary(burden.lmer)$coefficients[3, 1] + min.age.mslesion*(summary(burden.lmer)$coefficients[2, 1] + summary(burden.lmer)$coefficients[5, 1]), 
                                      summary(burden.lmer)$coefficients[1, 1] + summary(burden.lmer)$coefficients[3, 1] + max.age.mslesion*(summary(burden.lmer)$coefficients[2, 1] + summary(burden.lmer)$coefficients[5, 1])),
                             Tissue="MS lesion", 
                             lty="1"
)

p1 <- ggplot(MS.sample.info, 
             aes(Age, burden, fill=Tissue)) + 
  geom_jitter(width=0.7, height=0, colour="black", size=3, shape=21) +
  geom_line(aes(Age, burden, colour=Tissue, linetype=lty), data=rbind(ctl.line, ext.ctl.line, msnormal.line, mslesion.line)) +
  xlab("Age (years)") + ylab("sSNVs / neuron") +
  scale_fill_manual(labels=names(ms.tissue.colours), values=ms.tissue.colours) +
  scale_colour_manual(labels=names(ms.tissue.colours), values=ms.tissue.colours) +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(), legend.position=c(0.125, 0.9), 
        legend.key.size = unit(0.4, 'cm'), 
        legend.text = element_text(size=8)) +
  guides(linetype="none", shape="none", colour="none", fill=guide_legend(override.aes=list(shape=21)))

print(p1)

ggsave(p1, filename = file.path("..", "figures", "Fig_2A.pdf"), 
       width = unit(7, "cm"), height = unit(5, "cm"))


# Figure 2B - Age-adjusted sSNVs/neuron vs Tissue

# Fit model to get means and error bars for age-adjusted plot
age.adjusted.burden.lmer <- lmer(age.adjusted.burden ~ is.Control + is.MSnormal + is.MSlesion + (1 | Individual), 
                                 subset(MS.sample.info, WGA == "PTA"))
summary(age.adjusted.burden.lmer)

age.adjusted.burden.ci <- confint(age.adjusted.burden.lmer)

age.adjusted.burden.mean <- summary(age.adjusted.burden.lmer)$coefficients[, 1]
age.adjusted.burden.lb <- age.adjusted.burden.ci[3:6, 1]
age.adjusted.burden.ub <- age.adjusted.burden.ci[3:6, 2]

bar.width <- 0.4
age.adjusted.burden.bars <- data.frame(Tissue=c("External control", "Control", "MS normal", "MS lesion"), 
                                       x.start=1:4-bar.width, x.end=1:4+bar.width,
                                       Mean=age.adjusted.burden.mean, LB=age.adjusted.burden.lb, UB=age.adjusted.burden.ub)

# Coordinates for bars showing comparisons with statistically significant differences
h1 <- 2720
h2 <- 2870
h3 <- 3050
hw <- 80
sig.segments <- data.frame(   x = c(3.1, 2.1, 1,   1.1, 1.1, 2.9, 3.1, 3.9, 2.1, 4, 1, 4.1), 
                              xend = c(3.9, 4, 4.1, 2.9, 1.1, 2.9, 3.1, 3.9, 2.1, 4, 1, 4.1), 
                              y = c(h1, h2, h3, h1, h1,    h1,    h1,    h1,    h2,    h2,    h3,    h3), 
                              yend = c(h1, h2, h3, h1, h1-hw, h1-hw, h1-hw, h1-hw, h2-hw, h2-hw, h3-hw, h3-hw))
sig.labels <- data.frame(x=c(2.5, 3, 3.5, 1.9), 
                         y=c(h3+30, h2+30, h1+30, h1+30), 
                         label=c("***", "***", "***", "**"))

p2 <- ggplot(MS.sample.info, 
           aes(Tissue, age.adjusted.burden, fill = Tissue)) + 
  geom_jitter(shape = 21, width = 0.25, height=0, alpha=1, colour="black", size=3) +
  geom_segment(aes(x=x.start, y=Mean, xend=x.end, yend=Mean), age.adjusted.burden.bars, lwd=0.35, linetype="solid", colour="black") +
  geom_segment(aes(x=x.start, y=LB, xend=x.end, yend=LB), age.adjusted.burden.bars, lwd=0.15, linetype="solid", colour="black") +
  geom_segment(aes(x=x.start, y=UB, xend=x.end, yend=UB), age.adjusted.burden.bars, lwd=0.15, linetype="solid", colour="black") +
  geom_hline(yintercept = 0, linetype="dotted", colour="black") +
  geom_segment(aes(x, y, xend=xend, yend=yend), sig.segments, inherit.aes = F) +
  geom_text(aes(x, y, label=label), sig.labels, size=6, inherit.aes = F) +
  ylim(c(-420, 3080)) +
  ylab("sSNVs / neuron in excess of expectation for age") +
  scale_fill_manual(labels=names(ms.tissue.colours), values=ms.tissue.colours) +
  guides(fill="none") +
  theme_classic(base_size = 12) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme(axis.title.y = element_text(size=10)) +
  scale_x_discrete(labels=c("External control\n(n=33)", "Control\n(n=12)", "MS normal\n(n=29)", "MS lesion\n(n=32)"))

print(p2)

ggsave(p, filename = file.path("..", "figures", "Fig_2B.pdf"), 
       width = unit(5, "cm"), height = unit(5, "cm"))


# Figure 2C - sSNVs/neuron in MS cases

# Plot MS cases vs Tissue faceted by Individual
ms.tissue.shapes <- c(`MS normal`=21, `MS lesion`=24)
ind.sig.level <- c("MS10"="", "MS12"="", "MS13"="*", "MS14"="*", "MS4"="", "MS6"="*", "MS7"="", "MS8"="", "MS9"="", "MS11"="")

p3 <- ggplot(subset(MS.sample.info, Disease == "MS"), 
             aes(Individual, burden)) + 
  geom_jitter(aes(fill=Tissue, shape=Tissue), width = 0.125, height=0, alpha=1, colour="black", size=3) +
  stat_summary(aes(lty=Tissue), 
               fun = mean, 
               geom="crossbar", 
               position = "identity", lwd=0.25) +
  facet_grid(. ~ Individual, scale="free_x", space="free_x") +
  ylab("sSNVs / neuron") +
  scale_fill_manual(labels=names(ms.tissue.colours)[3:4], values=ms.tissue.colours[3:4]) +
  scale_shape_manual(values=ms.tissue.shapes) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(size=18, hjust=0.5, vjust=0.0)) + 
  scale_x_discrete(labels=ind.sig.level) +
  theme(axis.title.y = element_text(size=10)) + 
  theme(axis.title.x = element_blank()) + 
  theme(legend.title = element_blank(), legend.position=c(0.085, 0.9), legend.text = element_text(size=8)) +
  theme(strip.text.x = element_text(size = 8))

print(p3)
ggsave(p3, filename = file.path("..", "figures", "Fig_2C.pdf"), 
       width = unit(7, "cm"), height = unit(5, "cm"))


# Extended Data Figure 5 - sSNVs/neuron in controls

# Plot controls ordered by age
p4 <- ggplot(subset(MS.sample.info, Disease == "CTL"), 
            aes(Individual, burden)) + 
  geom_jitter(aes(fill=Tissue), shape=21, width = 0.075, height=0, alpha=1, colour="black", size=3) +
  stat_summary(fun = mean, 
               geom="crossbar", 
               position = "identity", lwd=0.25) +
  theme_classic(base_size = 12) +
  labs(x="", y="sSNVs / neuron") +
  scale_fill_manual(labels=names(ms.tissue.colours), values=ms.tissue.colours) +
  theme(axis.text.x = element_text(size=9, angle = 90, hjust=0.95, vjust=0.2)) + 
  theme(axis.title.y = element_text(size=10)) + 
  theme(legend.title = element_blank(), legend.position=c(0.132, 0.9), legend.text = element_text(size=8))
print(p4)

ggsave(p4, filename = file.path("..", "figures", "Ext_Data_Fig_5.pdf"), 
       width = unit(6, "cm"), height = unit(5, "cm"))
