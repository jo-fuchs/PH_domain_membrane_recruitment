####Description###############################################################################################################
### Membrane recruitment analysis Part 2
###
###   Goals:
###   - Normalizing line profiles to peaks of cell membrane marker (here: F-actin) to create merge
###   - Read out Membrane to Cytosol ratio of Intensities
###   - Create quality control graphs and summary Figures
###
###   requires data created through: ImageJ macro "Save_profiles.ijm"
###
###   Version 0.2 (18.06.2020)
###   - tidyversed
###   - account for cell sizes
###
###   Joachim Fuchs
###
##################################################################################################################


## Load packages, define functions

library(tidyverse)
library(ggsci)
library(ggbeeswarm)
source("merge_lines.R")


#####  Analysis  #################################################################################################################

# define Directory of Line profiles
folderdir <- dirname(file.choose(new = FALSE))

## Run functions
merge_lines(folderdir)


## add identifiers, set order of factors

## Load data
Stats <- read_delim(file.path(folderdir, "Results", "Membrane_recruitment_results.txt"), delim = "\t")
rawData <- read_delim(file.path(folderdir, "Results", "All_peaks.txt"), delim = "\t")


## Load group identifier Table & add to Stats and raw-Data tables
Groups <- read_delim(file.path(folderdir, "Results/Groups.txt"), delim = "\t")

Stats <- left_join(Stats, Groups, by = "Image")
rawData <- left_join(rawData, Groups, by = "Image")



### add responding/non-responding column by Stats$MemCentRatio_PH>=0.95
Stats$Responding[Stats$MemMinRatio >= 1] <- "Responding"
Stats$Responding[Stats$MemMinRatio < 1] <- "Not Responding"
Total <- merge(rawData, Stats)

#optional removing 2min Insulin which is only present in one experiment
Total <- Total %>% filter(Stimulation != "Insulin2")
Stats <- Stats %>% filter(Stimulation != "Insulin2")

# Set as Factors & set correct order
Stats$Stimulation <- factor(Stats$Stimulation, levels = c("Starve", "Insulin"))
Total$Stimulation <- factor(Total$Stimulation, levels = c("Starve", "Insulin"))
Stats$Responding <- factor(Stats$Responding, levels = c("Responding", "Not Responding"))
Total$Responding <- factor(Total$Responding, levels = c("Responding", "Not Responding"))


## Calculate fraction of responding cells
# Responders defined as Membrane/MinCell ration >=1
Percentage <- data.frame(Responding_fration = NA)
Responding_fration <- "Responding_fraction"
Percentage <- Total %>%
  group_by(Stimulation, Inhibition) %>%
  summarise(!!Responding_fration := sum(Responding == "Responding") / length(Responding))
rm(Responding_fration)

## Add responding percentage to each line-profile
Total <- merge(Total, Percentage)

## Change fraction to (1-Responding) in Non-responders
Total$Responding_fraction[Total$Responding == "Not Responding"] <- 1 - Total$Responding_fraction[Total$Responding == "Not Responding"]


# Summarize to Fields of view
MeanRatio <- "MeanRatio"
ByImage <- Total %>%
  group_by(Image, Experiment, Stimulation, Inhibition) %>%
  summarise(!!MeanRatio := mean(MemMinRatio))

ByImage$Stimulation <- factor(ByImage$Stimulation, levels = c("Starve", "Insulin"))


      
###  Plots  #######################################################################################################################

resultsdir <- file.path(folderdir, "Results")


### plot overlayed lines for Actin raw
plot_Actin <- ggplot(Total, aes(x = X, y = Actin)) +
  geom_line(aes(group = ROI), color = "grey", alpha = 0.3) +
 
  # Axes and Title
  xlab("Line position [µm]") +
  ylab("Actin Intensity\n") +
  ggtitle("Peak alignment") +

  # Theme
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text = element_text(face = "bold"), axis.title = element_text(size = rel(1.15)),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_startrek() +
  scale_fill_startrek()


plot_Actin

ggsave("ActinPeaks.jpg",
  plot = plot_Actin, device = "jpeg", path = resultsdir,
  scale = 1, width = 20, height = 15, units = "cm"
)
ggsave("ActinPeaks.pdf",
  plot = plot_Actin, device = "pdf", path = resultsdir,
  scale = 1, width = 20, height = 15, units = "cm"
)



### plot overlayed lines for Actin
plot_Actin <- ggplot(Total, aes(x = X_norm, y = Actin)) +
  geom_line(aes(group = ROI), color = "grey", alpha = 0.3) +
  coord_cartesian(xlim = c(0.3, 2.7)) +

  # Axes and Title
  xlab("Peak") +
  ylab("Actin Intensity\n") +
  ggtitle("Peak alignment") +

  # Theme
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text = element_text(face = "bold"), axis.title = element_text(size = rel(1.15)),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_startrek() +
  scale_fill_startrek()


plot_Actin

ggsave("ActinPeaks_Aligned.jpg",
  plot = plot_Actin, device = "jpeg", path = resultsdir,
  scale = 1, width = 20, height = 15, units = "cm"
)
ggsave("ActinPeaks_Aligned.pdf",
  plot = plot_Actin, device = "pdf", path = resultsdir,
  scale = 1, width = 20, height = 15, units = "cm"
)



# PH domain, Groups separate, PH normalized
plot_PH <- ggplot(Total, aes(x = X_norm, y = PH_norm)) +
  geom_line(aes(group = ROI), color = "grey", alpha = 0.2) +
  coord_cartesian(xlim = c(0.6, 2.4)) +
  scale_y_continuous(breaks = NULL) +
  
  facet_wrap(~Stimulation*Inhibition, ncol = 2) +
  geom_smooth(aes(fill = Stimulation), method = "loess", size = 0.5, color = "black", span = 0.1) +

  # Axes and Title
  xlab("PH line profiles aligned to Actin peaks") +
  ylab("Normalized PH-domain Intensity\n") +

  # Theme
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text = element_text(face = "bold"), axis.title = element_text(size = rel(1.15)),
    axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "grey90", color = NA)
  ) +
  scale_color_startrek() +
  scale_fill_startrek()


plot_PH

ggsave("PHPeaks_Aligned.jpg",
  plot = plot_PH, device = "jpeg", path = resultsdir,
  scale = 1, width = 20, height = 15, units = "cm"
)
ggsave("PHPeaks_Aligned.pdf",
  plot = plot_PH, device = "pdf", path = resultsdir,
  scale = 1, width = 20, height = 15, units = "cm"
)



# PH domain, Groups separate, PH normalized, Responders vs non
plot_PH <- ggplot(Total, aes(x = X_norm, y = PH_norm)) +
  geom_line(aes(group = ROI, color = Responding), alpha = 0.05) +
  coord_cartesian(xlim = c(0.3, 2.7), ylim = c(0, 1)) +
  facet_wrap(~Stimulation*Inhibition, ncol = 2) +
  geom_smooth(aes(fill = Responding, alpha = Responding_fraction), method = "loess", size = 0.5, color = "black", span = 0.1) +

  # Axes and Title
  xlab("PH line profiles aligned to Actin peaks") +
  ylab("Normalized PH-domain Intensity") +

  # Theme
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.title = element_text(size = rel(1.15)),
    axis.text = element_blank(), plot.title = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "grey90", color = NA)
  ) +
  scale_color_startrek() +
  scale_fill_startrek()


plot_PH

ggsave("PHPeaks_Responding.jpg",
  plot = plot_PH, device = "jpeg", path = resultsdir,
  scale = 1, width = 20, height = 15, units = "cm"
)
ggsave("PHPeaks_Responding.pdf",
  plot = plot_PH, device = "pdf", path = resultsdir,
  scale = 1, width = 20, height = 15, units = "cm"
)



## stats plots Mem to MinCell Ratio
plot2 <- ggplot(Stats, aes(x = Stimulation, y = abs(MemMinRatio))) +
  # geom_boxplot(outlier.alpha = 0)+
  geom_quasirandom(aes(color = Responding, alpha = MeanCenter_PH), size = 2, width = 0.3, shape = 16) + ## shape 16 = point without borders
  facet_wrap(~Inhibition, strip.position = "bottom") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5) + 
  stat_summary(fun.data = "mean_se", geom = "errorbar", color = "black", width = 0.3) + ## error bars (mean_se = SEM)
  coord_cartesian(ylim = c(0, 3)) +

  # Axes and Title
  ggtitle("BTK-PH membrane localization as PIP3-sensor") +
  xlab("") +
  ylab("Membrane to Minima in Cell ratio\n") +
  labs(alpha = "Expression level", color = "") +

  # Theme
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = rel(1.15)),
    axis.ticks.y = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, size = rel(1.3), face = "bold")
  ) +
  scale_color_startrek() +
  scale_fill_startrek()


plot2

ggsave("MemMinRatio.jpg",
  plot = plot2, device = "jpeg", path = resultsdir,
  scale = 1, width = 15, height = 15, units = "cm"
)
ggsave("MemMinRatio.pdf",
  plot = plot2, device = "pdf", path = resultsdir,
  scale = 1, width = 12, height = 15, units = "cm"
)



# stats plot Mem to Center Ratio
plot3 <- ggplot(Stats, aes(x = Stimulation, y = abs(MemCentRatio_PH))) +
  # geom_boxplot(outlier.alpha = 0)+
  geom_quasirandom(aes(color = Responding, alpha = MeanCenter_PH), size = 2, width = 0.3, shape = 16) + ## shape 16 = point without borders
  facet_wrap(~Inhibition, strip.position = "bottom") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5) + 
  stat_summary(fun.data = "mean_se", geom = "errorbar", color = "black", width = 0.3) + ## error bars (mean_se = SEM)
  coord_cartesian(ylim = c(0, 3)) +
  
  # Axes and Title
  ggtitle("BTK-PH membrane localization as PIP3-sensor") +
  xlab("") +
  ylab("Membrane to Center of Cell ratio\n") +
  labs(alpha = "Expression level", color = "") +
  
  # Theme
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = rel(1.15)),
    axis.ticks.y = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, size = rel(1.3), face = "bold")
  ) +
  scale_color_startrek() +
  scale_fill_startrek()



plot3

ggsave("MemCentRatio.jpg",
  plot = plot3, device = "jpeg", path = resultsdir,
  scale = 1, width = 15, height = 15, units = "cm"
)
ggsave("MemCentRatio.pdf",
  plot = plot3, device = "pdf", path = resultsdir,
  scale = 1, width = 12, height = 15, units = "cm"
)



## ByImage plots Mem to MinCell Ratio
plot2 <- ggplot(ByImage, aes(x = Stimulation, y = abs(MeanRatio))) +
  # geom_boxplot(outlier.alpha = 0)+
  geom_quasirandom(aes(color = Experiment), size = 2, width = 0.3, shape = 16) + ## shape 16 = point without borders
  facet_grid(~Inhibition) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5) + 
  stat_summary(fun.data = "mean_se", geom = "errorbar", color = "black", width = 0.3) + ## error bars (mean_se = SEM)
  coord_cartesian(ylim = c(0, 3)) +
  
  # Axes and Title
  ggtitle("BTK-PH membrane localization as PIP3-sensor") +
  xlab("") +
  ylab("Membrane to Minima in Cell ratio\n") +
  labs(alpha = "Expression level", color = "") +

  # Theme
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.title = element_text(size = rel(1.15)),
    axis.ticks.y = element_line(color = "black"),
    strip.background = element_rect(fill = "grey90", color = NA),
    plot.title = element_text(hjust = 0.5, size = rel(1.3), face = "bold")
  ) +
  scale_color_startrek() +
  scale_fill_startrek()


plot2

ggsave("MemMinRatio.jpg",
  plot = plot2, device = "jpeg", path = resultsdir,
  scale = 1, width = 15, height = 15, units = "cm"
)
ggsave("MemMinRatio.pdf",
  plot = plot2, device = "pdf", path = resultsdir,
  scale = 1, width = 15, height = 15, units = "cm"
)


