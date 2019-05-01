## set working directory
getwd()

## call required packages
library(Rcmdr)
library(ggplot2)
library(RColorBrewer)

pkgs <- c("tidyverse", "ggthemes", "skimr")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

## load source code
source("geom_flat_violin.R")

dirmain="XXXXXXX"
data <- read.csv(file=paste0(dirmain,"data/inVivo/DataResult_Anqi_Cal.csv"),header=TRUE)
plasma <- data[which(data$tissue %in% "Plasma"),"value"]
tumor <- data[which(data$tissue %in% "Tumor"),"value"]


MS1943 <- wilcox.test(plasma,tumor,paired=TRUE) 


colors2 <- brewer.pal(4, "Set2")[3:4]
p2_2<- ggplot(data, aes(tissue, value)) +
  geom_flat_violin(scale = "count", trim = FALSE, aes(fill = tissue)) +
  scale_fill_manual(values=brewer.pal(4, "Set2")[3:4])+
  geom_boxplot(width=0.3, size=1, fatten=1, colour="black") +
  geom_point(colour="red", size=2, alpha=0.5) +
  geom_line(aes(group=Sample.ID.in.Raw.Data), colour="grey", linetype="11") +
  labs(x = "", y = "MS1943 \n Concentration (nM)") +
  geom_text(aes(x=1.5, y=600000),
            label= paste0("P-value= ",
                          format(MS1943$p.value, digits=3, scientific = TRUE))) +
  geom_segment(aes(x=1, y=510000, xend=2, yend=510000), size=0.7) +
  geom_segment(aes(x=1, y=430000, xend=1, yend=510000), size=0.7) +
  geom_segment(aes(x=2, y=430000, xend=2, yend=510000), size=0.7) +
  ylim(1500, 1400000) +
  theme_light() + 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))
p2_2

