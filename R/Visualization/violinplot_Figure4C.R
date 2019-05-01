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

dirmain="XXXXXXXXX"
bodyweight <- read.csv(file=paste0(dirmain,"data/inVivo/bodyweight_mice.csv"),header=TRUE)
bodyweight$group <- paste(bodyweight$cell_line,bodyweight$time_point)

vehicule_time <- wilcox.test(bodyweight[which(bodyweight$time_point %in% "Day 0" & 
                                           bodyweight$cell_line %in% "1_Vehicle"), 
                                   "body_weight"], 
                        bodyweight[which(bodyweight$time_point %in% "Day 36" & 
                                           bodyweight$cell_line %in% "1_Vehicle"), 
                                   "body_weight"]) 

MS19_43_time <- wilcox.test(bodyweight[which(bodyweight$time_point %in% "Day 0" & 
                                                bodyweight$cell_line %in% "2_MS19-43"), 
                                        "body_weight"], 
                             bodyweight[which(bodyweight$time_point %in% "Day 36" & 
                                                bodyweight$cell_line %in% "2_MS19-43"), 
                                        "body_weight"]) 

day0 <- wilcox.test(bodyweight[which(bodyweight$time_point %in% "Day 0" & 
                                                bodyweight$cell_line %in% "1_Vehicle"), 
                                        "body_weight"], 
                             bodyweight[which(bodyweight$time_point %in% "Day 0" & 
                                                bodyweight$cell_line %in% "2_MS19-43"), 
                                        "body_weight"]) 

day36 <- wilcox.test(bodyweight[which(bodyweight$time_point %in% "Day 36" & 
                                       bodyweight$cell_line %in% "1_Vehicle"), 
                               "body_weight"], 
                    bodyweight[which(bodyweight$time_point %in% "Day 36" & 
                                       bodyweight$cell_line %in% "2_MS19-43"), 
                               "body_weight"]) 


p12 <- ggplot(data = bodyweight,
              mapping = aes(x = group, y = body_weight)) +
  geom_flat_violin(scale = "count", trim = FALSE, aes(fill = time_point)) +
  geom_boxplot(width=0.2,color="black") +
  geom_dotplot(binaxis = "y", dotsize = 0.5, stackdir = "down", binwidth = 0.1,
               position = position_nudge(-0.025)) +
  geom_line(aes(group=ID), colour="grey", linetype="11") +
  theme_light()+
  theme( axis.text.x = element_blank()) +  
  scale_fill_manual(values=brewer.pal(4, "Set1")[3:4])+
  labs(x = "    Vehicule                     MS19-43", 
       y = "body weight (gr)")+ 
  ylim(17, 24) + 
  geom_segment(aes(x=3, y=23, xend=4, yend=23), size=0.7) +
  geom_segment(aes(x=3, y=22.5, xend=3, yend=23), size=0.7) +
  geom_segment(aes(x=4, y=22.5, xend=4, yend=23), size=0.7) +
  geom_text(x=3.5, y=23.5, size=4.0,  
            label=paste0(format(MS19_43_time$p.value, 
                                digits=3, scientific = TRUE)))+
  geom_segment(aes(x=1, y=23, xend=2, yend=23), size=0.7) +
  geom_segment(aes(x=1, y=22.5, xend=1, yend=23), size=0.7) +
  geom_segment(aes(x=2, y=22.5, xend=2, yend=23), size=0.7) +
  geom_text(x=1.5, y=23.5, size=4.0,label=paste0("NS")) +
  guides(fill=guide_legend(""))
p12
