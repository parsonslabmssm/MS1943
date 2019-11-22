library("ggplot2")

black.italic.12.text <- element_text(face = "italic", color = "black", size = 12)

Dir_project=""

#### Fig 2a
CellViability <- read.csv(paste0(Dir_project,"data/data_fig2a.csv"),header=TRUE)
CellViability$dose <- as.factor(CellViability$dose)

## stat CellViability C24
c24 <- CellViability[which(CellViability$treatment %in% "C24"),]
c24_stat1 <- t.test(c24[which(c24$dose %in% "0 μM"), "value"], 
                    c24[which(c24$dose %in% "0.625 μM"), "value"],
                    alternative = "two.sided",paired = FALSE) 

c24_stat2 <- t.test(c24[which(c24$dose %in% "0 μM"), "value"], 
                    c24[which(c24$dose %in% "1.25 μM"), "value"],
                    alternative = "two.sided",paired = FALSE)
c24_stat3 <- t.test(c24[which(c24$dose %in% "0 μM"), "value"], 
                    c24[which(c24$dose %in% "2.5 μM"), "value"],
                    alternative = "two.sided",paired = FALSE) 
c24_stat4 <- t.test(c24[which(c24$dose %in% "0 μM"), "value"], 
                    c24[which(c24$dose %in% "5 μM"), "value"],
                    alternative = "two.sided",paired = FALSE) 

## stat CellViability CPI-1205
CPI <- CellViability[which(CellViability$treatment %in% "CPI-1205"),]
CPI_stat1 <- t.test(CPI[which(CPI$dose %in% "0 μM"), "value"], 
                    CPI[which(CPI$dose %in% "0.625 μM"), "value"],
                    alternative = "two.sided",paired = FALSE) 

CPI_stat2 <- t.test(CPI[which(CPI$dose %in% "0 μM"), "value"], 
                    CPI[which(CPI$dose %in% "1.25 μM"), "value"],
                    alternative = "two.sided",paired = FALSE)
CPI_stat3 <- t.test(CPI[which(CPI$dose %in% "0 μM"), "value"], 
                    CPI[which(CPI$dose %in% "2.5 μM"), "value"],
                    alternative = "two.sided",paired = FALSE) 
CPI_stat4 <- t.test(CPI[which(CPI$dose %in% "0 μM"), "value"], 
                    CPI[which(CPI$dose %in% "5 μM"), "value"],
                    alternative = "two.sided",paired = FALSE) 

## stat CellViability EPZ6438
EPZ <- CellViability[which(CellViability$treatment %in% "EPZ6438"),]
EPZ_stat1 <- t.test(EPZ[which(EPZ$dose %in% "0 μM"), "value"], 
                    EPZ[which(EPZ$dose %in% "0.625 μM"), "value"],
                    alternative = "two.sided",paired = FALSE) 

EPZ_stat2 <- t.test(EPZ[which(EPZ$dose %in% "0 μM"), "value"], 
                    EPZ[which(EPZ$dose %in% "1.25 μM"), "value"],
                    alternative = "two.sided",paired = FALSE)
EPZ_stat3 <- t.test(EPZ[which(EPZ$dose %in% "0 μM"), "value"], 
                    EPZ[which(EPZ$dose %in% "2.5 μM"), "value"],
                    alternative = "two.sided",paired = FALSE) 
EPZ_stat4 <- t.test(EPZ[which(EPZ$dose %in% "0 μM"), "value"], 
                    EPZ[which(EPZ$dose %in% "5 μM"), "value"],
                    alternative = "two.sided",paired = FALSE) 

## stat CellViability GSK126
GSK <- CellViability[which(CellViability$treatment %in% "GSK126"),]
GSK_stat1 <- t.test(GSK[which(GSK$dose %in% "0 μM"), "value"], 
                    GSK[which(GSK$dose %in% "0.625 μM"), "value"],
                    alternative = "two.sided",paired = FALSE) 

GSK_stat2 <- t.test(GSK[which(GSK$dose %in% "0 μM"), "value"], 
                    GSK[which(GSK$dose %in% "1.25 μM"), "value"],
                    alternative = "two.sided",paired = FALSE)
GSK_stat3 <- t.test(GSK[which(GSK$dose %in% "0 μM"), "value"], 
                    GSK[which(GSK$dose %in% "2.5 μM"), "value"],
                    alternative = "two.sided",paired = FALSE) 
GSK_stat4 <- t.test(GSK[which(GSK$dose %in% "0 μM"), "value"], 
                    GSK[which(GSK$dose %in% "5 μM"), "value"],
                    alternative = "two.sided",paired = FALSE) 

## stat CellViability MS1943
MS <- CellViability[which(CellViability$treatment %in% "MS1943"),]
MS_stat1 <- t.test(MS[which(MS$dose %in% "0 μM"), "value"], 
                   MS[which(MS$dose %in% "0.625 μM"), "value"],
                   alternative = "two.sided",paired = FALSE) 

MS_stat2 <- t.test(MS[which(MS$dose %in% "0 μM"), "value"], 
                   MS[which(MS$dose %in% "1.25 μM"), "value"],
                   alternative = "two.sided",paired = FALSE)
MS_stat3 <- t.test(MS[which(MS$dose %in% "0 μM"), "value"], 
                   MS[which(MS$dose %in% "2.5 μM"), "value"],
                   alternative = "two.sided",paired = FALSE) 
MS_stat4 <- t.test(MS[which(MS$dose %in% "0 μM"), "value"], 
                   MS[which(MS$dose %in% "5 μM"), "value"],
                   alternative = "two.sided",paired = FALSE) 

#"#d9f0d3",
p <- ggplot(CellViability, aes(x=treatment, y=value)) + 
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.7), dotsize = 0.6,
               aes(fill=dose)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
               geom="pointrange", color="black",
               mapping = aes(group = dose),
               position=position_dodge(0.7)) +
  ylim(0, 175) +
  scale_fill_manual(values=c("#762a83","#af8dc3", "#e7d4e8",
                             "#7fbf7b","#1b7837"))+
  theme_linedraw(base_size=12) +
  theme(axis.text.x = element_text(color = "black", size = 12),
        plot.title = element_text(hjust = 0.5)) +
  xlab("") + ylab("Normalized Cell Viability (%)") +labs(fill = "") +
  ggtitle("MDA-MB-468") 

pc24 <- p +  geom_segment(aes(x=0.7, y=106, xend=0.85, yend=106), size=0.7) +
  geom_segment(aes(x=0.7, y=102, xend=0.7, yend=106), size=0.7) +
  geom_segment(aes(x=0.85, y=102, xend=0.85, yend=106), size=0.7) +
  geom_text(x=0.8, y=111, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(c24_stat1$p.value, 
                                digits=2, scientific = FALSE))) +
  geom_segment(aes(x=0.7, y=120, xend=1, yend=120), size=0.7) +
  geom_segment(aes(x=0.7, y=116, xend=0.7, yend=120), size=0.7) +
  geom_segment(aes(x=1, y=116, xend=1, yend=120), size=0.7) +
  geom_text(x=0.85, y=127, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(c24_stat2$p.value, 
                                digits=2, scientific = FALSE)))+
  geom_segment(aes(x=0.7, y=152, xend=1.15, yend=152), size=0.7) +
  geom_segment(aes(x=0.7, y=148, xend=0.7, yend=152), size=0.7) +
  geom_segment(aes(x=1.15, y=148, xend=1.15, yend=152), size=0.7) +
  geom_text(x=0.9, y=159, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(c24_stat3$p.value, 
                                digits=2, scientific = FALSE)))+
  geom_segment(aes(x=0.7, y=167, xend=1.3, yend=167), size=0.7) +
  geom_segment(aes(x=0.7, y=162, xend=0.7, yend=167), size=0.7) +
  geom_segment(aes(x=1.3, y=162, xend=1.3, yend=167), size=0.7) +
  geom_text(x=0.9, y=173, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(c24_stat4$p.value, 
                                digits=2, scientific = FALSE)))

#  CPI-1205
pCPI <- pc24 +geom_segment(aes(x=1.7, y=111, xend=1.85, yend=111), size=0.7) +
  geom_segment(aes(x=1.7, y=107, xend=1.7, yend=111), size=0.7) +
  geom_segment(aes(x=1.85, y=107, xend=1.85, yend=111), size=0.7) +
  geom_text(x=1.8, y=116, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(CPI_stat1$p.value, 
                                digits=2, scientific = FALSE))) +
  geom_segment(aes(x=1.7, y=125, xend=2, yend=125), size=0.7) +
  geom_segment(aes(x=1.7, y=121, xend=1.7, yend=125), size=0.7) +
  geom_segment(aes(x=2, y=121, xend=2, yend=125), size=0.7) +
  geom_text(x=1.85, y=132, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(CPI_stat2$p.value, 
                                digits=2, scientific = FALSE)))+
  geom_segment(aes(x=1.7, y=139, xend=2.15, yend=139), size=0.7) +
  geom_segment(aes(x=1.7, y=135, xend=1.7, yend=139), size=0.7) +
  geom_segment(aes(x=2.15, y=135, xend=2.15, yend=139), size=0.7) +
  geom_text(x=1.9, y=146, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(CPI_stat3$p.value, 
                                digits=2, scientific = FALSE))) +
  geom_segment(aes(x=1.7, y=152, xend=2.3, yend=152), size=0.7) +
  geom_segment(aes(x=1.7, y=148, xend=1.7, yend=152), size=0.7) +
  geom_segment(aes(x=2.3, y=148, xend=2.3, yend=152), size=0.7) +
  geom_text(x=1.9, y=159, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(CPI_stat4$p.value, 
                                digits=2, scientific = FALSE)))

#  EPZ6438
pEPZ <- pCPI +geom_segment(aes(x=2.7, y=111, xend=2.85, yend=111), size=0.7) +
  geom_segment(aes(x=2.7, y=107, xend=2.7, yend=111), size=0.7) +
  geom_segment(aes(x=2.85, y=107, xend=2.85, yend=111), size=0.7) +
  geom_text(x=2.8, y=116, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(EPZ_stat1$p.value, 
                                digits=2, scientific = FALSE))) +
  geom_segment(aes(x=2.7, y=125, xend=3, yend=125), size=0.7) +
  geom_segment(aes(x=2.7, y=121, xend=2.7, yend=125), size=0.7) +
  geom_segment(aes(x=3, y=121, xend=3, yend=125), size=0.7) +
  geom_text(x=2.85, y=132, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(EPZ_stat2$p.value, 
                                digits=2, scientific = FALSE)))+
  geom_segment(aes(x=2.7, y=139, xend=3.15, yend=139), size=0.7) +
  geom_segment(aes(x=2.7, y=135, xend=2.7, yend=139), size=0.7) +
  geom_segment(aes(x=3.15, y=135, xend=3.15, yend=139), size=0.7) +
  geom_text(x=2.9, y=146, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(EPZ_stat3$p.value, 
                                digits=2, scientific = FALSE))) +
  geom_segment(aes(x=2.7, y=152, xend=3.3, yend=152), size=0.7) +
  geom_segment(aes(x=2.7, y=148, xend=2.7, yend=152), size=0.7) +
  geom_segment(aes(x=3.3, y=148, xend=3.3, yend=152), size=0.7) +
  geom_text(x=2.9, y=159, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(EPZ_stat4$p.value, 
                                digits=2, scientific = FALSE)))
#  GSK126
pGSK <- pEPZ +geom_segment(aes(x=3.7, y=116, xend=3.85, yend=116), size=0.7) +
  geom_segment(aes(x=3.7, y=112, xend=3.7, yend=116), size=0.7) +
  geom_segment(aes(x=3.85, y=112, xend=3.85, yend=116), size=0.7) +
  geom_text(x=3.85, y=121, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(GSK_stat1$p.value, 
                                digits=2, scientific = FALSE))) +
  geom_segment(aes(x=3.7, y=127, xend=4, yend=127), size=0.7) +
  geom_segment(aes(x=3.7, y=123, xend=3.7, yend=127), size=0.7) +
  geom_segment(aes(x=4, y=123, xend=4, yend=127), size=0.7) +
  geom_text(x=3.85, y=133, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(GSK_stat2$p.value, 
                                digits=2, scientific = FALSE)))+
  geom_segment(aes(x=3.7, y=139, xend=4.15, yend=139), size=0.7) +
  geom_segment(aes(x=3.7, y=135, xend=3.7, yend=139), size=0.7) +
  geom_segment(aes(x=4.15, y=135, xend=4.15, yend=139), size=0.7) +
  geom_text(x=3.9, y=146, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(GSK_stat3$p.value, 
                                digits=2, scientific = FALSE))) +
  geom_segment(aes(x=3.7, y=152, xend=4.3, yend=152), size=0.7) +
  geom_segment(aes(x=3.7, y=148, xend=3.7, yend=152), size=0.7) +
  geom_segment(aes(x=4.3, y=148, xend=4.3, yend=152), size=0.7) +
  geom_text(x=3.9, y=159, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(GSK_stat4$p.value, 
                                digits=2, scientific = FALSE)))

#  MS1943
pMS <- pGSK +geom_segment(aes(x=4.7, y=111, xend=4.85, yend=111), size=0.7) +
  geom_segment(aes(x=4.7, y=107, xend=4.7, yend=111), size=0.7) +
  geom_segment(aes(x=4.85, y=107, xend=4.85, yend=111), size=0.7) +
  geom_text(x=4.8, y=116, size=4.2,  angle = 0,
            colour="gray40", family='sans',fontface="plain",
            label=paste0(format(MS_stat1$p.value, 
                                digits=2, scientific = FALSE))) +
  geom_segment(aes(x=4.7, y=125, xend=5, yend=125), size=0.7) +
  geom_segment(aes(x=4.7, y=121, xend=4.7, yend=125), size=0.7) +
  geom_segment(aes(x=5, y=121, xend=5, yend=125), size=0.7) +
  geom_text(x=4.9, y=132, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(MS_stat2$p.value, 
                                digits=2, scientific = TRUE)))+
  geom_segment(aes(x=4.7, y=139, xend=5.15, yend=139), size=0.7) +
  geom_segment(aes(x=4.7, y=135, xend=4.7, yend=139), size=0.7) +
  geom_segment(aes(x=5.15, y=135, xend=5.15, yend=139), size=0.7) +
  geom_text(x=4.9, y=146, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(MS_stat3$p.value, 
                                digits=2, scientific = FALSE))) +
  geom_segment(aes(x=4.7, y=152, xend=5.3, yend=152), size=0.7) +
  geom_segment(aes(x=4.7, y=148, xend=4.7, yend=152), size=0.7) +
  geom_segment(aes(x=5.3, y=148, xend=5.3, yend=152), size=0.7) +
  geom_text(x=4.9, y=159, angle = 0,size =4.2, 
            colour="gray40", family="Arial",fontface="plain",
            label=paste0(format(MS_stat4$p.value, 
                                digits=2, scientific = TRUE)))
pMS


filename_apotose=paste0(Dir_project,"figure/Fig2_a.svg")
svg(filename=filename_apotose,width=14,height=9,pointsize=12)
pMS
dev.off()

#### Fig 2f
apopose <- read.csv(paste0(Dir_project,"data/data_fig2f.csv"),header=TRUE)
apopose$day <- as.factor(apopose$day)

## stat D1
D1 <- apopose[which(apopose$day == 1),]
D1_stat1 <- t.test(D1[which(D1$drug %in% "DMSO"), "value"], 
                   D1[which(D1$drug %in% "MS1943 (4 uM)"), "value"],
                   alternative = "two.sided",paired = FALSE) 

D1_stat2 <- t.test(D1[which(D1$drug %in% "DMSO"), "value"], 
                   D1[which(D1$drug %in% "C24 (4 uM)"), "value"],
                   alternative = "two.sided",paired = FALSE)

D1_stat3 <- t.test(D1[which(D1$drug %in% "MS1943 (4 uM)"), "value"], 
                   D1[which(D1$drug %in% "C24 (4 uM)"), "value"],
                   alternative = "two.sided",paired = FALSE)

## stat D2
D2 <- apopose[which(apopose$day == 2),]
D2_stat1 <- t.test(D2[which(D2$drug %in% "DMSO"), "value"], 
                   D2[which(D2$drug %in% "MS1943 (4 uM)"), "value"],
                   alternative = "two.sided",paired = FALSE) 

D2_stat2 <- t.test(D2[which(D2$drug %in% "DMSO"), "value"], 
                   D2[which(D2$drug %in% "C24 (4 uM)"), "value"],
                   alternative = "two.sided",paired = FALSE)

D2_stat3 <- t.test(D2[which(D2$drug %in% "MS1943 (4 uM)"), "value"], 
                   D2[which(D2$drug %in% "C24 (4 uM)"), "value"],
                   alternative = "two.sided",paired = FALSE)

## stat D3
D3 <- apopose[which(apopose$day == 3),]
D3_stat1 <- t.test(D3[which(D3$drug %in% "DMSO"), "value"], 
                   D3[which(D3$drug %in% "MS1943 (4 uM)"), "value"],
                   alternative = "two.sided",paired = FALSE) 

D3_stat2 <- t.test(D3[which(D3$drug %in% "DMSO"), "value"], 
                   D3[which(D3$drug %in% "C24 (4 uM)"), "value"],
                   alternative = "two.sided",paired = FALSE)

D3_stat3 <- t.test(D3[which(D3$drug %in% "MS1943 (4 uM)"), "value"], 
                   D3[which(D3$drug %in% "C24 (4 uM)"), "value"],
                   alternative = "two.sided",paired = FALSE)

## stat D4
D4 <- apopose[which(apopose$day == 4),]
D4_stat1 <- t.test(D4[which(D4$drug %in% "DMSO"), "value"], 
                   D4[which(D4$drug %in% "MS1943 (4 uM)"), "value"],
                   alternative = "two.sided",paired = FALSE) 

D4_stat2 <- t.test(D4[which(D4$drug %in% "DMSO"), "value"], 
                   D4[which(D4$drug %in% "C24 (4 uM)"), "value"],
                   alternative = "two.sided",paired = FALSE)

D4_stat3 <- t.test(D4[which(D4$drug %in% "MS1943 (4 uM)"), "value"], 
                   D4[which(D4$drug %in% "C24 (4 uM)"), "value"],
                   alternative = "two.sided",paired = FALSE)

p <- ggplot(apopose, aes(x=day, y=value)) + 
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.8), dotsize = 0.7,
               aes(fill=drug)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
               geom="pointrange", color="black",
               mapping = aes(group = drug),
               position=position_dodge(0.8)) +
  ylim(0, 7) +
  scale_fill_manual(values=c("#1b9e77","#d95f02", "#7570b3"))+
  theme_linedraw(base_size=12) +
  theme(axis.text.x =black.italic.12.text) +
  xlab("Days after treatment") + ylab("Normalized apoptotic counts") +labs(fill = "")
p

pD1 <- p + geom_segment(aes(x=0.75, y=1.8, xend=1, yend=1.8), size=0.7) +
  geom_segment(aes(x=0.75, y=1.6, xend=0.75, yend=1.8), size=0.7) +
  geom_segment(aes(x=1, y=1.6, xend=1, yend=1.8), size=0.7) +
  geom_text(x=0.9, y=1.95, size=4.2,  angle = 0,
            colour="gray40", family='sans',fontface="plain",
            label=paste0(format(D1_stat2$p.value, 
                                digits=2, scientific = FALSE))) +
  geom_segment(aes(x=1, y=3, xend=1.25, yend=3), size=0.7) +
  geom_segment(aes(x=1.25, y=2.8, xend=1.25, yend=3), size=0.7) +
  geom_segment(aes(x=1, y=2.8, xend=1, yend=3), size=0.7) +
  geom_text(x=1.15, y=3.2, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(D1_stat1$p.value, 
                                digits=2, scientific = FALSE))) +
  geom_segment(aes(x=0.75, y=3.6, xend=1.25, yend=3.6), size=0.7) +
  geom_segment(aes(x=1.25, y=3.4, xend=1.25, yend=3.6), size=0.7) +
  geom_segment(aes(x=0.75, y=3.4, xend=0.75, yend=3.6), size=0.7) +
  geom_text(x=1, y=3.8, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(D1_stat3$p.value, 
                                digits=2, scientific = FALSE)))
pD1

pD2 <- pD1 + geom_segment(aes(x=1.75, y=2, xend=2, yend=2), size=0.7) +
  geom_segment(aes(x=1.75, y=1.8, xend=1.75, yend=2), size=0.7) +
  geom_segment(aes(x=2, y=1.8, xend=2, yend=2), size=0.7) +
  geom_text(x=1.9, y=2.15, size=4.2,  angle = 0,
            colour="gray40", family='sans',fontface="plain",
            label=paste0(format(D2_stat2$p.value, 
                                digits=2, scientific = FALSE))) +
  geom_segment(aes(x=2, y=3.3, xend=2.25, yend=3.3), size=0.7) +
  geom_segment(aes(x=2.25, y=3.1, xend=2.25, yend=3.3), size=0.7) +
  geom_segment(aes(x=2, y=3.1, xend=2, yend=3.3), size=0.7) +
  geom_text(x=2.15, y=3.5, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(D2_stat1$p.value, 
                                digits=2, scientific = FALSE))) +
  geom_segment(aes(x=1.75, y=3.9, xend=2.25, yend=3.9), size=0.7) +
  geom_segment(aes(x=2.25, y=3.7, xend=2.25, yend=3.9), size=0.7) +
  geom_segment(aes(x=1.75, y=3.7, xend=1.75, yend=3.9), size=0.7) +
  geom_text(x=2, y=4.1, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(D2_stat3$p.value, 
                                digits=2, scientific = FALSE)))
pD2  

pD3 <- pD2 + geom_segment(aes(x=2.75, y=2.7, xend=3, yend=2.7), size=0.7) +
  geom_segment(aes(x=2.75, y=2.5, xend=2.75, yend=2.7), size=0.7) +
  geom_segment(aes(x=3, y=2.5, xend=3, yend=2.7), size=0.7) +
  geom_text(x=2.9, y=2.9, size=4.2,  angle = 0,
            colour="gray40", family='sans',fontface="plain",
            label=paste0(format(D3_stat2$p.value, 
                                digits=2, scientific = FALSE))) +
  geom_segment(aes(x=3, y=4.6, xend=3.25, yend=4.6), size=0.7) +
  geom_segment(aes(x=3.25, y=4.4, xend=3.25, yend=4.6), size=0.7) +
  geom_segment(aes(x=3, y=4.4, xend=3, yend=4.6), size=0.7) +
  geom_text(x=3.15, y=4.8, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(D3_stat1$p.value, 
                                digits=2, scientific = FALSE))) +
  geom_segment(aes(x=2.75, y=5.2, xend=3.25, yend=5.2), size=0.7) +
  geom_segment(aes(x=3.25, y=5, xend=3.25, yend=5.2), size=0.7) +
  geom_segment(aes(x=2.75, y=5, xend=2.75, yend=5.2), size=0.7) +
  geom_text(x=3, y=5.4, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(D3_stat3$p.value, 
                                digits=2, scientific = FALSE)))
pD3

pD4 <- pD3 + geom_segment(aes(x=3.75, y=3.1, xend=4, yend=3.1), size=0.7) +
  geom_segment(aes(x=3.75, y=2.9, xend=3.75, yend=3.1), size=0.7) +
  geom_segment(aes(x=4, y=2.9, xend=4, yend=3.1), size=0.7) +
  geom_text(x=3.9, y=3.3, size=4.2,  angle = 0,
            colour="gray40", family='sans',fontface="plain",
            label=paste0(format(D4_stat2$p.value, 
                                digits=2, scientific = FALSE))) +
  geom_segment(aes(x=4, y=5.8, xend=4.25, yend=5.8), size=0.7) +
  geom_segment(aes(x=4.25, y=5.6, xend=4.25, yend=5.8), size=0.7) +
  geom_segment(aes(x=4, y=5.6, xend=4, yend=5.8), size=0.7) +
  geom_text(x=4.15, y=6, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(D4_stat1$p.value, 
                                digits=2, scientific = FALSE))) +
  geom_segment(aes(x=3.75, y=6.4, xend=4.25, yend=6.4), size=0.7) +
  geom_segment(aes(x=4.25, y=6.2, xend=4.25, yend=6.4), size=0.7) +
  geom_segment(aes(x=3.75, y=6.2, xend=3.75, yend=6.4), size=0.7) +
  geom_text(x=4, y=6.6, size=4.2,  angle = 0,
            colour="gray40", family='Arial',fontface="plain",
            label=paste0(format(D4_stat3$p.value, 
                                digits=2, scientific = FALSE)))
pD4

filename_apotose=paste0(Dir_project,"figure/Fig2_f.svg")
svg(filename=filename_apotose,width=7,height=7,pointsize=12)
pD4
dev.off()

#### Fig 5b
IHC_data <- read.csv(paste0(Dir_project,"data/IHC_data.csv"),header=TRUE)
IHC_data$type <- as.factor(IHC_data$type)

## Ezh2
IHC_Ezh2 <- IHC_data[which(IHC_data$cell %in% "Ezh2"),]
IHC_Ezh2_stat <- t.test(IHC_Ezh2[which(IHC_Ezh2$type %in% "Vehicle"), 
                                 "values"], 
                        IHC_Ezh2[which(IHC_Ezh2$type %in% "MS1943"), 
                                 "values"],
                        alternative = "two.sided",paired = FALSE) 

p <- ggplot(IHC_Ezh2,
            aes(x=type, y=values)) + 
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.8), dotsize = 0.8,
               aes(fill=type)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
               geom="pointrange", color="black", size = 0.5,
               mapping = aes(group = type),
               position=position_dodge(0.8)) +
  ylim(0, 110) +
  scale_fill_manual(values=c("#2c7bb6","#d7191c"))+
  theme_linedraw(base_size=12) +
  theme(axis.text.x = element_text(color = "black", size = 12,
                                   angle = 45, hjust = 1)) +
  xlab("") + ylab("EZH2 positive cells (%)") +labs(fill = "") +
  geom_segment(aes(x=1, y=103, xend=2, yend=103), size=0.7) +
  geom_segment(aes(x=1, y=99, xend=1, yend=103), size=0.7) +
  geom_segment(aes(x=2, y=99, xend=2, yend=103), size=0.7) +
  geom_text(x=1.5, y=107, size=4.0,  
            colour="gray40", family='Arial',fontface="plain",
            label=paste0("p=",format(IHC_Ezh2_stat$p.value, 
                                     digits=2, scientific = TRUE)))
p

filename_apotose=paste0(Dir_project,"figure/Fig5_b_EZH2.svg")
svg(filename=filename_apotose,width=4,height=4,pointsize=12)
p
dev.off()

## H3K27me3
IHC_H3K27me3 <- IHC_data[which(IHC_data$cell %in% "H3K27me3"),]
IHC_H3K27me3_stat <- t.test(IHC_H3K27me3[which(IHC_H3K27me3$type %in% "Vehicle"), 
                                         "values"], 
                            IHC_H3K27me3[which(IHC_H3K27me3$type %in% "MS1943"), 
                                         "values"],
                            alternative = "two.sided",paired = FALSE) 

p <- ggplot(IHC_H3K27me3,
            aes(x=type, y=values)) + 
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.8), dotsize = 0.8,
               aes(fill=type)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
               geom="pointrange", color="black", size = 0.5,
               mapping = aes(group = type),
               position=position_dodge(0.8)) +
  ylim(0, 110) +
  scale_fill_manual(values=c("#2c7bb6","#d7191c"))+
  theme_linedraw(base_size=12) +
  theme(axis.text.x = element_text(color = "black", size = 12,
                                   angle = 45, hjust = 1)) +
  xlab("") + ylab("H3K27me3 positive cells (%)") +labs(fill = "") +
  geom_segment(aes(x=1, y=103, xend=2, yend=103), size=0.7) +
  geom_segment(aes(x=1, y=99, xend=1, yend=103), size=0.7) +
  geom_segment(aes(x=2, y=99, xend=2, yend=103), size=0.7) +
  geom_text(x=1.5, y=107, size=4.0,  
            colour="gray40", family='Arial',fontface="plain",
            label=paste0("p=",format(IHC_H3K27me3_stat$p.value, 
                                     digits=2, scientific = FALSE)))
p

filename_apotose=paste0(Dir_project,"figure/Fig5_b_H3K27me3.svg")
svg(filename=filename_apotose,width=4,height=4,pointsize=12)
p
dev.off()

## Cleaved Caspase-3
IHC_Casp3 <- IHC_data[which(IHC_data$cell %in% "Cleaved caspase-3"),]
IHC_Casp3_stat <- t.test(IHC_Casp3[which(IHC_Casp3$type %in% "Vehicle"), 
                                   "values"], 
                         IHC_Casp3[which(IHC_Casp3$type %in% "MS1943"), 
                                   "values"],
                         alternative = "two.sided",paired = FALSE) 

p <- ggplot(IHC_Casp3,
            aes(x=type, y=values)) + 
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.8), dotsize = 0.8,
               aes(fill=type)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
               geom="pointrange", color="black", size = 0.5,
               mapping = aes(group = type),
               position=position_dodge(0.8)) +
  ylim(0, 20) +
  scale_fill_manual(values=c("#2c7bb6","#d7191c"))+
  theme_linedraw(base_size=12) +
  theme(axis.text.x = element_text(color = "black", size = 12,
                                   angle = 45, hjust = 1)) +
  xlab("") + ylab("Cleaved Caspase-3 positive cells (%)") +labs(fill = "") +
  geom_segment(aes(x=1, y=15, xend=2, yend=15), size=0.7) +
  geom_segment(aes(x=1, y=13, xend=1, yend=15), size=0.7) +
  geom_segment(aes(x=2, y=13, xend=2, yend=15), size=0.7) +
  geom_text(x=1.5, y=17, size=4.0,  
            colour="gray40", family='Arial',fontface="plain",
            label=paste0("p=",format(IHC_Casp3_stat$p.value, 
                                     digits=2, scientific = FALSE)))
p

filename_apotose=paste0(Dir_project,"figure/Fig5_b_Cas3.svg")
svg(filename=filename_apotose,width=4,height=4,pointsize=12)
p
dev.off()

## Ki-67
IHC_Ki67 <- IHC_data[which(IHC_data$cell %in% "Ki-67"),]
IHC_Ki67_stat <- t.test(IHC_Ki67[which(IHC_Ki67$type %in% "Vehicle"), 
                                 "values"], 
                        IHC_Ki67[which(IHC_Ki67$type %in% "MS1943"), 
                                 "values"],
                        alternative = "two.sided",paired = FALSE) 

p <- ggplot(IHC_Ki67,
            aes(x=type, y=values)) + 
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.8), dotsize = 0.8,
               aes(fill=type)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
               geom="pointrange", color="black", size = 0.5,
               mapping = aes(group = type),
               position=position_dodge(0.8)) +
  ylim(0, 110) +
  scale_fill_manual(values=c("#2c7bb6","#d7191c"))+
  theme_linedraw(base_size=12) +
  theme(axis.text.x = element_text(color = "black", size = 12,
                                   angle = 45, hjust = 1)) +
  xlab("") + ylab("Ki-67 positive cells (%)") +labs(fill = "") +
  geom_segment(aes(x=1, y=103, xend=2, yend=103), size=0.7) +
  geom_segment(aes(x=1, y=99, xend=1, yend=103), size=0.7) +
  geom_segment(aes(x=2, y=99, xend=2, yend=103), size=0.7) +
  geom_text(x=1.5, y=107, size=4.0,  
            colour="gray40", family='Arial',fontface="plain",
            label=paste0("p=",format(IHC_Ki67_stat$p.value, 
                                     digits=2, scientific = FALSE)))
p

filename_apotose=paste0(Dir_project,"figure/Fig5_b_Ki67.svg")
svg(filename=filename_apotose,width=4,height=4,pointsize=12)
p
dev.off()


#### RNAseq data
dir <- "results/RNAseq"
samples <- read.table(file.path(dir,"meta_EZH2.csv"), header=TRUE,sep=",")

xbp_data <- read.table(file.path(dir,"salmon/XBP_transcript.csv"),
                       header = TRUE, sep=",")

xbp_data_all <- merge(xbp_data,samples,by="run")
xbp_data_all$run <- as.factor(xbp_data_all$run)
xbp_data_all$treated <- as.character(xbp_data_all$treated)
xbp_data_all$treated[which(xbp_data_all$treated %in% "degrader")] <- "MS1943"
xbp_data_all$treated <- as.factor(xbp_data_all$treated)
xbp_data_all$transcript <- as.factor(xbp_data_all$transcript)
dim(xbp_data_all)
dtu <- xbp_data_all

dtu_summary <- data_summary(xbp_data_all, varname="TPM", 
                    groupnames=c("name", "treated"))
# Convert treated to a factor variable
dtu_summary$treated=as.factor(dtu$treated)
head(dtu_summary)
dim(dtu_summary)
#

#
p <- ggplot(dtu, aes(x=name, y=TPM)) + 
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.5), dotsize = 0.7,
               aes(fill=treated)) +
  ylim(0, 125) +
  scale_fill_manual(values=c("#d7191c","#2c7bb6"))+
  theme_linedraw(base_size=12) +
  theme(axis.text.x =black.italic.12.text) +
  xlab("") + ylab("Transcripts per Million") +labs(fill = "")
p

filename_DTU=paste0(Dir_project,"Fig6_c.svg")
svg(filename=filename_DTU,width=9,height=7,pointsize=12)
p
dev.off()


#### qPCR results
#### Figure 6 e
ms1943 <- read.csv(paste0(Dir_project,"Fig6_qRT.PCR_XY1943.csv"),header=TRUE)

p <- ggplot(ms1943, aes(x=gene, y=fold.change)) + 
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.5), dotsize = 0.7,
               aes(fill=cell_line)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
               geom="pointrange", color="black",
               mapping = aes(group = cell_line),
               position=position_dodge(0.5)) +
  ylim(0, 3.25) +
  scale_fill_manual(values=c("#d01c8b","#b8e186", "#f1b6da",  "#4dac26"))+
  theme_linedraw(base_size=12) +
  theme(axis.text.x =black.italic.12.text) +
  xlab("") + ylab("Fold change (MS1943/DMSO)") +labs(fill = "")
p

filename_ms1943=paste0(Dir_project,"Fig6_e.svg")
svg(filename=filename_ms1943,width=7,height=7,pointsize=12)
p
dev.off()

#### Figure 6 f
c24 <- read.csv(paste0(Dir_project,"Fig6_qRT.PCR_c24.csv"),header=TRUE)

p <- ggplot(c24, aes(x=gene, y=fold.change)) + 
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.7), dotsize = 0.7,
               aes(fill=cell_line)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
               geom="pointrange", color="black",
               mapping = aes(group = cell_line),
               position=position_dodge(0.7)) +
  ylim(0, 3.25) +
  scale_fill_manual(values=c("#d01c8b","#b8e186", "#f1b6da",  "#4dac26"))+
  theme_linedraw(base_size=12) +
  theme(axis.text.x =black.italic.12.text) +
  xlab("") + ylab("Fold change (C24/DMSO)") +labs(fill = "")
p

filename_c24=paste0(Dir_project,"Fig6_f.svg")
svg(filename=filename_c24,width=7,height=7,pointsize=12)
p
dev.off()

#### Fig 11
gene_time <- read.csv(paste0(Dir_project,"data/data_fig11.csv"),header=TRUE)
gene_time$time <- as.factor(gene_time$time)

p <- ggplot(gene_time, aes(x=gene, y=value)) + 
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.8), dotsize = 0.7,
               aes(fill=time)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
               geom="pointrange", color="black",
               mapping = aes(group = time),
               position=position_dodge(0.8)) +
  ylim(0, 16) +
  scale_fill_manual(labels = c("4 hrs", "24 hrs", "48 hrs"),
                    values=c("#ffffcc","#41b6c4", "#253494"),
                    name ="") +
  theme_linedraw(base_size=12) +
  theme(axis.text.x =black.italic.12.text) +
  xlab("") + ylab("Fold change (MS1943 vs DMSO)") 
p

filename_apotose=paste0(Dir_project,"figure/Fig11.svg")
svg(filename=filename_apotose,width=7,height=7,pointsize=12)
p
dev.off()



