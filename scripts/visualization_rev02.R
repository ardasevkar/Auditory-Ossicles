# Setup -------------------------------------------------------------------
Sys.setlocale(category = 'LC_ALL','en_US.UTF-8')
options(timeout=1000, scipen = 999)
`%nin%` = Negate(`%in%`)
sign <- function(sayı){
  
  x <- ifelse(sayı <= 0.05 & sayı > 0.01, "*",
              ifelse(sayı <= 0.01 & sayı > 0.001 , "**",
                     ifelse(sayı <= 0.001, "*** (p ≤ 0.001)", "n.s.")))
  return(x)}


setwd("../")

# Loading Libraries -------------------------------------------------------
library(raster)
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)
library(ggpubr)
library(ggforce)
library(rnaturalearth)
library(rnaturalearthdata)
library(superb)
library(readxl)
library(reshape2)
library(tidyverse)
library(png)
library(grid)
library(ggh4x)
library(DescTools)
library(tidyverse)
#palette <- c("#e63946", "#130870","#2b35af", "#4360ed", "#4895ef", "#faf0ca")
#palette <- c("#264653", "#2a9d8f","#e9c46a", "#f4a261", "#e76f51")


palette <- c("#e63845", "#07213f","#048BA8", "#FAF3DD", "#E0913E",  "#F194B4")
#https://coolors.co/e63845-07213f-048ba8-faf3dd-e0913e


col_p <- palette[1]
col_s <- palette[2]
col_m <- palette[3]
col_i <- palette[4]
col_ui <- palette[5]
col_t <- palette[6]


# Loading data ------------------------------------------------------------
all <- read.table("data/seqstats.txt", header = T, sep = "\t")
all$Bone <- factor(all$Bone, levels = c("Petrous", "Stapes", "Malleus", "Incus", "Unidentified Ossicle", "Tooth"))
all$bone_other <- ifelse(all$Bone %nin% "Petrous", "Non-Petrous", "Petrous")
all$bone_other <- factor(all$bone_other, levels = c("Petrous", "Non-Petrous"))
all$temp <- paste0(all$comb, all$Sample)

seqstats <- read.table("data/seqstats.txt", header = T, sep = "\t")
seqstats$Bone[seqstats$Bone %in% "Unidentified Ossicle"] <- "UI Ossicle"
seqstats$Bone <- factor(seqstats$Bone, levels = c("Petrous", "Stapes", "Malleus", "Incus", "UI Ossicle", "Tooth"))
seqstats$comb[seqstats$comb %in% "Petrous|Unidentified Ossicle"] <- "Petrous|UI Ossicle"


# Figure 2 --------------------------------------
# Code below contains scripts for visualization of Main Figure 2
fig2 <- all %>%
  filter(Bone %nin% "Petrous") %>%
  filter(variable %in% c("Human Proportion (%)", "Clonality (%)", 
                         "Avg Read Length (bp)", "PMD (%)", "Cont. Est."))

fig2$variable <- as.character(fig2$variable)
fig2$variable[fig2$variable %in% "Human Proportion (%)"] <- "Human Proportion"
fig2$variable[fig2$variable %in% "Clonality (%)"] <- "Clonality"
fig2$variable[fig2$variable %in% "Avg Read Length (bp)"] <- "Mean Fragment Length"
fig2$variable[fig2$variable %in% "Cont. Est."] <- "Contamination"
fig2$variable[fig2$variable %in% "PMD (%)"] <- "PMD"

fig2$variable <- factor(fig2$variable, levels = c("Human Proportion", "Mean Fragment Length", 
                                                    "PMD", "Clonality", "Contamination"))

###### ____First Panel -----------
p_1a <- fig2 %>%
  filter(variable %in% "Human Proportion") %>%
  ggplot(aes(x = Bone, y = diff))+
  geom_hline(yintercept = 0, color = "#ff2800", linetype = "dashed", size = 1)+
  geom_boxplot(aes(fill = Bone, color = Bone), outliers = FALSE, width = .5)+
  labs(y = "Change Relative to Petrous (%)",
       x = "")+
  geom_quasirandom(aes(color = Bone, fill = Bone), shape = 21, size = 3, stroke = .5, alpha = .8,
                   show.legend = F, width = 0.3) +
  facet_grid(.~variable, scales = "free")+
  #facetted_pos_scales(y = list(scale_y_continuous(breaks = c(0, 100, 200, 300, 500, 750, 1500))))+
  scale_fill_manual(values = c(col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("grey70","black", "black", "black", "black"))+
  scale_x_discrete(labels = c("Stapes", "Malleus", "Incus","UI Ossicle", "Tooth"))+
  scale_y_continuous(breaks = c(0, 100, 200, 300, 500, 750, 1500)) +
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p_1a

# Distribution of Human Proportion  Scatter Plot 
hp <- seqstats %>%
  filter(variable %in% "Human Proportion (%)") %>%
  as.data.frame()

hp$Bone <- factor(hp$Bone, levels = c("Petrous", "Stapes", "Malleus", "Incus", "UI Ossicle"))

hp$comb[hp$comb %in% "Petrous|Stapes"] <- 1
hp$comb[hp$comb %in% "Petrous|Malleus"] <- 2
hp$comb[hp$comb %in% "Petrous|Incus"] <- 3
hp$comb[hp$comb %in% "Petrous|UI Ossicle"] <- 4

for(i in 1:4){
  temp <- hp %>%
    ungroup() %>%
    filter(comb %in% i) %>%
    select(Sample, Bone, value) %>%
    spread(Bone, value)
  
  test <- wilcox.test(temp[,2], temp[,3], paired = T)$p.value
  star <- paste0(sign(test), " (p=", round(test, digits = 3), ")");star
  assign(paste0("star_", i), star)
}

star_1;star_2;star_3;star_4

p_1b <- hp %>%
  ggplot(aes(x = Bone, y = value))+
  ylim(0,90)+
  labs(x = "", y = "Human Proportion (%)")+
  geom_line(aes(group = Sample), color = "#121212", size = .75, show.legend = F)+
  geom_point(aes(fill = Bone, color = Bone), shape = 21, size = 5, show.legend = T, stroke = .75)+
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui, col_t))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  showSignificance(c(1,2), 85, -0.0001, star_1,
                   textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
                   panel = list(comb = 1))+
  showSignificance(c(1,2), 85, -0.0001, star_3,
                   textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
                   panel = list(comb = 3))+
  facet_grid(.~comb, space = "free", scales = "free_x")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_blank(),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p_1b


p_A <- ggarrange(p_1a, p_1b, ncol = 2, align = "v", widths = c(.3, 1));p_A



###### ____Second Panel ------------

p_2a <- fig2 %>%
  filter(variable %in% "Mean Fragment Length") %>%
  ggplot(aes(x = Bone, y = diff))+
  geom_hline(yintercept = 0, color = "#ff2800", linetype = "dashed", size = 1)+
  geom_boxplot(aes(fill = Bone, color = Bone), outliers = FALSE, width = .5)+
  labs(y = "Change Relative to Petrous (%)",
       x = "")+
  geom_quasirandom(aes(color = Bone, fill = Bone), shape = 21, size = 3, stroke = .5, alpha = .8,
                   show.legend = F, width = 0.3) +
  facet_grid(.~variable, scales = "free")+
  scale_fill_manual(values = c(col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("grey70","black", "black", "black", "black"))+
  scale_x_discrete(labels = c("Stapes", "Malleus", "Incus","UI Ossicle", "Tooth"))+
  scale_y_continuous(limits = c(-20,20))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p_2a



# Scatter Plot: Distribution of mean fragment length
avg <- seqstats %>%
  filter(variable %in% "Avg Read Length (bp)") %>%
  as.data.frame()

avg$Bone <- factor(avg$Bone, levels = c("Petrous", "Stapes", "Malleus", "Incus", "UI Ossicle"))

avg$comb[avg$comb %in% "Petrous|Stapes"] <- 1
avg$comb[avg$comb %in% "Petrous|Malleus"] <- 2
avg$comb[avg$comb %in% "Petrous|Incus"] <- 3
avg$comb[avg$comb %in% "Petrous|UI Ossicle"] <- 4


for(i in 1:4){
  temp <- avg %>%
    ungroup() %>%
    filter(comb %in% i) %>%
    select(Sample, Bone, value) %>%
    spread(Bone, value)
  
  test <- wilcox.test(temp[,2], temp[,3], paired = T)$p.value
  star <- paste0(sign(test), " (p=", round(test, digits = 3), ")");star
  assign(paste0("star_", i), star)
}

star_1;star_2;star_3;star_4


# Temporary script for the manuscript: Difference in fragment lenghts in Stapes vs Petrous
temp <- avg %>%
  filter(type %in% "Stapes")

mean(temp$value[temp$Bone %in% "Stapes"])-mean(temp$value[temp$Bone %in% "Petrous"])
wilcox.test(temp$value[temp$Bone %in% "Stapes"], temp$value[temp$Bone %in% "Petrous"],
            paired = T)$p.value

p_2b <- avg %>%
  ggplot(aes(x = Bone, y = value))+
  ylim(40,100)+
  labs(x = "", y = "Mean Fragment Length (bp)")+
  geom_line(aes(group = Sample), color = "#121212", size = .75, show.legend = F)+
  geom_point(aes(fill = Bone, color = Bone), shape = 21, size = 5, show.legend = T, stroke = .75)+
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  showSignificance(c(1,2), 95, -0.0001, star_1,
                   textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
                   panel = list(comb = 1))+
  showSignificance(c(1,2), 92, -0.0001, star_2,
                   textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
                   panel = list(comb = 2))+
  facet_grid(.~comb, space = "free", scales = "free_x")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_blank(),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p_2b


p_B <- ggarrange(p_2a, p_2b, ncol = 2, align = "v", widths = c(.3, 1));p_B



###### ____Third Panel ------------
p_3a <- fig2 %>%
  filter(variable %in% "PMD") %>%
  ggplot(aes(x = Bone, y = diff))+
  geom_hline(yintercept = 0, color = "#ff2800", linetype = "dashed", size = 1)+
  geom_boxplot(aes(fill = Bone, color = Bone), outliers = FALSE, width = .5)+
  labs(y = "Change Relative to Petrous (%)",
       x = "")+
  geom_quasirandom(aes(color = Bone, fill = Bone), shape = 21, size = 3, stroke = .5, alpha = .8,
                   show.legend = F, width = 0.3) +
  facet_grid(.~variable, scales = "free")+
  scale_fill_manual(values = c(col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("grey70","black", "black", "black", "black"))+
  scale_x_discrete(labels = c("Stapes", "Malleus", "Incus","UI Ossicle", "Tooth"))+
  scale_y_continuous(limits = c(-50,100))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p_3a


# Scatter Plot: Distribution of postmortem damage values
pmd <- seqstats %>%
  filter(variable %in% "PMD (%)") %>%
  as.data.frame()

pmd$Bone <- factor(pmd$Bone, levels = c("Petrous", "Stapes", "Malleus", "Incus", "UI Ossicle"))

pmd$comb[pmd$comb %in% "Petrous|Stapes"] <- 1
pmd$comb[pmd$comb %in% "Petrous|Malleus"] <- 2
pmd$comb[pmd$comb %in% "Petrous|Incus"] <- 3
pmd$comb[pmd$comb %in% "Petrous|UI Ossicle"] <- 4

for(i in 1:4){
  temp <- pmd %>%
    ungroup() %>%
    filter(comb %in% i) %>%
    select(Sample, Bone, value) %>%
    spread(Bone, value)
  
  test <- wilcox.test(temp[,2], temp[,3], paired = T)$p.value
  star <- paste0(sign(test), " (p=", round(test, digits = 3), ")");star
  assign(paste0("star_", i), star)
}
star_1;star_2;star_3;star_4

p_3b <- pmd %>%
  ggplot(aes(x = Bone, y = value))+
  ylim(0,60)+
  labs(x = "", y = "PMD (%)")+
  geom_line(aes(group = Sample), color = "#121212", size = .75, show.legend = F)+
  geom_point(aes(fill = Bone, color = Bone), shape = 21, size = 5, show.legend = T, stroke = .75)+
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui, col_t))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  facet_grid(.~comb, space = "free", scales = "free_x")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_blank(),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p_3b

p_C <- ggarrange(p_3a, p_3b, ncol = 2, align = "v", widths = c(.3, 1));p_C


###### ____Fourth Panel ------------
p_4a <- fig2 %>%
  filter(variable %in% "Clonality") %>%
  mutate(diff = ifelse(diff < -100, -200, diff)) %>%
  ggplot(aes(x = Bone, y = diff))+
  geom_hline(yintercept = 0, color = "#ff2800", linetype = "dashed", size = 1)+
  geom_boxplot(aes(fill = Bone, color = Bone), outliers = FALSE, width = .5)+
  labs(y = "Change Relative to Petrous (%)",
       x = "")+
  geom_quasirandom(aes(color = Bone, fill = Bone), shape = 21, size = 3, stroke = .5, alpha = .8,
                   show.legend = F, width = 0.3) +
  scale_fill_manual(values = c(col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("grey70","black", "black", "black", "black"))+
  scale_x_discrete(labels = c("Stapes", "Malleus", "Incus","UI Ossicle", "Tooth"))+
  scale_y_continuous(limits = c(-250,100), breaks = c(-250, 0, 50, 100), 
                     labels = c(-600, 0, 50, 100))+
  facet_grid(.~variable)+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank())+
  annotate("rect", xmin = 0.40, xmax = 0.60, ymin = -195, ymax = -140, 
           fill = "white", color = NA) +
  annotate("segment", x = 0.35, xend = 0.65, y = -170, yend = -150, 
           color = "black", linewidth = 1) +
  annotate("segment", x = 0.35, xend = 0.65, y = -185, yend = -165, 
           color = "black", linewidth = 1);p_4a
  




# Scatter Plot: Distribution of clonality values
clon <- seqstats %>%
  filter(variable %in% "Clonality (%)") %>%
  as.data.frame()

clon$Bone <- factor(clon$Bone, levels = c("Petrous", "Stapes", "Malleus", "Incus", "UI Ossicle"))

clon$comb[clon$comb %in% "Petrous|Stapes"] <- 1
clon$comb[clon$comb %in% "Petrous|Malleus"] <- 2
clon$comb[clon$comb %in% "Petrous|Incus"] <- 3
clon$comb[clon$comb %in% "Petrous|UI Ossicle"] <- 4

for(i in 1:4){
  temp <- clon %>%
    ungroup() %>%
    filter(comb %in% i) %>%
    select(Sample, Bone, value) %>%
    spread(Bone, value)
  
  test <- wilcox.test(temp[,2], temp[,3], paired = T)$p.value
  star <- paste0(sign(test), " (p=", round(test, digits = 3), ")");star
  assign(paste0("star_", i), star)
}

star_1;star_2;star_3;star_4


p_4b <- clon %>%
  ggplot(aes(x = Bone, y = value))+
  labs(x = "", y = "Clonality (%)")+
  geom_line(aes(group = Sample), color = "#121212", size = .75, show.legend = F)+
  geom_point(aes(fill = Bone, color = Bone), shape = 21, size = 5, show.legend = T, stroke = .75)+
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
facet_grid(.~comb, space = "free", scales = "free_x")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_blank(),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p_4b

p_D <- ggarrange(p_4a, p_4b, ncol = 2, align = "v", widths = c(.3, 1));p_D



###### ____Fifth Panel ------------
p_5a <- fig2 %>%
  filter(variable %in% "Contamination") %>%
  ggplot(aes(x = Bone, y = -diff))+
  geom_hline(yintercept = 0, color = "#ff2800", linetype = "dashed", size = 1)+
  geom_boxplot(aes(fill = Bone, color = Bone), outliers = FALSE, width = .5)+
  labs(y = "Change Relative to Petrous (%)",
       x = "")+
  geom_quasirandom(aes(color = Bone, fill = Bone), shape = 21, size = 3, stroke = .5, alpha = .8,
                   show.legend = F, width = 0.3) +
  facet_grid(.~variable, scales = "free")+
  scale_fill_manual(values = c(col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("grey70","black", "black", "black", "black"))+
  scale_x_discrete(labels = c("Stapes", "Malleus", "Incus","UI Ossicle", "Tooth"))+
  #scale_y_continuous(limits = c(-600,100), breaks = c(-600, 0, 25, 50, 100))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p_5a


# Scatter Plot: Distribution of contamination estimates
cont <- seqstats %>%
  filter(variable %in% "Cont. Est.") %>%
  as.data.frame()

cont$Bone <- factor(cont$Bone, levels = c("Petrous", "Stapes", "Malleus", "Incus", "UI Ossicle"))

cont$comb[cont$comb %in% "Petrous|Stapes"] <- 1
cont$comb[cont$comb %in% "Petrous|Malleus"] <- 2
cont$comb[cont$comb %in% "Petrous|Incus"] <- 3
cont$comb[cont$comb %in% "Petrous|UI Ossicle"] <- 4

for(i in 1:4){
  temp <- cont %>%
    ungroup() %>%
    filter(comb %in% i) %>%
    select(Sample, Bone, value) %>%
    spread(Bone, value)
  
  test <- wilcox.test(temp[,2], temp[,3], paired = T)$p.value
  star <- paste0(sign(test), " (p=", round(test, digits = 3), ")");star
  assign(paste0("star_", i), star)
}

star_1;star_2;star_3;star_4


p_5b <- cont %>%
  ggplot(aes(x = Bone, y = value))+
  # ylim(0,23)+
  labs(x = "", y = "Contamination Estimate")+
  geom_line(aes(group = Sample), color = "#121212", size = .75, show.legend = F)+
  geom_point(aes(fill = Bone, color = Bone), shape = 21, size = 5, show.legend = T, stroke = .75)+
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  facet_grid(.~comb, space = "free", scales = "free_x")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_blank(),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p_5b

p_E <- ggarrange(p_5a, p_5b, ncol = 2, align = "v", widths = c(.3, 1));p_E



# Combining all panels for Main Figure 2
p2 <- ggarrange(p_A, p_B, p_C, p_D, p_E, ncol = 1, align = "h", labels = LETTERS[1:5]);p2
ggsave("figures/png/figure2.png", dpi = 600, device = png, height = 17, width = 15)
ggsave("figures/svg/figure2.svg", dpi = 600, device = svg, height = 17, width = 15)



# Figure 3 --------------------------------------

##### ____ First panel ------
# Loading dataset
data <- read_xlsx("data/Supplemental_Tables_rev02.xlsx", sheet = "S1", skip = 1)

# Remove tooth
data <- data %>%
  filter(Bone %nin% "Tooth") %>%
  filter(`Library Type` %nin% "petrous&tooth")

table(data$Bone)

x <- match(x = c("Mean Read Length (BAM, bp)", 
                 "PMD 5'", "Contamination Estimate (contamMix)"), table = colnames(data))

colnames(data)[x] <- c("Avg Read Length (bp)", "PMD (%)", "Cont. Est.")

# Editting column values
data$Bone <- paste0(toupper(substr(data$Bone, 1, 1)), substr(data$Bone, 2, nchar(data$Bone)))
data$bone_other <- paste0(toupper(substr(data$`Bone Type`, 1, 1)), 
                          substr(data$`Bone Type`, 2, nchar(data$`Bone Type`)))
data$`Human Proportion (%)` <- data$`Human Proportion (%)`*100
data$`PMD (%)` <- data$`PMD (%)`*100

added_data <- data

added_data$Bone[added_data$Bone %in% "Unidentified Ossicle"] <- "UI Ossicle"

added_data$Bone <- factor(added_data$Bone, levels = c("Petrous", "Stapes", "Malleus", "Incus", "UI Ossicle"))

# Distribution of Human Proportion
hp <- added_data %>%
  select(Individual, Bone, `Human Proportion (%)`)

hp %>%
  ggplot(aes(x = Bone, y = `Human Proportion (%)`))+
  ylim(0, 100)+
  geom_boxplot()+
  geom_point()

# Running pairwise wilcoxon tests on bone types
x <- as.data.frame(pairwise.wilcox.test(hp$`Human Proportion (%)`, 
                                        hp$Bone, p.adjust.method =  "BH", paired = F)$p.value, paired = F)
x$comb <- row.names(x)

x <- melt(x, id.vars = "comb")
x <- x %>%
  filter(!is.na(value)) %>%
  filter(value <= 0.05) %>%
  mutate(comb = paste0(substr(comb, 1, 1), ".", substr(variable, 1, 1)))

x$comb


p.s <- paste0(sign(x$value[x$comb %in% "S.P"]));p.s
p.s <- paste0("***");p.s

p.u <- paste0(sign(x$value[x$comb %in% "U.P"]), " (p=", round(x$value[x$comb %in% "U.P"], 2), ")");p.u
p.u <- paste0(sign(x$value[x$comb %in% "U.P"]));p.u

s.m <- paste0(sign(x$value[x$comb %in% "M.S"]), " (p=", round(x$value[x$comb %in% "M.S"], 2), ")");s.m
s.m <- paste0(sign(x$value[x$comb %in% "M.S"]));s.m

p.i <- paste0(sign(x$value[x$comb %in% "I.P"]), " (p=", round(x$value[x$comb %in% "I.P"], 2), ")");p.i
p.i <- paste0(sign(x$value[x$comb %in% "I.P"]));p.i

s.i <- paste0(sign(x$value[x$comb %in% "I.S"]), " (p=", round(x$value[x$comb %in% "I.S"], 2), ")");s.i
s.i <- paste0(sign(x$value[x$comb %in% "I.S"]));s.i


# Temporary: Incus vs Malleus&Petrous
temp_test <- hp %>%
  filter(Bone %in% c("Incus", "Malleus", "Petrous"))

temp_test$comp <- ifelse(temp_test$Bone %nin% c("Incus"), "Non-Incus", "Incus")

wilcox.test(temp_test$`Human Proportion (%)`[temp_test$comp %in% "Incus"], 
            temp_test$`Human Proportion (%)`[temp_test$comp %in% "Non-Incus"], paired = F)$p.value


# Plot
library(ggbeeswarm)
p1 <- hp %>%
  mutate(facet = "Human Proportion (%)",
         facet2 = "All Data") %>%
  ggplot(aes(x = Bone, y = `Human Proportion (%)`, fill = Bone, color = Bone))+
  labs(y = "",
       x = "",
       title = "Change in sequence statistics after preparing library from other bones")+
  scale_y_continuous(limits = c(0,105), breaks = seq(0,100,20))+
  geom_boxplot(size = .5, outlier.color = NA, alpha = .85, width = .5, outliers = FALSE)+
  geom_quasirandom(aes(color = Bone), shape = 21, size = 2, stroke = .5, alpha = .6,
                   show.legend = F, width = 0.3) +
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  showSignificance(c(2,3), 80, -0.0001, s.m,
                   textParams = list(size = 4, fontface = "bold", vjust = .3), segmentParams = list(size = .75))+
  showSignificance(c(2,4), 85, -0.0001, s.i,
                   textParams = list(size = 4, fontface = "bold", vjust = .3), segmentParams = list(size = .75))+
  showSignificance(c(1,2), 90, -0.0001, p.s,
                   textParams = list(size = 4, fontface = "bold", vjust = .3), segmentParams = list(size = .75))+
  showSignificance(c(1,4), 95, -0.0001, p.i,
                   textParams = list(size = 4, fontface = "bold", vjust = .3), segmentParams = list(size = .75))+
  showSignificance(c(1,5), 100, -0.0001, p.u,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75))+
  facet_grid(facet~facet2, switch = "y")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p1


hp$comp <- ifelse(hp$Bone %in% "Stapes", "Stapes", "Non-Stapes")

wilcox.test(hp$`Human Proportion (%)`[hp$comp %in% "Stapes"], 
            hp$`Human Proportion (%)`[hp$comp %in% "Non-Stapes"])$p.value


# Distribution of Fragment  Length
avg <- added_data %>%
  select(Individual, Bone, `Avg Read Length (bp)`)

# Temporary script for the manuscript: Difference in fragment lengths stapes vs petrous
temp <- avg %>%
  filter(Bone %in% c("Stapes", "Malleus", "Incus", "Petrous"))

mean(temp$`Avg Read Length (bp)`[temp$Bone %in% "Stapes"])-mean(temp$`Avg Read Length (bp)`[temp$Bone %in% "Petrous"])
mean(temp$`Avg Read Length (bp)`[temp$Bone %in% "Incus"])-mean(temp$`Avg Read Length (bp)`[temp$Bone %in% "Petrous"])
mean(temp$`Avg Read Length (bp)`[temp$Bone %in% "Malleus"])-mean(temp$`Avg Read Length (bp)`[temp$Bone %in% "Petrous"])


wilcox.test(temp$`Avg Read Length (bp)`[temp$Bone %in% "Stapes"], temp$`Avg Read Length (bp)`[temp$Bone %in% "Petrous"])$p.value


avg %>%
  ggplot(aes(x = Bone, y = `Avg Read Length (bp)`))+
  ylim(0, 100)+
  geom_boxplot()+
  geom_point()

# Running pairwise wilcoxon tests on bone types
x <-  as.data.frame(pairwise.wilcox.test(avg$`Avg Read Length (bp)`, avg$Bone, 
                                         p.adj =  "BH", paired = F)$p.value)

x$comb <- row.names(x)

x <- melt(x, id.vars = "comb")
x <- x %>%
  filter(!is.na(value)) %>%
  filter(value <= 0.05) %>%
  mutate(comb = paste0(substr(comb, 1, 1), ".", substr(variable, 1, 1)))

x$comb


p2 <- avg %>%
  mutate(facet = "Mean Fragment Length (bp)") %>%
  ggplot(aes(x = Bone, y = `Avg Read Length (bp)`, fill = Bone, color = Bone))+
  labs(y = "",
       x = "",
       title = "Change in sequence statistics after preparing library from other bones")+
  ylim(40,100)+
  geom_boxplot(size = .5, outlier.color = NA, alpha = .85, width = .5, outliers = FALSE)+
  geom_boxplot(size = .5, outlier.color = NA, alpha = .85, width = .5, outliers = FALSE)+
  geom_quasirandom(aes(color = Bone), shape = 21, size = 2, stroke = .5, alpha = .6,
                   show.legend = F, width = 0.3) +
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  facet_grid(facet~., switch = "y")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank())+
  guides(fill = guide_legend(override.aes = list(size = .5)));p2


pmd <- added_data %>%
  select(Individual, Bone, `PMD (%)`)

pmd %>%
  ggplot(aes(x = Bone, y = `PMD (%)`))+
  ylim(0, 60)+
  geom_boxplot()+
  geom_point()

# Running pairwise wilcoxon tests on bone types
x <-  as.data.frame(pairwise.wilcox.test(pmd$`PMD (%)`, pmd$Bone, p.adjust.method =  "BH")$p.value)

x$comb <- row.names(x)

x <- melt(x, id.vars = "comb")
x <- x %>%
  filter(!is.na(value)) %>%
  filter(value <= 0.05) %>%
  mutate(comb = paste0(substr(comb, 1, 1), ".", substr(variable, 1, 1)))

p3 <- pmd %>%
  mutate(label = "PMD (%)") %>%
  ggplot(aes(x = Bone, y = `PMD (%)`, fill = Bone, color = Bone))+
  labs(y = "",
       x = "",
       title = "Change in sequence statistics after preparing library from other bones")+
  ylim(0,60)+
  geom_boxplot(size = .5, outlier.color = NA, alpha = .85, width = .5, outliers = FALSE)+
  geom_quasirandom(aes(color = Bone), shape = 21, size = 2, stroke = .5, alpha = .6,
                   show.legend = F, width = 0.3) +
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  facet_grid(label~., switch = "y")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank())+
  guides(fill = guide_legend(override.aes = list(size = .5)));p3


# Combine plots for D panel
p_3A <- ggarrange(p1, p2, p3, ncol = 1, nrow = 3, labels = c("A", "B", "C"));p_3A


##### _____Second panel ------
data <- read_xlsx("data/Supplemental_Tables_rev02.xlsx", sheet = "S1", skip = 1)
x <- match(x = c("Mean Read Length (BAM, bp)", 
                 "PMD 5'", "Contamination Estimate (contamMix)"), table = colnames(data))

colnames(data)[x] <- c("Avg Read Length (bp)", "PMD (%)", "Cont. Est.")

# Editting column values
data$Bone <- paste0(toupper(substr(data$Bone, 1, 1)), substr(data$Bone, 2, nchar(data$Bone)))
data$bone_other <- paste0(toupper(substr(data$`Bone Type`, 1, 1)), substr(data$`Bone Type`, 2, nchar(data$`Bone Type`)))
data$`Human Proportion (%)` <- data$`Human Proportion (%)`*100
data$`PMD (%)` <- data$`PMD (%)`*100

# HP ≥ 0
data_all <- data %>%
  filter(`Anthropological Age` %in% "Infant") %>%
  mutate(category = "Site 7 HP ≥ 0%")

# HP ≥ 25
data_25 <-  data %>%
  filter(`Anthropological Age` %in% "Infant") %>%
  filter(`Human Proportion (%)` >= 25) %>%
  mutate(category = "Site 7 HP ≥ 25%")

# HP ≥ 50
data_50 <- data %>%
  filter(`Anthropological Age` %in% "Infant") %>%
  filter(`Human Proportion (%)` >= 50) %>%
  mutate(category = "Site 7 HP ≥ 50%")

# Combine
data_all <- rbind(data_all, data_25, data_50)

added_data <- data_all

added_data$Bone[added_data$Bone %in% "Unidentified Ossicle"] <- "UI Ossicle"

added_data$Bone <- factor(added_data$Bone, levels = c("Petrous", "Stapes", "Malleus", "Incus", "UI Ossicle"))

# Distribution of Human Proportion
hp <- added_data %>%
  select(Individual, Bone, `Human Proportion (%)`, category) %>%
  mutate(yaxis = "Human Proportion (%)")

hp %>%
  ggplot(aes(x = Bone, y = `Human Proportion (%)`))+
  ylim(0, 100)+
  geom_boxplot()+
  geom_point()


# Running Statistical Tests for HP ≥ 0%
hp_0 <- hp %>%
  filter(category %in% "Site 7 HP ≥ 0%")

# Running pairwise wilcoxon tests on bone types
x <- as.data.frame(pairwise.wilcox.test(hp_0$`Human Proportion (%)`, 
                                        hp_0$Bone, p.adjust.method =  "BH", paired = F)$p.value, paired = F)
x$comb <- row.names(x)

x <- melt(x, id.vars = "comb")
x <- x %>%
  filter(!is.na(value)) %>%
  filter(value <= 0.05) %>%
  mutate(comb = paste0(substr(comb, 1, 1), ".", substr(variable, 1, 1)))

x$comb

# Assigning significance stars
p.s_hp0 <- paste0(sign(x$value[x$comb %in% "S.P"]));p.s_hp0
p.s_hp0 <- paste0("***");p.s_hp0

s.m_hp0 <- paste0(sign(x$value[x$comb %in% "M.S"]), " (p=", round(x$value[x$comb %in% "M.S"], 2), ")");s.m_hp0
s.m_hp0 <- paste0(sign(x$value[x$comb %in% "M.S"]));s.m_hp0


# Running Statistical Tests for HP ≥ 25% 
hp_25 <- hp %>%
  filter(category %in% "Site 7 HP ≥ 25%")

# Running pairwise wilcoxon tests on bone types
x <- as.data.frame(pairwise.wilcox.test(hp_25$`Human Proportion (%)`, 
                                        hp_25$Bone, p.adjust.method =  "BH", paired = F)$p.value, paired = F)
x$comb <- row.names(x)

x <- melt(x, id.vars = "comb")
x <- x %>%
  filter(!is.na(value)) %>%
  filter(value <= 0.05) %>%
  mutate(comb = paste0(substr(comb, 1, 1), ".", substr(variable, 1, 1)))

x$comb

# Assigning significance stars
p.s_hp25 <- paste0(sign(x$value[x$comb %in% "S.P"]));p.s_hp25
p.s_hp25 <- paste0("***");p.s_hp25


# Running Statistical Tests for HP ≥ 50%
hp_50 <- hp %>%
  filter(category %in% "Site 7 HP ≥ 50%")

# Running pairwise wilcoxon tests on bone types
x <- as.data.frame(pairwise.wilcox.test(hp_50$`Human Proportion (%)`, 
                                        hp_50$Bone, p.adjust.method =  "BH", paired = F)$p.value, paired = F)
x$comb <- row.names(x)

x <- melt(x, id.vars = "comb")
x <- x %>%
  filter(!is.na(value)) %>%
  filter(value <= 0.05) %>%
  mutate(comb = paste0(substr(comb, 1, 1), ".", substr(variable, 1, 1)))

x$comb



p1 <- hp %>%
  ggplot(aes(x = Bone, y = `Human Proportion (%)`, fill = Bone, color = Bone))+
  labs(y = "",
       x = "",
       title = "Change in sequence statistics after preparing library from other bones")+
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20))+
  geom_boxplot(size = .5, outlier.color = NA, alpha = .85, width = .5)+
  geom_quasirandom(aes(color = Bone), shape = 21, size = 2, stroke = .5, alpha = .6,
                   show.legend = F, width = 0.3) +
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  showSignificance(c(1,2), 80, -0.0001, p.s_hp0,
                   textParams = list(size = 4, fontface = "bold", vjust = .25), segmentParams = list(size = .75), 
                   panel = list(category = "Site 7 HP ≥ 0%"))+
  showSignificance(c(2,3), 85, -0.0001, s.m_hp0,
                   textParams = list(size = 4, fontface = "bold", vjust = .25), segmentParams = list(size = .75), 
                   panel = list(category = "Site 7 HP ≥ 0%"))+
  showSignificance(c(1,2), 80, -0.0001, p.s_hp25,
                   textParams = list(size = 4, fontface = "bold", vjust = .25), segmentParams = list(size = .75), 
                   panel = list(category = "Site 7 HP ≥ 25%"))+
  facet_grid(yaxis~category, switch = "y")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p1




# Comparison on Mean Fragment Length 

# Distribution of Human Proportion
avg <- added_data %>%
  select(Individual, Bone, `Avg Read Length (bp)`, category) %>%
  mutate(yaxis = "Mean Fragment Length (bp)")

avg %>%
  ggplot(aes(x = Bone, y = `Avg Read Length (bp)`))+
  ylim(0, 100)+
  geom_boxplot()+
  geom_point()


# Running Statistical Tests for HP ≥ 0%
hp_0 <- avg %>%
  filter(category %in% "Site 7 HP ≥ 0%")

# Running pairwise wilcoxon tests on bone types
x <- as.data.frame(pairwise.wilcox.test(hp_0$`Avg Read Length (bp)`, 
                                        hp_0$Bone, p.adjust.method =  "BH", paired = F)$p.value, paired = F)
x$comb <- row.names(x)

x <- melt(x, id.vars = "comb")
x <- x %>%
  filter(!is.na(value)) %>%
  filter(value <= 0.05) %>%
  mutate(comb = paste0(substr(comb, 1, 1), ".", substr(variable, 1, 1)))

x$comb

# Running Statistical Tests for HP ≥ 25%
hp_25 <- avg %>%
  filter(category %in% "Site 7 HP ≥ 25%")

# Running pairwise wilcoxon tests on bone types
x <- as.data.frame(pairwise.wilcox.test(hp_25$`Avg Read Length (bp)`, 
                                        hp_25$Bone, p.adjust.method =  "BH", paired = F)$p.value, paired = F)
x$comb <- row.names(x)

x <- melt(x, id.vars = "comb")
x <- x %>%
  filter(!is.na(value)) %>%
  filter(value <= 0.05) %>%
  mutate(comb = paste0(substr(comb, 1, 1), ".", substr(variable, 1, 1)))

x$comb


# Running Statistical Tests for HP ≥ 50%
hp_50 <- avg %>%
  filter(category %in% "Site 7 HP ≥ 50%")

# Running pairwise wilcoxon tests on bone types
x <- as.data.frame(pairwise.wilcox.test(hp_50$`Avg Read Length (bp)`, 
                                        hp_50$Bone, p.adjust.method =  "BH", paired = F)$p.value, paired = F)
x$comb <- row.names(x)

x <- melt(x, id.vars = "comb")
x <- x %>%
  filter(!is.na(value)) %>%
  filter(value <= 0.05) %>%
  mutate(comb = paste0(substr(comb, 1, 1), ".", substr(variable, 1, 1)))

x$comb


p2 <- avg %>%
  ggplot(aes(x = Bone, y = `Avg Read Length (bp)`, fill = Bone, color = Bone))+
  labs(y = "",
       x = "",
       title = "Change in sequence statistics after preparing library from other bones")+
  scale_y_continuous(limits = c(40,100), breaks = seq(40,100,20))+
  geom_boxplot(size = .5, outlier.color = NA, alpha = .85, width = .5)+
  geom_quasirandom(aes(color = Bone), shape = 21, size = 2, stroke = .5, alpha = .6,
                   show.legend = F, width = 0.3) +
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  facet_grid(yaxis~category, switch = "y")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text.y = element_text(face = "bold", size = 11,
                                    color = "white"),
        strip.text.x = element_blank(),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p2


# Comparison on Post-mortem Damage 

# Distribution of Human Proportion
pmd <- added_data %>%
  select(Individual, Bone, `PMD (%)`, category) %>%
  mutate(yaxis = "PMD (%)")

pmd %>%
  ggplot(aes(x = Bone, y = `PMD (%)`))+
  ylim(0, 100)+
  geom_boxplot()+
  geom_point()


# Running Statistical Tests for HP ≥ 0%
hp_0 <- pmd %>%
  filter(category %in% "Site 7 HP ≥ 0%")

# Running pairwise wilcoxon tests on bone types
x <- as.data.frame(pairwise.wilcox.test(hp_0$`PMD (%)`, 
                                        hp_0$Bone, p.adjust.method =  "BH", paired = F)$p.value, paired = F)
x$comb <- row.names(x)

x <- melt(x, id.vars = "comb")
x <- x %>%
  filter(!is.na(value)) %>%
  filter(value <= 0.05) %>%
  mutate(comb = paste0(substr(comb, 1, 1), ".", substr(variable, 1, 1)))

x$comb

# Running Statistical Tests for HP ≥ 25%
hp_25 <- pmd %>%
  filter(category %in% "Site 7 HP ≥ 25%")

# Running pairwise wilcoxon tests on bone types
x <- as.data.frame(pairwise.wilcox.test(hp_25$`PMD (%)`, 
                                        hp_25$Bone, p.adjust.method =  "BH", paired = F)$p.value, paired = F)
x$comb <- row.names(x)

x <- melt(x, id.vars = "comb")
x <- x %>%
  filter(!is.na(value)) %>%
  filter(value <= 0.05) %>%
  mutate(comb = paste0(substr(comb, 1, 1), ".", substr(variable, 1, 1)))

x$comb


# Running Statistical Tests for HP ≥ 50%
hp_50 <- pmd %>%
  filter(category %in% "Site 7 HP ≥ 50%")

# Running pairwise wilcoxon tests on bone types
x <- as.data.frame(pairwise.wilcox.test(hp_50$`PMD (%)`, 
                                        hp_50$Bone, p.adjust.method =  "BH", paired = F)$p.value, paired = F)
x$comb <- row.names(x)

x <- melt(x, id.vars = "comb")
x <- x %>%
  filter(!is.na(value)) %>%
  filter(value <= 0.05) %>%
  mutate(comb = paste0(substr(comb, 1, 1), ".", substr(variable, 1, 1)))

x$comb


p3 <- pmd %>%
  ggplot(aes(x = Bone, y = `PMD (%)`, fill = Bone, color = Bone))+
  labs(y = "",
       x = "",
       title = "Change in sequence statistics after preparing library from other bones")+
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,20))+
  geom_boxplot(size = .5, outlier.color = NA, alpha = .85, width = .5)+
  geom_quasirandom(aes(color = Bone), shape = 21, size = 2, stroke = .5, alpha = .6,
                   show.legend = F, width = 0.3) +
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  facet_grid(yaxis~category, switch = "y")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text.y = element_text(face = "bold", size = 11,
                                    color = "white"),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p3




# Saving
p_3B <- ggarrange(p1, p2, p3, ncol = 1, nrow = 3, labels = c("D", "E", "F"));p_3B


# Combining both panels for Main Figure 3
p_3 <- ggarrange(p_3A, p_3B, ncol = 2, nrow = 1, align = "v", widths = c(.5,1));p_3
ggsave("figures/png/figure3.png", p_3, device = png, height = 10, width = 20, dpi = 600)
ggsave("figures/svg/figure3.svg", p_3, device = svg, height = 10, width = 20, dpi = 600)


# Supplemental Figure 2 --------------------------------------------------
# Microbial Analysis

# Comparison of Non-Oral and Oral relative abundances between bone types

# Reading data
biplot <- read.table("data/source_comparison", header = T, sep = "\t", check.names = F)
biplot$`Tissue Type` <- factor(biplot$`Tissue Type`, 
                               levels = c("Petrous", "Stapes", "Malleus", "Incus", "Unidentified Ossicle", "Tooth"))

#### ____First Panel -----
##### Oral 
x <- as.data.frame(pairwise.wilcox.test(biplot$Oral, biplot$`Tissue Type`, 
                                        p.adjust.method =  "BH", paired = F)$p.value, paired = F)
x$comb <- row.names(x)

#Filtering p values higher than 0.05
x <- melt(x, id.vars = "comb")
x <- x %>%
  filter(!is.na(value)) %>%
  filter(value <= 0.05) %>%
  mutate(comb = paste0(substr(comb, 1, 1), ".", substr(variable, 1, 1)));x$comb

# Adjusting labels
t.p <- paste0(sign(x$value[x$comb %in% "T.P"]));t.p
t.p <- "***"
t.s <- paste0(sign(x$value[x$comb %in% "T.P"]));t.s
t.s <- "***"
t.m <- paste0(sign(x$value[x$comb %in% "T.P"]));t.m
t.m <- "***"
t.i <- paste0(sign(x$value[x$comb %in% "T.P"]));t.i
t.i <- "***"



# Plotting
pb_1 <- biplot %>%
  ggplot()+
  geom_boxplot(aes(x = `Tissue Type`, y = Oral, fill = `Tissue Type`, color = `Tissue Type`), 
               outliers = FALSE, show.legend = F)+
  geom_quasirandom(aes(x = `Tissue Type`, y = Oral, fill = `Tissue Type`),
                   shape = 21, size = 3, stroke = .75, width = 0.2, alpha = .8)+
  #geom_point(position = position_jitter(width = .2, height = .01))+
  scale_fill_manual(values = palette, labels = c("Petrous (n=105)", "Stapes (n=14)", "Malleus (n=13)",
                                                 "Incus (n=10)", "UI Ossicle (n=5)", "Tooth (n=5)"))+
  scale_color_manual(values = c("black","black","black", "black", "black",  "black"),
                     labels = c("Petrous (n=105)", "Stapes (n=14)", "Malleus (n=13)",
                                "Incus (n=10)", "UI Ossicle (n=5)", "Tooth (n=5)"))+
  ylim(0,10)+
  labs(x = "",
       y = "Oral Reads (%)")+
  showSignificance(c(1,6), 6, -0.0001, t.p,
                   textParams = list(size = 4, fontface = "bold", vjust = .25), segmentParams = list(size = .75))+
  showSignificance(c(2,6), 7, -0.0001, t.s,
                   textParams = list(size = 4, fontface = "bold", vjust = .25), segmentParams = list(size = .75))+
  showSignificance(c(3,6), 8, -0.0001, t.m,
                   textParams = list(size = 4, fontface = "bold", vjust = .25), segmentParams = list(size = .75))+
  showSignificance(c(4,6), 9, -0.0001, t.i,
                   textParams = list(size = 4, fontface = "bold", vjust = .25), segmentParams = list(size = .75))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(),
        legend.direction = "horizontal",
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank())+
  guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 1)));pb_1


#### ____Second Panel -----
# Environmental 
median(biplot$`Non-oral`[biplot$`Tissue Type` %in% "Tooth"])
median(biplot$`Non-oral`[biplot$`Tissue Type` %in% "Petrous"])
median(biplot$`Non-oral`[biplot$`Tissue Type` %in% "Malleus"])
median(biplot$`Non-oral`[biplot$`Tissue Type` %in% "Stapes"])
median(biplot$`Non-oral`[biplot$`Tissue Type` %in% "Incus"])

median(biplot$`Non-oral`[biplot$`Tissue Type` %nin% c("Petrous", "Tooth")])


x <- as.data.frame(pairwise.wilcox.test(biplot$`Non-oral`, biplot$`Tissue Type`, 
                                        p.adjust.method =  "BH", paired = F)$p.value, paired = F)
x$comb <- row.names(x)

#Filtering p values higher than 0.05
x <- melt(x, id.vars = "comb")
x <- x %>%
  filter(!is.na(value)) %>%
  filter(value <= 0.05) %>%
  mutate(comb = paste0(substr(comb, 1, 1), ".", substr(variable, 1, 1)));x$comb


# Adjusting labels
p.s <- paste0(sign(x$value[x$comb %in% "S.P"]));p.s
p.s <- "***"
p.i <- paste0(sign(x$value[x$comb %in% "I.P"]));p.i
t.p <- paste0(sign(x$value[x$comb %in% "T.P"]));t.p
t.s <- paste0(sign(x$value[x$comb %in% "T.S"]));t.s
t.m <- paste0(sign(x$value[x$comb %in% "T.M"]));t.m
t.i <- paste0(sign(x$value[x$comb %in% "T.I"]));t.i
t.u <- paste0(sign(x$value[x$comb %in% "T.U"]));t.u

# Plotting
table(biplot$`Tissue Type`)

pb_2 <- biplot %>%
  ggplot()+
  geom_boxplot(aes(x = `Tissue Type`, y = `Non-oral`, fill = `Tissue Type`, color = `Tissue Type`), 
               outliers = FALSE, show.legend = F)+
  geom_quasirandom(aes(x = `Tissue Type`, y = `Non-oral`, fill = `Tissue Type`),
                   shape = 21, size = 3, stroke = .75, width = 0.2, alpha = .8)+
  scale_fill_manual(values = palette, labels = c("Petrous (n=101)", "Stapes (n=14)", "Malleus (n=13)",
                                                 "Incus (n=10)", "UI Ossicle (n=5)", "Tooth (n=5)"))+
  scale_color_manual(values = c("black","black","black", "black", "black",  "black"),
                     labels = c("Petrous (n=101)", "Stapes (n=14)", "Malleus (n=13)",
                                "Incus (n=10)", "UI Ossicle (n=5)", "Tooth (n=5)"))+
  ylim(0,34)+
  labs(x = "",
       y = "Non-Oral Reads (%)")+
  showSignificance(c(1,2), 22, -0.0001, p.s,
                   textParams = list(size = 4, fontface = "bold", vjust = .25), segmentParams = list(size = .75))+
  showSignificance(c(1,4), 24, -0.0001, p.i,
                   textParams = list(size = 4, fontface = "bold", vjust = .25), segmentParams = list(size = .75))+
  showSignificance(c(1,6), 26, -0.0001, t.p,
                   textParams = list(size = 4, fontface = "bold", vjust = .25), segmentParams = list(size = .75))+
  showSignificance(c(2,6), 28, -0.0001, t.s,
                   textParams = list(size = 4, fontface = "bold", vjust = .25), segmentParams = list(size = .75))+
  showSignificance(c(3,6), 30, -0.0001, t.m,
                   textParams = list(size = 4, fontface = "bold", vjust = .25), segmentParams = list(size = .75))+
  showSignificance(c(4,6), 32, -0.0001, t.i,
                   textParams = list(size = 4, fontface = "bold", vjust = .25), segmentParams = list(size = .75))+
  showSignificance(c(5,6), 34, -0.0001, t.u,
                   textParams = list(size = 4, fontface = "bold", vjust = .25), segmentParams = list(size = .75))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank())+
  guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 1)));pb_2


combined <- ggarrange(pb_1, pb_2, ncol = 2, nrow = 1, labels = c("A", "B"), common.legend = T);combined

# Saving
ggsave("figures/svg/supplemental_figure_2.svg", device = svg, width = 15, height = 7)
ggsave("figures/png/supplemental_figure_2.png", device = "png", height = 7, width = 15, dpi = 600)


# Supplemental Figure 3 --------------------------------------------------
# Comparison of mean fragment lengths of raw and mapped reads
# Loading data
data <- read_xlsx("data/Supplemental_Tables_rev02.xlsx", sheet = "S1", skip = 1)

data_f <- data %>%
  filter(Region %in% "Anatolia") %>%
  select(Individual, Bone, `Mean Read Length (BAM, bp)`, `Mean Read Length (FASTQ, bp)`)

# Change column names
colnames(data_f) <- c("sample", "bone", "bam", "fastq")

df_long <- data_f %>%
  select(sample, bone, bam, fastq) %>%
  pivot_longer(
    cols = c(bam, fastq),
    names_to = "file_type",
    values_to = "avg_len")


# Define bone order
bone_order <- c("Stapes", "Malleus", "Incus", "Unidentified Ossicle", "Petrous", "Tooth")
df_long$bone <- factor(df_long$bone, levels = bone_order)



ggplot(df_long, aes(x = avg_len, fill = file_type)) +
  geom_density(alpha = 0.5, color = "black") +
  facet_wrap(~bone, scales = "free") +
  scale_fill_manual(values = c("fastq" = "#0d3b66", "bam" = "#faf0ca"),
                    name = "File Type",
                    labels = c("Mapped (BAM)", "Raw (FASTQ)")) +
  labs(x = "Mean Fragment Length (bp)",
       y = "Density") +
  scale_x_continuous(limits = c(20,130), breaks = seq(20,130,15))+
  theme_bw() +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 11, color = "white"),
        strip.background = element_rect(fill = "#121212"),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 11),
        axis.text = element_text(face = "bold", size = 10),
        axis.title = element_text(face = "bold", size = 11),
        plot.caption = element_text(face = "italic", size = 8))

# Saving
ggsave(device = png, filename = "figures/png/supplemental_figure_3.png", height = 7, width = 12, dpi = 600)
ggsave(device = svg, filename = "figures/svg/supplemental_figure_3.svg", height = 7, width = 12)