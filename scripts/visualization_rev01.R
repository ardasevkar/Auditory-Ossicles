#This code is for visualization of summary statistics for ossicle and petrous libraries

# Setup -------------------------------------------------------------------
Sys.setlocale(category = 'LC_ALL','en_US.UTF-8')
options(timeout=1000, scipen = 999)
`%nin%` = Negate(`%in%`)
sign <- function(sayı){
  
  x <- ifelse(sayı <= 0.05 & sayı > 0.01, "*",
              ifelse(sayı <= 0.01 & sayı > 0.001 , "**",
                     ifelse(sayı <= 0.001, "*** (p ≤ 0.001)", "n.s.")))
  return(x)}


# Loading Libraries -------------------------------------------------------
library(raster)
library(ggplot2)
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





# Figure 1 ----------------------------------------------------------------
# The code below includes the scripts used to generate each panel of Figure 1 in the manuscript.


##### First Panel ----------------------------------------------------------------
img <- readPNG("figures/ossicles.png")

# Convert to rasterGrob
grob_img <- rasterGrob(img, interpolate = TRUE)

# Wrap into a ggplot object
p_A <- ggplot() +
  annotation_custom(grob_img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void()+
  theme(panel.background = element_rect(fill = "white", color = "white"),
        plot.margin=unit(c(.5,0,0,0), 'cm'));p_A

#Second Panel
hp <- read.table("data/seqstats.txt", header = T, sep = "\t")
hp <- hp %>%
  filter(variable %in% "Human Proportion (%)") %>%
  as.data.frame()

hp$Bone[hp$Bone %in% "Unidentified Ossicle"] <- "UI Ossicle"
hp$comb[hp$comb %in% "Petrous|Unidentified Ossicle"] <- "Petrous|UI Ossicle"

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
  #star <- paste0(sign(test), " (p=", round(test, digits = 3), ")");star
  star <- sign(test);star
  assign(paste0("star_", i), star)
}

star_1;star_2;star_3;star_4

##### Second Panel ----------------------------------------------------------------
p_B <- hp %>%
  ggplot(aes(x = Bone, y = value))+
  ylim(0,85)+
  labs(x = "", y = "Human Proportion (%)")+
  geom_line(aes(group = Sample), color = "#121212", size = .75, show.legend = F)+
  geom_point(aes(fill = Bone, color = Bone), shape = 21, size = 5.5, show.legend = T, stroke = .75)+
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("#121212","grey70","#121212", "#121212", "#121212", "#121212"))+
  showSignificance(c(1,2), 82, -0.0001, star_1,
                   textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
                   panel = list(comb = 1))+
  showSignificance(c(1,2), 82, -0.0001, star_3,
                   textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
                   panel = list(comb = 3))+
  facet_grid(.~comb, space = "free", scales = "free_x")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_blank(),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank(),
        plot.margin=unit(c(.5,0,0,0), 'cm'));p_B

##### Third Panel ----------------------------------------------------------------
third <- all %>%
  filter(Bone %nin% "Petrous") %>%
  filter(variable %in% c("Human Proportion (%)", "Clonality (%)", 
                         "Avg Read Length (bp)", "PMD (%)", "Cont. Est."))

third$variable <- as.character(third$variable)
third$variable[third$variable %in% "Human Proportion (%)"] <- "Human Proportion"
third$variable[third$variable %in% "Clonality (%)"] <- "Clonality"
third$variable[third$variable %in% "Avg Read Length (bp)"] <- "Mean Fragment Length"
third$variable[third$variable %in% "Cont. Est."] <- "Contamination"
third$variable[third$variable %in% "PMD (%)"] <- "PMD"

third$variable <- factor(third$variable, levels = c("Human Proportion", "Mean Fragment Length", 
                                                    "PMD", "Clonality", "Contamination"))
p_C <- third %>%
  ggplot(aes(x = Bone, y = diff))+
  geom_hline(yintercept = 0, color = "#ff2800", linetype = "dashed", size = 1)+
  geom_boxplot(aes(fill = Bone, color = Bone), outliers = FALSE, width = .5)+
  labs(y = "Change Relative to Petrous (%)",
       x = "")+
  geom_point(aes(fill = Bone, color = Bone), shape = 21, size = 4, position = position_jitterdodge(jitter.width = 1))+
  facet_wrap(.~variable, scales = "free", nrow = 1)+
  facetted_pos_scales(y = list(scale_y_continuous(breaks = c(0, 100, 200, 300, 500, 750, 1500))))+
  scale_fill_manual(values = c(col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("grey70","black", "black", "black", "black"))+
  scale_x_discrete(labels = c("Stapes", "Malleus", "Incus","UI Ossicle", "Tooth"))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 13,
                                  color = "white"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p_C



hp$comp <- ifelse(hp$Bone %nin% c("Petrous"), "Ossicle", "Non-Ossicle")

wilcox.test(hp$value[hp$comp %in% "Ossicle"], 
            hp$value[hp$comp %in% "Non-Ossicle"], paired = T)$p.value


##### Fourth Panel ------------------------------------------------------------
# Loading dataset
data <- read_xlsx("data/Supplementary_Tables_rev01.xlsx", sheet = "S1")

# Remove tooth
data <- data %>%
  filter(Bone %nin% "Tooth") %>%
  filter(`Library Type` %nin% "petrous&tooth")

table(data$Bone)

x <- match(x = c("Mean Read Length (bp)", 
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
p1 <- hp %>%
  ggplot(aes(x = Bone, y = `Human Proportion (%)`, fill = Bone, color = Bone))+
  labs(y = "Human Proportion (%)",
       x = "",
       title = "Change in sequence statistics after preparing library from other bones")+
  scale_y_continuous(limits = c(0,105), breaks = seq(0,100,20))+
  geom_boxplot(size = .5, outlier.color = NA, alpha = .85, width = .5, outliers = FALSE)+
  geom_point(aes(color = Bone), shape = 21, position = position_jitterdodge(jitter.width = 1), 
             size = 3, stroke = .75, show.legend = F)+
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  showSignificance(c(2,3), 80, -0.0001, s.m,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75))+
  showSignificance(c(2,4), 85, -0.0001, s.i,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75))+
  showSignificance(c(1,2), 90, -0.0001, p.s,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75))+
  showSignificance(c(1,4), 95, -0.0001, p.i,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75))+
  showSignificance(c(1,5), 100, -0.0001, p.u,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
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
  ggplot(aes(x = Bone, y = `Avg Read Length (bp)`, fill = Bone, color = Bone))+
  labs(y = "Mean Fragment Length (bp)",
       x = "",
       title = "Change in sequence statistics after preparing library from other bones")+
  ylim(40,100)+
  geom_boxplot(size = .5, outlier.color = NA, alpha = .85, width = .5, outliers = FALSE)+
  geom_point(aes(color = Bone), shape = 21, position = position_jitterdodge(jitter.width = 1), 
             size = 3, stroke = .75, show.legend = F)+
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
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
  ggplot(aes(x = Bone, y = `PMD (%)`, fill = Bone, color = Bone))+
  labs(y = "PMD (%)",
       x = "",
       title = "Change in sequence statistics after preparing library from other bones")+
  ylim(0,60)+
  geom_boxplot(size = .5, outlier.color = NA, alpha = .85, width = .5, outliers = FALSE )+
  geom_point(aes(color = Bone), shape = 21, position = position_jitterdodge(jitter.width = 1), 
             size = 3, stroke = .75, show.legend = F)+
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank())+
  guides(fill = guide_legend(override.aes = list(size = .5)));p3


# Combine plots for D panel
p_D <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1);p_D


# Merging plots all panels together
top_panel <- ggarrange(p_A, NULL, p_B, ncol = 3, nrow = 1, widths = c(.35,0.005,1), heights = c(.35,0.005,1),
                       labels = c("A","", "B"), align = "v", font.label = list(size = 20, face = "bold"));top_panel

bottom_panel <- ggarrange(p_C, p_D, ncol = 1, nrow = 2, labels = c("C", "D"),
                          font.label = list(size = 20, face = "bold"));bottom_panel

all_plot <- ggarrange(top_panel, NULL, bottom_panel, nrow = 3, ncol = 1, heights = c(.45, 0.02 ,1));all_plot

ggsave("figures/png/figure_1.png", device = "png", width = 23, height = 15, dpi = 300)
ggsave("figures/svg/figure_1.svg", device = svg, width = 23, height = 15)



# Supplementary Figure 2 --------------------------------------------------
# Code below contains summary plots of sequence statistics of ossicle and petrous libraries prepared from the same individual.

# Loading data
seqstats <- read.table("data/seqstats.txt", header = T, sep = "\t")
seqstats$Bone[seqstats$Bone %in% "Unidentified Ossicle"] <- "UI Ossicle"
seqstats$Bone <- factor(seqstats$Bone, levels = c("Petrous", "Stapes", "Malleus", "Incus", "UI Ossicle", "Tooth"))
seqstats$comb[seqstats$comb %in% "Petrous|Unidentified Ossicle"] <- "Petrous|UI Ossicle"

##### Distribution of Human Proportion ----------------------------------------
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



p1 <- hp %>%
  ggplot(aes(x = Bone, y = value))+
  ylim(0,90)+
  labs(x = "", y = "Human Proportion (%)")+
  geom_line(aes(group = Sample), color = "#121212", size = 1, show.legend = F)+
  geom_point(aes(fill = Bone, color = Bone), shape = 21, size = 5.5, show.legend = T, stroke = .75)+
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
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p1


##### Distribution of Mean Fragment Length ----------------------------------------
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


p2 <- avg %>%
  ggplot(aes(x = Bone, y = value))+
  ylim(40,100)+
  labs(x = "", y = "Mean Fragment Length (bp)")+
  geom_line(aes(group = Sample), color = "#121212", size = 1, show.legend = F)+
  geom_point(aes(fill = Bone, color = Bone), shape = 21, size = 5.5, show.legend = T, stroke = .75)+
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  showSignificance(c(1,2), 95, -0.0001, star_1,
                   textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
                   panel = list(comb = 1))+
  showSignificance(c(1,2), 92, -0.0001, star_2,
                   textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
                   panel = list(comb = 2))+
  # showSignificance(c(1,2), 92, -0.0001, star_3,
  #                  textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
  #                  panel = list(comb = 3))+
  # showSignificance(c(1,2), 92, -0.0001, star_4,
  #                  textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
  #                  panel = list(comb = 4))+
  facet_grid(.~comb, space = "free", scales = "free_x")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_blank(),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p2


##### Distribution of Postmortem Damage----------------------------------------
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

p3 <- pmd %>%
  ggplot(aes(x = Bone, y = value))+
  ylim(0,55)+
  labs(x = "", y = "PMD (%)")+
  geom_line(aes(group = Sample), color = "#121212", size = 1, show.legend = F)+
  geom_point(aes(fill = Bone, color = Bone), shape = 21, size = 5.5, show.legend = T, stroke = .75)+
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui, col_t))+
  # showSignificance(c(1,2), 50, -0.0001, star_1,
  #                  textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
  #                  panel = list(comb = 1))+
  # showSignificance(c(1,2), 50, -0.0001, star_2,
  #                  textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
  #                  panel = list(comb = 2))+
  # showSignificance(c(1,2), 50, -0.0001, star_3,
  #                  textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
  #                  panel = list(comb = 3))+
  # showSignificance(c(1,2), 50, -0.0001, star_4,
  #                  textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
  #                  panel = list(comb = 4))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
facet_grid(.~comb, space = "free", scales = "free_x")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_blank(),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p3


##### Distribution of Clonality Levels of Downsampled BAM Files ----------------------------------------------------
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


p4 <- clon %>%
  ggplot(aes(x = Bone, y = value))+
  # ylim(0,23)+
  labs(x = "", y = "Clonality (%)")+
  geom_line(aes(group = Sample), color = "#121212", size = 1, show.legend = F)+
  geom_point(aes(fill = Bone, color = Bone), shape = 21, size = 5.5, show.legend = T, stroke = .75)+
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  # showSignificance(c(1,2), 37, -0.0001, star_1,
  #                  textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
  #                  panel = list(comb = 1))+
  # showSignificance(c(1,2), 37, -0.0001, star_2,
  #                  textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
  #                  panel = list(comb = 2))+
  # showSignificance(c(1,2), 37, -0.0001, star_3,
  #                  textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
  #                  panel = list(comb = 3))+
  # showSignificance(c(1,2), 37, -0.0001, star_4,
  #                  textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
  #                  panel = list(comb = 4))+
facet_grid(.~comb, space = "free", scales = "free_x")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_blank(),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p4

ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 1, nrow = 4)

ggsave("figures/png/supplementary_figure_2.png", device = "png", width = 10, height = 15, dpi = 600)
ggsave("figures/svg/supplementary_figure_2.svg", device = svg, width = 10, height = 15)


# Supplementary Figure 3 --------------------------------------------------
# Figure for filtering only Site 7 individuals
data <- read_xlsx("data/Supplementary_Tables_rev01.xlsx", sheet = "S1")
x <- match(x = c("Mean Read Length (bp)", 
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
  mutate(category = "HP ≥ 0%")

# HP ≥ 25
data_25 <-  data %>%
  filter(`Anthropological Age` %in% "Infant") %>%
  filter(`Human Proportion (%)` >= 25) %>%
  mutate(category = "HP ≥ 25%")

# HP ≥ 50
data_50 <- data %>%
  filter(`Anthropological Age` %in% "Infant") %>%
  filter(`Human Proportion (%)` >= 50) %>%
  mutate(category = "HP ≥ 50%")

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


##### Running Statistical Tests for HP ≥ 0% ------
hp_0 <- hp %>%
  filter(category %in% "HP ≥ 0%")

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


##### Running Statistical Tests for HP ≥ 25% -----
hp_25 <- hp %>%
  filter(category %in% "HP ≥ 25%")

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


##### Running Statistical Tests for HP ≥ 50% -----
hp_50 <- hp %>%
  filter(category %in% "HP ≥ 50%")

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
  scale_y_continuous(limits = c(0,90), breaks = seq(0,100,20))+
  geom_boxplot(size = .5, outlier.color = NA, alpha = .85, width = .5)+
  geom_point(aes(color = Bone), shape = 21, position = position_jitterdodge(jitter.width = 1), 
             size = 3, stroke = .75, show.legend = F)+
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  showSignificance(c(1,2), 80, -0.0001, p.s_hp0,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75), 
                   panel = list(category = "HP ≥ 0%"))+
  showSignificance(c(2,3), 85, -0.0001, s.m_hp0,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75), 
                   panel = list(category = "HP ≥ 0%"))+
  showSignificance(c(1,2), 80, -0.0001, p.s_hp25,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75), 
                   panel = list(category = "HP ≥ 25%"))+
  facet_grid(yaxis~category, switch = "y")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p1




##### Comparison on Mean Fragment Length -----

# Distribution of Human Proportion
avg <- added_data %>%
  select(Individual, Bone, `Avg Read Length (bp)`, category) %>%
  mutate(yaxis = "Mean Fragment Length")

avg %>%
  ggplot(aes(x = Bone, y = `Avg Read Length (bp)`))+
  ylim(0, 100)+
  geom_boxplot()+
  geom_point()


# Running Statistical Tests for HP ≥ 0%
hp_0 <- avg %>%
  filter(category %in% "HP ≥ 0%")

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
  filter(category %in% "HP ≥ 25%")

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
  filter(category %in% "HP ≥ 50%")

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
  geom_point(aes(color = Bone), shape = 21, position = position_jitterdodge(jitter.width = 1), 
             size = 3, stroke = .75, show.legend = F)+
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  facet_grid(yaxis~category, switch = "y")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text.y = element_text(face = "bold", size = 11,
                                    color = "white"),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p2


# Comparison on Post-mortem Damage--------------------------------------

# Distribution of Human Proportion
pmd <- added_data %>%
  select(Individual, Bone, `PMD (%)`, category) %>%
  mutate(yaxis = "PMD")

pmd %>%
  ggplot(aes(x = Bone, y = `PMD (%)`))+
  ylim(0, 100)+
  geom_boxplot()+
  geom_point()


# Running Statistical Tests for HP ≥ 0%
hp_0 <- pmd %>%
  filter(category %in% "HP ≥ 0%")

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
  filter(category %in% "HP ≥ 25%")

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
  filter(category %in% "HP ≥ 50%")

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
  geom_point(aes(color = Bone), shape = 21, position = position_jitterdodge(jitter.width = 1), 
             size = 3, stroke = .75, show.legend = F)+
  scale_fill_manual(values = c(col_p, col_s, col_m, col_i, col_ui))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black", "black"))+
  facet_grid(yaxis~category, switch = "y")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text.y = element_text(face = "bold", size = 11,
                                    color = "white"),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        strip.placement = "outside",
        legend.key.height = unit(1, "cm"),
        panel.grid.major.x = element_blank());p3




# Saving
p_site7_all <- ggarrange(p1, p2, p3, ncol = 1, nrow = 3, labels = c("A", "B", "C"));p_site7_all
ggsave("figures/png/supplementary_figure_3.png", device = "png", height = 10, width = 15, dpi = 300)
ggsave("figures/svg/supplementary_figure_3.svg", device = svg, height = 10, width = 15)




# Supplementary Figure 4--------------------------------------------------
# Microbial Analysis

# Comparison of Non-Oral and Oral relative abundances between bone types

# Reading data
biplot <- read.table("data/source_comparison", header = T, sep = "\t", check.names = F)
biplot$`Tissue Type` <- factor(biplot$`Tissue Type`, 
                               levels = c("Petrous", "Stapes", "Malleus", "Incus", "Unidentified Ossicle", "Tooth"))


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


library(ggbeeswarm)

# Plotting
pb_1 <- biplot %>%
  ggplot()+
  geom_boxplot(aes(x = `Tissue Type`, y = Oral, fill = `Tissue Type`, color = `Tissue Type`), 
               outliers = FALSE, show.legend = F)+
  geom_quasirandom(aes(x = `Tissue Type`, y = Oral, fill = `Tissue Type`),
                   shape = 21, size = 3, stroke = .75, width = 0.2)+
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
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75))+
  showSignificance(c(2,6), 7, -0.0001, t.s,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75))+
  showSignificance(c(3,6), 8, -0.0001, t.m,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75))+
  showSignificance(c(4,6), 9, -0.0001, t.i,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75))+
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
  guides(fill = guide_legend(nrow = 1));pb_1


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
                   shape = 21, size = 3, stroke = .75, width = 0.2)+
  scale_fill_manual(values = palette, labels = c("Petrous (n=101)", "Stapes (n=14)", "Malleus (n=13)",
                                                 "Incus (n=10)", "UI Ossicle (n=5)", "Tooth (n=5)"))+
  scale_color_manual(values = c("black","black","black", "black", "black",  "black"),
                     labels = c("Petrous (n=101)", "Stapes (n=14)", "Malleus (n=13)",
                                "Incus (n=10)", "UI Ossicle (n=5)", "Tooth (n=5)"))+
  ylim(0,34)+
  labs(x = "",
       y = "Non-Oral Reads (%)")+
  showSignificance(c(1,2), 22, -0.0001, p.s,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75))+
  showSignificance(c(1,4), 24, -0.0001, p.i,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75))+
  showSignificance(c(1,6), 26, -0.0001, t.p,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75))+
  showSignificance(c(2,6), 28, -0.0001, t.s,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75))+
  showSignificance(c(3,6), 30, -0.0001, t.m,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75))+
  showSignificance(c(4,6), 32, -0.0001, t.i,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75))+
  showSignificance(c(5,6), 34, -0.0001, t.u,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = .75))+
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
  guides(fill = guide_legend(nrow = 1));pb_2


combined <- ggarrange(pb_1, pb_2, ncol = 2, nrow = 1, labels = c("A", "B"), common.legend = T);combined

# Saving
ggsave("figures/svg/supplementary_figure_4.svg", device = svg, width = 15, height = 7)
ggsave("figures/png/supplementary_figure_4.png", device = "png", height = 7, width = 15, dpi = 600)
