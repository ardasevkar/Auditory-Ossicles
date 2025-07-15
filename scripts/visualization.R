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
library(tidyverse)
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



# Loading data ------------------------------------------------------------
all <- read.table("data/seqstats.txt", header = T, sep = "\t")
all$Bone <- factor(all$Bone, levels = c("Petrous", "Stapes", "Malleus", "Incus", "Unidentifiable Ossicle"))
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

hp$Bone[hp$Bone %in% "Unidentifiable Ossicle"] <- "Unidentifiable\nOssicle"
hp$comb[hp$comb %in% "Petrous|Unidentifiable Ossicle"] <- "Petrous|Unidentifiable\nOssicle"

hp$Bone <- factor(hp$Bone, levels = c("Petrous", "Stapes", "Malleus", "Incus", "Unidentifiable\nOssicle"))

hp$comb[hp$comb %in% "Petrous|Stapes"] <- 1
hp$comb[hp$comb %in% "Petrous|Malleus"] <- 2
hp$comb[hp$comb %in% "Petrous|Incus"] <- 3
hp$comb[hp$comb %in% "Petrous|Unidentifiable\nOssicle"] <- 4 

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

##### Second Panel ----------------------------------------------------------------
p_B <- hp %>%
  ggplot(aes(x = Bone, y = value))+
  ylim(0,85)+
  labs(x = "", y = "Human Proportion (%)")+
  geom_line(aes(group = Sample), color = "#121212", size = 1, show.legend = F)+
  geom_point(aes(fill = Bone, color = Bone), shape = 21, size = 5.5, show.legend = T, stroke = .75)+
  scale_fill_manual(values = c("#e63946", "#0a2342","#087e8a", "#ebf2fa", "#F58F29"))+
  scale_color_manual(values = c("black","grey70","black", "black", "black"))+
  showSignificance(c(1,2), 82, -0.0001, "*",
                   textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
                   panel = list(comb = 1))+
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
  filter(variable %in% c("Human Proportion (%)", "Downsampled Clonality (%)", 
                         "Avg Read Length (bp)", "PMD (%)", "Cont. Est."))

third$variable <- as.character(third$variable)
third$variable[third$variable %in% "Human Proportion (%)"] <- "Human Proportion"
third$variable[third$variable %in% "Downsampled Clonality (%)"] <- "Clonality"
third$variable[third$variable %in% "Avg Read Length (bp)"] <- "Mean Fragment Length"
third$variable[third$variable %in% "Cont. Est."] <- "Contamination"
third$variable[third$variable %in% "PMD (%)"] <- "PMD"


third$variable <- factor(third$variable, levels = c("Human Proportion", "Mean Fragment Length", 
                                                    "PMD", "Clonality", "Contamination"))
p_C <- third %>%
  filter(diff <= 500) %>%
  ggplot(aes(x = Bone, y = diff))+
  geom_hline(yintercept = 0, color = "#ff2800", linetype = "dashed", size = 1)+
  geom_boxplot(aes(fill = Bone))+
  labs(y = "Change Relative to Petrous (%)",
       x = "")+
  geom_point(aes(fill = Bone, color = Bone), shape = 21, size = 4)+
  facet_wrap(.~variable, scales = "free", nrow = 1)+
  scale_fill_manual(values = c("#0a2342","#087e8a", "#ebf2fa", "#F58F29"))+
  scale_color_manual(values = c("grey70","black", "black", "black"))+
  scale_x_discrete(labels = c("Stapes", "Malleus", "Incus","Unidentifiable\nOssicle"))+
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



##### Fourth Panel ------------------------------------------------------------
# Loading dataset
data <- read_xlsx("data/Supplementary Table 1.xlsx")
x <- match(x = c("Mean Read Length (bp)", 
                 "PMD 5'", "Contamination Estimate (contamMix)"), table = colnames(data))

colnames(data)[x] <- c("Avg Read Length (bp)", "PMD (%)", "Cont. Est.")

# Editting column values
data$Bone <- paste0(toupper(substr(data$Bone, 1, 1)), substr(data$Bone, 2, nchar(data$Bone)))
data$bone_other <- paste0(toupper(substr(data$`Bone Type`, 1, 1)), substr(data$`Bone Type`, 2, nchar(data$`Bone Type`)))
data$`Human Proportion (%)` <- data$`Human Proportion (%)`*100
data$`PMD (%)` <- data$`PMD (%)`*100

added_data <- data

added_data$Bone[added_data$Bone %in% "Unidentifiable Ossicle"] <- "Unidentifiable\nOssicle"

added_data$Bone <- factor(added_data$Bone, levels = c("Petrous", "Stapes", "Malleus", 
                                                      "Incus", "Unidentifiable\nOssicle"))

# Distribution of Human Proportion
hp <- added_data %>%
  select(Sample, Bone, `Human Proportion (%)`)

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

p.s <- paste0(sign(x$value[x$comb %in% "S.P"]));p.s
p.s <- paste0("***");p.s
p.u <- paste0(sign(x$value[x$comb %in% "U.P"]), " (p=", round(x$value[x$comb %in% "U.P"], 2), ")");p.u
p.u <- paste0(sign(x$value[x$comb %in% "U.P"]));p.u
s.m <- paste0(sign(x$value[x$comb %in% "M.S"]), " (p=", round(x$value[x$comb %in% "M.S"], 2), ")");s.m
s.m <- paste0(sign(x$value[x$comb %in% "M.S"]));s.m

p1 <- hp %>%
  ggplot(aes(x = Bone, y = `Human Proportion (%)`, fill = Bone, color = Bone))+
  labs(y = "Human Proportion (%)",
       x = "",
       title = "Change in sequence statistics after preparing library from other bones")+
  ylim(0,100)+
  geom_boxplot(size = .5, outlier.color = NA, alpha = .85, width = .5)+
  geom_point(aes(color = Bone), shape = 21, position = position_jitterdodge(jitter.width = 1), 
             size = 3, stroke = .75, show.legend = F)+
  scale_fill_manual(values = c("#e63946", "#0a2342","#087e8a", "#ebf2fa", "#F58F29"))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black"))+
  showSignificance(c(1,2), 80, -0.0001, p.s,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = 1))+
  showSignificance(c(2,3), 85, -0.0001, s.m,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = 1))+
  showSignificance(c(1,5), 95, -0.0001, p.u,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = 1))+
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
  select(Sample, Bone, `Avg Read Length (bp)`)

avg %>%
  ggplot(aes(x = Bone, y = `Avg Read Length (bp)`))+
  ylim(0, 100)+
  geom_boxplot()+
  geom_point()

# Running pairwise wilcoxon tests on bone types
x <-  as.data.frame(pairwise.wilcox.test(avg$`Avg Read Length (bp)`, avg$Bone, 
                                         p.adjust.method =  "BH", paired = F)$p.value)

x$comb <- row.names(x)

x <- melt(x, id.vars = "comb")
x <- x %>%
  filter(!is.na(value)) %>%
  filter(value <= 0.05) %>%
  mutate(comb = paste0(substr(comb, 1, 1), ".", substr(variable, 1, 1)))

p.s <- paste0(sign(x$value[x$comb %in% "S.P"]), " (p=", round(x$value[x$comb %in% "S.P"], 3), ")");p.s
p.s <- paste0(sign(x$value[x$comb %in% "S.P"]));p.s

p2 <- avg %>%
  ggplot(aes(x = Bone, y = `Avg Read Length (bp)`, fill = Bone, color = Bone))+
  labs(y = "Mean Fragment Length (bp)",
       x = "",
       title = "Change in sequence statistics after preparing library from other bones")+
  ylim(40,100)+
  geom_boxplot(size = .5, outlier.color = NA, alpha = .85, width = .5)+
  geom_point(aes(color = Bone), shape = 21, position = position_jitterdodge(jitter.width = 1), 
             size = 3, stroke = .75, show.legend = F)+
  scale_fill_manual(values = c("#e63946", "#0a2342","#087e8a", "#ebf2fa", "#F58F29"))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black"))+
  showSignificance(c(1,2), 92, -0.0001, p.s,
                   textParams = list(size = 4, fontface = "bold"), segmentParams = list(size = 1))+
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


# Distribution of post-mortem damage levels
pmd <- added_data %>%
  select(Sample, Bone, `PMD (%)`)

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
  geom_boxplot(size = .5, outlier.color = NA, alpha = .85, width = .5)+
  geom_point(aes(color = Bone), shape = 21, position = position_jitterdodge(jitter.width = 1), 
             size = 3, stroke = .75, show.legend = F)+
  scale_fill_manual(values = c("#e63946", "#0a2342","#087e8a", "#ebf2fa", "#F58F29"))+
  scale_color_manual(values = c("black","grey70", "black", "black", "black"))+
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
                       labels = c("A","", "B"), align = "v");top_panel

bottom_panel <- ggarrange(p_C, p_D, ncol = 1, nrow = 2, labels = c("C", "D"));bottom_panel

all_plot <- ggarrange(top_panel, NULL, bottom_panel, nrow = 3, ncol = 1, heights = c(.45, 0.02 ,1));all_plot

ggsave("figures/png/main_plot.png", device = png, width = 21, height = 15, dpi = 300)
ggsave("figures/svg/main_plot.svg", device = svg, width = 21, height = 15)







# Supplementary Figure 1 --------------------------------------------------
# The code below contains scripts for Supplementary Figure 1 of the manuscript, which summarizes the distribution of samples that have both ossicle and petrous libraries.

data <- read_xlsx("data/Supplementary Table 1.xlsx")
data <- data %>%
  filter(`Library Type` %in% "ossicle&petrous")

# Calculating sample size for each location
map_data <- data %>%
  select(Sample, Site, Latitude, Longitude) %>%
  group_by(Site) %>%
  distinct(Sample, .keep_all = T) %>%
  mutate(n_sample = n(),
         label = paste0(Site, " (n=", n_sample, ")")) %>%
  distinct(Site, .keep_all = T)


# To download map data in raster format, use the link below:
# https://github.com/nvkelso/natural-earth-raster/tree/master/10m_rasters/NE1_HR_LC_SR_W_DR

# Loading map data
raster50 <- ne_load(
  destdir = "../data/",
  scale = 10,
  category = "raster",
  type = "NE1_HR_LC_SR_W",
  returnclass = "sf")

# Cropping map data
map_extent <- raster::extent(25, 46, 35, 43)
nat.crop <- raster::stack(raster::crop(raster50, y = map_extent))

rast.table <- data.frame(xyFromCell(nat.crop, 1:ncell(nat.crop)),
                         getValues(nat.crop/255))


# Assign rgb colors according to elevation
rast.table$rgb <- with(rast.table, rgb(NE1_HR_LC_SR_W_1,
                                       NE1_HR_LC_SR_W_2,
                                       NE1_HR_LC_SR_W_3,
                                       1))

# Converting colors on black and white scale
rast.table$rgb <- ColToGray(rast.table$rgb)

ggplot() +
  geom_raster(data = rast.table, aes(x = x, y = y),fill = rast.table$rgb) +
  coord_fixed(ratio = 1.3) +
  geom_point(data = map_data, aes(x = Longitude, y = Latitude), shape = 25, fill = "#121212", size = 3, stroke = .1)+
  geom_label_repel(data = map_data, aes(x = Longitude, y = Latitude, label = label),
                   direction = "y", nudge_y = .1, fontface = "bold", size = 4, fill = alpha("white", .75))+
  theme_minimal()+
  theme(legend.position=c(.93,.15),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(color = "black",size = 15.5, face = "bold"),
        #legend.key = element_rect(fill = "grey60"),
        legend.key.size = unit(.8, "cm"),
        legend.box.background = element_rect(color = NA, fill = alpha("white", 0.8), size = 10),
        legend.title = element_text(color = "black", size = 16, face = "bold"))

ggsave("figures/png/sample_map.png", device = png, dpi = 300, height = 12, width = 14)
ggsave("figures/svg/sample_map.svg", device = svg, height = 12, width = 14)





# Supplementary Figure 2 --------------------------------------------------
# Code below contains summary plots of sequence statistics of ossicle and petrous libraries prepared from the same individual.

# Loading data
seqstats <- read.table("data/seqstats.txt", header = T, sep = "\t")
seqstats$Bone[seqstats$Bone %in% "Unidentifiable Ossicle"] <- "Unidentifiable\nOssicle"
seqstats$Bone <- factor(seqstats$Bone, levels = c("Petrous", "Stapes", "Malleus", "Incus", "Unidentifiable\nOssicle"))
seqstats$comb[seqstats$comb %in% "Petrous|Unidentifiable Ossicle"] <- "Petrous|Unidentifiable\nOssicle"

##### Distribution of Human Proportion ----------------------------------------
hp <- seqstats %>%
  filter(variable %in% "Human Proportion (%)") %>%
  as.data.frame()

hp$Bone <- factor(hp$Bone, levels = c("Petrous", "Stapes", "Malleus", "Incus", "Unidentifiable\nOssicle"))

hp$comb[hp$comb %in% "Petrous|Stapes"] <- 1
hp$comb[hp$comb %in% "Petrous|Malleus"] <- 2
hp$comb[hp$comb %in% "Petrous|Incus"] <- 3
hp$comb[hp$comb %in% "Petrous|Unidentifiable\nOssicle"] <- 4

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
  scale_fill_manual(values = c("#e63946", "#0a2342","#087e8a", "#ebf2fa", "#F58F29"))+
  scale_color_manual(values = c("black","grey70","black", "black", "black"))+
  showSignificance(c(1,2), 85, -0.0001, star_1,
                   textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
                   panel = list(comb = 1))+
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

avg$Bone <- factor(avg$Bone, levels = c("Petrous", "Stapes", "Malleus", "Incus", "Unidentifiable\nOssicle"))

avg$comb[avg$comb %in% "Petrous|Stapes"] <- 1
avg$comb[avg$comb %in% "Petrous|Malleus"] <- 2
avg$comb[avg$comb %in% "Petrous|Incus"] <- 3
avg$comb[avg$comb %in% "Petrous|Unidentifiable\nOssicle"] <- 4

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

p2 <- avg %>%
  ggplot(aes(x = Bone, y = value))+
  ylim(50,95)+
  labs(x = "", y = "Mean Fragment Length (bp)")+
  geom_line(aes(group = Sample), color = "#121212", size = 1, show.legend = F)+
  geom_point(aes(fill = Bone, color = Bone), shape = 21, size = 5.5, show.legend = T, stroke = .75)+
  scale_fill_manual(values = c("#e63946", "#0a2342","#087e8a", "#ebf2fa", "#F58F29"))+
  scale_color_manual(values = c("black","grey70","black", "black", "black"))+
  showSignificance(c(1,2), 92, -0.0001, star_1,
                   textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1),
                   panel = list(comb = 1))+
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

pmd$Bone <- factor(pmd$Bone, levels = c("Petrous", "Stapes", "Malleus", "Incus", "Unidentifiable\nOssicle"))

pmd$comb[pmd$comb %in% "Petrous|Stapes"] <- 1
pmd$comb[pmd$comb %in% "Petrous|Malleus"] <- 2
pmd$comb[pmd$comb %in% "Petrous|Incus"] <- 3
pmd$comb[pmd$comb %in% "Petrous|Unidentifiable\nOssicle"] <- 4

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


p3 <- pmd %>%
  ggplot(aes(x = Bone, y = value))+
  ylim(15,55)+
  labs(x = "", y = "PMD (%)")+
  geom_line(aes(group = Sample), color = "#121212", size = 1, show.legend = F)+
  geom_point(aes(fill = Bone, color = Bone), shape = 21, size = 5.5, show.legend = T, stroke = .75)+
  scale_fill_manual(values = c("#e63946", "#0a2342","#087e8a", "#ebf2fa", "#F58F29"))+
  scale_color_manual(values = c("black","grey70","black", "black", "black"))+
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
  filter(variable %in% "Downsampled Clonality (%)") %>%
  as.data.frame()

clon$Bone <- factor(clon$Bone, levels = c("Petrous", "Stapes", "Malleus", "Incus", "Unidentifiable\nOssicle"))

clon$comb[clon$comb %in% "Petrous|Stapes"] <- 1
clon$comb[clon$comb %in% "Petrous|Malleus"] <- 2
clon$comb[clon$comb %in% "Petrous|Incus"] <- 3
clon$comb[clon$comb %in% "Petrous|Unidentifiable\nOssicle"] <- 4

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

p4 <- clon %>%
  ggplot(aes(x = Bone, y = value))+
  # ylim(0,23)+
  labs(x = "", y = "Downsampled Clonality (%)")+
  geom_line(aes(group = Sample), color = "#121212", size = 1, show.legend = F)+
  geom_point(aes(fill = Bone, color = Bone), shape = 21, size = 5.5, show.legend = T, stroke = .75)+
  scale_fill_manual(values = c("#e63946", "#0a2342","#087e8a", "#ebf2fa", "#F58F29"))+
  scale_color_manual(values = c("black","grey70","black", "black", "black"))+
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

ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "C"), ncol = 1, nrow = 4)
ggsave("figures/png/seqstats_same.png", device = png, width = 10, height = 15, dpi = 300)
ggsave("figures/svg/seqstats_same.svg", device = png, width = 10, height = 15)







