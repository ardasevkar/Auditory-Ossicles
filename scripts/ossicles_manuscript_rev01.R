# This script contains all codes for Sağlıcan & Sevkar et. al, 2025 "A comparison of ancient DNA yields across ossicles and the petrous bone reveals the highest preservation in the stapes and incus" manuscript.

setwd("../")

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
library(ggplot2)
library(dplyr)
library(ggtern)
library(ggpubr)
library(ggforce)
library(rnaturalearth)
library(rnaturalearthdata)
library(superb)
library(readxl)
library(reshape2)
library(tidyverse)
library(vegan)

# Color palette
palette <- c("#e63845", "#07213f","#048BA8", "#FAF3DD", "#E0913E", "#F194B4")
#palette <- c("#e63946", "#130870","#2b35af", "#4360ed", "#4895ef", "#faf0ca")

col_p <- palette[1]
col_s <- palette[2]
col_m <- palette[3]
col_i <- palette[4]
col_ui <- palette[5]
col_t <- palette[6]

# Loading Dataset
raw <- read_xlsx("data/Supplementary_Tables_rev01.xlsx", na = "S1")

# Statistics on Petrous and Ossicle Libraries Prepared From the Same Indivudals-------------------------

#Filtering dataset
data <- raw %>%
  filter(`Library Type` %nin% c("only ossicle", "only petrous"))

# Changing column names
x <- match(x = c("Mean Read Length (bp)", 
                 "PMD 5'", "Contamination Estimate (contamMix)"), table = colnames(data))

colnames(data)[x] <- c("Avg Read Length (bp)", "PMD (%)", "Cont. Est.")
colnames(data)[which(colnames(data) %in% "Individual")] <- "Sample"

# Editting column values
data$Bone <- paste0(toupper(substr(data$Bone, 1, 1)), substr(data$Bone, 2, nchar(data$Bone)))
data$bone_other <- paste0(toupper(substr(data$`Bone Type`, 1, 1)), substr(data$`Bone Type`, 2, nchar(data$`Bone Type`)))
data$`Human Proportion (%)` <- data$`Human Proportion (%)`*100
data$`PMD (%)` <- data$`PMD (%)`*100




###### Petrous vs Stapes -------------------------------------------------------
test <- data %>%
  filter(Bone %in% c("Stapes", "Petrous")) %>%
  group_by(Sample)%>%
  mutate(n=n())%>%
  filter(n>1)%>%
  select(Sample, Bone, `Library ID`, `Human Proportion (%)`, `Clonality (%)`, 
         `Avg Read Length (bp)`, `PMD (%)`, `Cont. Est.`)%>%
  ungroup()

test_m <- melt(test, id.vars = c("Sample", "Bone", "Library ID"))
test_m$value <- as.numeric(test_m$value)

p1 <- test_m %>%
  filter(variable %nin% c( "Cont. Est.", "Avg Read Length (bp)")) %>%
  ggplot(aes(x = Sample, y = value, fill = Bone))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_grid(variable~., scales = "free", switch = "y")+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c(col_p, col_s))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 10,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        panel.grid.major.x = element_blank());p1


p2_1 <-  test_m %>%
  filter(variable %in% c("Avg Read Length (bp)")) %>%
  ggplot(aes(x = Sample, y = value, fill = Bone))+
  geom_point(shape = 21, size = 5, stroke = .25, position = position_dodge(width = .2))+
  facet_grid(variable~., scales = "free", switch = "y")+
  ylim(40,95)+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c(col_p, col_s))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 10,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside");p2_1


contamMix <- test_m %>%
  filter(variable %in% c("Cont. Est."))

p2_2 <- contamMix %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = Sample, y = value, fill = Bone))+
  geom_point(shape = 21, size = 5, stroke = .25, position = position_dodge(width = .2))+
  facet_grid(variable~., scales = "free", switch = "y")+
  ylim(.7,1)+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c(col_p, col_s))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 10,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside");p2_2


p2 <- ggarrange(p2_1, p2_2, nrow = 2, ncol = 1, align = "hv", heights = c(1,1), labels = c("B", "C"));p2
ggarrange(p1, p2, labels = c("A"), nrow = 1)

ggsave(filename = "figures/png/other/p_vs_s_1.png", height = 8, width = 15, device = "png", dpi = 300)
ggsave(filename = "figures/svg/other/p_vs_s_1.svg", height = 8, width = 15, device = svg)


# Calculating changing percentages
change <- test_m %>%
  filter(variable %in% "Human Proportion (%)")

temp <- change %>%
  select(Sample, Bone, value) %>%
  spread(Bone, value)

# Running Wilcoxon signed-rank test
test <- wilcox.test(temp$Petrous, temp$Stapes, paired = T)$p.value
#star <- paste0(sign(test), " (p=", round(test, digits = 3), ")")
star <- sign(test)

# Plotting first panel
p1 <- change %>%
  ggplot(aes(x = Bone, y = value))+
  ylim(0,90)+
  labs(y = "Human Proportion (%)",
       x = "")+
  geom_line(aes(group = Sample, color = Sample), size = 2, show.legend = F)+
  geom_point(aes(fill = Sample), shape = 21, size = 5, show.legend = T, stroke = .75)+
  scale_color_viridis_d(option = "rocket")+
  scale_fill_viridis_d(option = "rocket")+
  showSignificance(c(1,2), 85, -0.1, star,
                   textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 12, face = "bold"),
        plot.caption = element_text(face = "italic"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside")+
  guides(color = guide_legend(ncol = 7, size = 5),
         fill = guide_legend(ncol = 7, size = 5));p1


stapes <- test_m %>%
  group_by(Sample, variable) %>%
  mutate(diff = ((value[Bone %in% "Stapes"]-value[Bone %in% "Petrous"])/value[Bone %in% "Petrous"])*100,
         type = "Stapes",
         comb = "Petrous|Stapes") %>%
  filter(variable %in% c("Human Proportion (%)", "mt Proportion (%)", "Avg Read Length (bp)", 
                         "PMD (%)", "Cont. Est.", "mt_coverage", "Clonality (%)")) %>%
  select(type, Bone, Sample, value, diff, variable, comb) %>%
  distinct(Sample, Bone, .keep_all = T)

# Drawing second panel
p2_1 <- stapes %>%
  filter(Bone %in% "Stapes") %>%
  filter(variable %in% "Human Proportion (%)") %>%
  ggplot(aes(x = variable, y = diff))+
  facet_zoom(y = diff <= 350, ylim = c(0,310), zoom.size = 1.5, horizontal = F)+
  geom_boxplot(lwd = .50, fill = "#fbfbf9", width = .5)+
  geom_hline(yintercept = 0, linetype = "dashed", size = .75, color = "red")+
  labs(x = "",
       y = "Change (%)",
       title = "Change in sequence statistics after preparing library from stapes")+
  geom_point(aes(fill = Sample), shape = 21, size = 5, show.legend = T, stroke = .75)+
  scale_fill_viridis_d(option = "rocket")+
  scale_x_discrete(labels = "Human Proportion")+
  theme_bw()+
  theme(
    zoom.y = element_rect(fill = "grey90", color = "grey70"),
    axis.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold"),
    legend.background = element_rect(fill = "white"))+
  guides(color = guide_legend(ncol = 7, size = 5),
         fill = guide_legend(ncol = 7, size = 5));p2_1

summary(stapes$diff[stapes$variable %in% "Human Proportion (%)"])
summary(stapes$diff[stapes$variable %in% "Avg Read Length (bp)"])
summary(stapes$diff[stapes$variable %in% "PMD (%)"])
summary(stapes$diff[stapes$variable %in% "Cont. Est."])



x <- stapes %>%
  filter(variable %in% "Avg Read Length (bp)") %>%
  select(Sample, value, Bone) %>%
  arrange(Sample)

mean(x$value[x$Bone %in% "Petrous"]) - mean(x$value[x$Bone %in% "Stapes"])
wilcox.test(x$value[x$Bone %in% "Petrous"], x$value[x$Bone %in% "Stapes"], paired = T)



#Drawing Third Panel
p2_2 <- stapes %>%
  filter(!is.na(diff)) %>%
  filter(Bone %in% "Stapes") %>%
  filter(variable %nin% c("Human Proportion (%)", "mt_coverage", "mt Proportion (%)")) %>%
  ggplot(aes(x = variable, y = diff))+
  geom_boxplot(lwd = .50, fill = "#fbfbf9", width = .75)+
  geom_hline(yintercept = 0, linetype = "dashed", size = .75, color = "red")+
  scale_x_discrete(labels = c("Clonality", "Avg Read Length", "PMD", "Cont. Est."))+
  #geom_line(aes(group = Sample, color = Sample), size = 2)+
  labs(x = "",
       y = "Change (%)",
       title = "Change in sequence statistics after preparing library from stapes")+
  geom_point(aes(fill = Sample), shape = 21, size = 5, show.legend = T, stroke = .75)+
  scale_fill_viridis_d(option = "rocket")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        legend.background = element_rect(fill = "white"))+
  guides(color = guide_legend(ncol = 7, size = 5),
         fill = guide_legend(ncol = 7, size = 5));p2_2

p2 <- ggarrange(p2_1, p2_2, ncol = 2, nrow = 1, legend = F, labels = c("B", "C"));p2

ggarrange(p1, p2_1, p2_2, ncol = 3, nrow = 1, common.legend = T, labels = c("A", "B", "C"), legend = F)

ggsave(filename = "figures/png/other/p_vs_s_2.png", height = 5, width = 15, device = "png", dpi = 300)
ggsave(filename = "figures/svg/other/p_vs_s_2.svg", height = 5, width = 15, device = svg)


###### Petrous vs Malleus ------------------------------------------------------
test <- data %>%
  filter(Bone %in% c("Petrous", "Malleus")) %>%
  group_by(Sample)%>%
  mutate(n=n())%>%
  filter(n>1)%>%
  select(Sample, Bone, `Library ID`, `Human Proportion (%)`, `Clonality (%)`, 
         `Avg Read Length (bp)`, `PMD (%)`, `Cont. Est.`)%>%
  ungroup()

test_m <- melt(test, id.vars = c("Sample", "Bone", "Library ID"))
test_m$value <- as.numeric(test_m$value)

test_m$Bone <- factor(test$Bone, levels = c("Petrous", "Malleus"))

p1 <- test_m %>%
  filter(variable %nin% c( "Cont. Est.", "Avg Read Length (bp)")) %>%
  ggplot(aes(x = Sample, y = value, fill = Bone))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_grid(variable~., scales = "free", switch = "y")+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c(col_p, col_m))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 10,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        panel.grid.major.x = element_blank());p1


p2_1 <-  test_m %>%
  filter(variable %in% c("Avg Read Length (bp)")) %>%
  ggplot(aes(x = Sample, y = value, fill = Bone))+
  geom_point(shape = 21, size = 5, stroke = .25, position = position_dodge(width = .2))+
  facet_grid(variable~., scales = "free", switch = "y")+
  ylim(40,95)+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c(col_p, col_m))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 10,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside");p2_1

contamMix <- test_m %>%
  filter(variable %in% c("Cont. Est."))

p2_2 <- contamMix %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = Sample, y = value, fill = Bone))+
  geom_point(shape = 21, size = 5, stroke = .25, position = position_dodge(width = .2))+
  facet_grid(variable~., scales = "free", switch = "y")+
  ylim(.80,1)+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c(col_p, col_m))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 10,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside");p2_2


p2 <- ggarrange(p2_1, p2_2, nrow = 2, ncol = 1, align = "hv", heights = c(1,1), labels = c("B", "C"));p2
ggarrange(p1, p2, labels = c("A"), nrow = 1)

ggsave(filename = "figures/png/other/p_vs_m_1.png", height = 8, width = 15, device = "png", dpi = 300)
ggsave(filename = "figures/svg/other/p_vs_m_1.svg", height = 8, width = 15, device = svg)


# Calculating changing percentages
change <- test_m %>%
  filter(variable %in% "Human Proportion (%)")

temp <- change %>%
  select(Sample, Bone, value) %>%
  spread(Bone, value)

# Running Wilcoxon signed-rank test
test <- wilcox.test(temp$Petrous, temp$Malleus, paired = T)$p.value
star <- paste0(sign(test), " (p=", round(test, digits = 3), ")")

# Plotting first panel
p1 <- change %>%
  ggplot(aes(x = Bone, y = value))+
  ylim(0,90)+
  labs(y = "Human Proportion (%)",
       x = "")+
  geom_line(aes(group = Sample, color = Sample), size = 2, show.legend = F)+
  geom_point(aes(fill = Sample), shape = 21, size = 5, show.legend = T, stroke = .75)+
  scale_color_viridis_d(option = "rocket")+
  scale_fill_viridis_d(option = "rocket")+
  showSignificance(c(1,2), 85, -0.1, star,
                   textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 12, face = "bold"),
        plot.caption = element_text(face = "italic"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside")+
  guides(color = guide_legend(ncol = 7, size = 5),
         fill = guide_legend(ncol = 7, size = 5));p1


malleus <- test_m %>%
  group_by(Sample, variable) %>%
  mutate(diff = ((value[Bone %in% "Malleus"]-value[Bone %in% "Petrous"])/value[Bone %in% "Malleus"])*100,
         type = "Malleus",
         comb = "Petrous|Malleus") %>%
  filter(variable %in% c("Human Proportion (%)", "mt Proportion (%)", "Avg Read Length (bp)", 
                         "PMD (%)", "Cont. Est.", "mt_coverage", "Clonality (%)")) %>%
  select(type, Bone, Sample, value, diff, variable, comb) %>%
  distinct(Sample, Bone, .keep_all = T)

# Drawing second panel
p2_1 <- malleus %>%
  filter(Bone %in% "Malleus") %>%
  filter(variable %in% "Human Proportion (%)") %>%
  ggplot(aes(x = variable, y = diff))+
  #facet_zoom(y = diff <= 350, ylim = c(0,310), zoom.size = 1.5, horizontal = F)+
  geom_boxplot(lwd = .50, fill = "#fbfbf9", width = .5)+
  geom_hline(yintercept = 0, linetype = "dashed", size = .75, color = "red")+
  labs(x = "",
       y = "Change (%)",
       title = "Change in sequence statistics after preparing library from stapes")+
  geom_point(aes(fill = Sample), shape = 21, size = 5, show.legend = T, stroke = .75)+
  scale_fill_viridis_d(option = "rocket")+
  scale_x_discrete(labels = "Human Proportion")+
  theme_bw()+
  theme(
    zoom.y = element_rect(fill = "grey90", color = "grey70"),
    axis.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold"),
    legend.background = element_rect(fill = "white"))+
  guides(color = guide_legend(ncol = 7, size = 5),
         fill = guide_legend(ncol = 7, size = 5));p2_1

summary(malleus$diff[malleus$variable %in% "Human Proportion (%)"])
summary(malleus$diff[malleus$variable %in% "Avg Read Length (bp)"])
summary(malleus$diff[malleus$variable %in% "PMD (%)"])
summary(malleus$diff[malleus$variable %in% "Cont. Est."])

x <- malleus %>%
  filter(variable %in% "Avg Read Length (bp)") %>%
  select(Sample, value, Bone) %>%
  arrange(Sample)

mean(x$value[x$Bone %in% "Petrous"]) - mean(x$value[x$Bone %in% "Malleus"])
wilcox.test(x$value[x$Bone %in% "Petrous"], x$value[x$Bone %in% "Malleus"], paired = T)

#Drawing Third Panel
p2_2 <- malleus %>%
  filter(Sample %nin% "Ind004") %>%
  filter(!is.na(diff)) %>%
  filter(Bone %in% "Malleus") %>%
  filter(variable %nin% c("Human Proportion (%)", "mt_coverage", "mt Proportion (%)")) %>%
  ggplot(aes(x = variable, y = diff))+
  geom_boxplot(lwd = .50, fill = "#fbfbf9", width = .75)+
  geom_hline(yintercept = 0, linetype = "dashed", size = .75, color = "red")+
  scale_x_discrete(labels = c("Clonality", "Avg Read Length", "PMD", "Cont. Est."))+
  #geom_line(aes(group = Sample, color = Sample), size = 2)+
  labs(x = "",
       y = "Change (%)",
       title = "Change in sequence statistics after preparing library from stapes")+
  geom_point(aes(fill = Sample), shape = 21, size = 5, show.legend = T, stroke = .75)+
  scale_fill_viridis_d(option = "rocket")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        legend.background = element_rect(fill = "white"))+
  guides(color = guide_legend(ncol = 7, size = 5),
         fill = guide_legend(ncol = 7, size = 5));p2_2

p2 <- ggarrange(p2_1, p2_2, ncol = 2, nrow = 1, legend = F, labels = c("B", "C"));p2

ggarrange(p1, p2_1, p2_2, ncol = 3, nrow = 1, common.legend = T, labels = c("A", "B", "C"), legend = F)

ggsave(filename = "figures/png/other/p_vs_m_2.png", height = 5, width = 15, device = "png", dpi = 300)
ggsave(filename = "figures/svg/other/p_vs_m_2.svg", height = 5, width = 15, device = svg)




###### Petrous vs Incus --------------------------------------------------------
test <- data %>%
  filter(Bone %in% c("Incus", "Petrous")) %>%
  group_by(Sample)%>%
  mutate(n=n())%>%
  filter(n>1)%>%
  select(Sample, Bone, `Library ID`, `Human Proportion (%)`, `Clonality (%)`, 
         `Avg Read Length (bp)`, `PMD (%)`, `Cont. Est.`)%>%
  ungroup()

test_m <- melt(test, id.vars = c("Sample", "Bone", "Library ID"))
test_m$Bone <- factor(test_m$Bone, levels = c("Petrous", "Incus"))
test_m$value <- as.numeric(test_m$value)

p1 <- test_m %>%
  filter(variable %nin% c("mt_coverage", "Cont. Est.", "Avg Read Length (bp)")) %>%
  ggplot(aes(x = Sample, y = value, fill = Bone))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_grid(variable~., scales = "free", switch = "y")+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c(col_p, col_i))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 10,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        panel.grid.major.x = element_blank());p1

p2_1 <-  test_m %>%
  filter(variable %in% c("Avg Read Length (bp)")) %>%
  ggplot(aes(x = Sample, y = value, fill = Bone))+
  geom_point(shape = 21, size = 5, stroke = .25, position = position_dodge(width = .2))+
  facet_grid(variable~., scales = "free", switch = "y")+
  ylim(40,80)+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c(col_p, col_i))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 10,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside");p2_1


contamMix <- test_m %>%
  filter(variable %in% c("Cont. Est."))

p2_2 <- contamMix %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = Sample, y = value, fill = Bone))+
  geom_point(shape = 21, size = 5, stroke = .25, position = position_dodge(width = .2))+
  facet_grid(variable~., scales = "free", switch = "y")+
  ylim(.6,1)+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c(col_p, col_i))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 10,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside");p2_2


p2 <- ggarrange(p2_1, p2_2, nrow = 2, ncol = 1, align = "hv", heights = c(1,1), labels = c("B", "C"));p2
ggarrange(p1, p2, labels = c("A"), nrow = 1)

ggsave(filename = "figures/png/other/p_vs_i_1.png", height = 6, width = 15, device = "png", dpi = 300)
ggsave(filename = "figures/svg/other/p_vs_i_1.svg", height = 6, width = 15, device = svg)


# Calculating changing percentages
change <- test_m %>%
  filter(variable %in% "Human Proportion (%)")

temp <- change %>%
  select(Sample, Bone, value) %>%
  spread(Bone, value)

test <- wilcox.test(temp$Petrous, temp$Incus, paired = T)$p.value

star <- paste0(sign(test), " (p=", round(test, digits = 3), ")")

p1 <- change %>%
  ggplot(aes(x = Bone, y = value))+
  ylim(0,90)+
  labs(y = "Human Proportion (%)",
       x = "")+
  geom_line(aes(group = Sample, color = Sample), size = 2, show.legend = F)+
  geom_point(aes(fill = Sample), shape = 21, size = 5, show.legend = T, stroke = .75)+
  scale_color_viridis_d(option = "rocket")+
  scale_fill_viridis_d(option = "rocket")+
  showSignificance(c(1,2), 85, -0.1, star,
                   textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 12, face = "bold"),
        plot.caption = element_text(face = "italic"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside");p1


incus <- test_m %>%
  group_by(Sample, variable) %>%
  mutate(diff = ((value[Bone %in% "Incus"]-value[Bone %in% "Petrous"])/value[Bone %in% "Petrous"])*100,
         type = "Incus",
         comb = "Petrous|Incus") %>%
  filter(variable %in% c("Human Proportion (%)", "mt Proportion (%)",
                         "Avg Read Length (bp)", "PMD (%)", "Cont. Est.", "mt_coverage", "Clonality (%)")) %>%
  select(type, Bone, Sample, value, diff, variable, comb) %>%
  distinct(Sample, Bone, .keep_all = T)

# Drawing Second Panel
p2_1 <- incus %>%
  filter(Bone %in% "Incus") %>%
  filter(variable %in% "Human Proportion (%)") %>%
  ggplot(aes(x = variable, y = diff))+
  # ylim(0,80)+
  #facet_zoom(y = diff <= 350, ylim = c(0,310), zoom.size = 1.5, horizontal = F)+
  geom_boxplot(lwd = .50, fill = "#fbfbf9", width = .5)+
  geom_hline(yintercept = 0, linetype = "dashed", size = .75, color = "red")+
  labs(x = "",
       y = "Change (%)",
       title = "Change in sequence statistics after preparing library from incus")+
  geom_point(aes(fill = Sample), shape = 21, size = 5, show.legend = T, stroke = .75)+
  scale_fill_viridis_d(option = "rocket")+
  scale_x_discrete(labels = "Human Proportion")+
  theme_bw()+
  theme(
    zoom.y = element_rect(fill = "grey90", color = "grey70"),
    axis.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold"),
    legend.background = element_rect(fill = "white"));p2_1


# Drawing Third Panel
p2_2 <- incus %>%
  filter(Bone %in% "Incus") %>%
  filter(variable %nin% c("Human Proportion (%)", "mt_coverage", "mt Proportion (%)")) %>%
  ggplot(aes(x = variable, y = diff))+
  geom_hline(yintercept = 0, linetype = "dashed", size = .75, color = "red")+
  #ylim(-10,20)+
  geom_boxplot(lwd = .50, fill = "#fbfbf9", width = .75)+
  #geom_line(aes(group = Sample, color = Sample), size = 2)+
  labs(x = "",
       y = "Change (%)",
       title = "Change in sequence statistics after preparing library from stapes")+
  geom_point(aes(fill = Sample), shape = 21, size = 5, show.legend = T, stroke = .75)+
  scale_fill_viridis_d(option = "rocket")+
  scale_x_discrete(labels = c("Clonality", "Avg Read Length", "PMD", "Cont. Est."))+
  #facet_grid(.~variable, scales = "free", space = "free")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        legend.background = element_rect(fill = "white"));p2_2


p2 <- ggarrange(p2_1, p2_2, ncol = 2, nrow = 1, legend = F, labels = c("B", "C"));p2

ggarrange(p1, p2_1, p2_2, ncol = 3, nrow = 1, common.legend = T, labels = c("A", "B", "C"), legend = F)

ggsave(filename = "figures/png/other/p_vs_i_2.png", height = 5, width = 15, device = "png", dpi = 300)
ggsave(filename = "figures/svg/other/p_vs_i_2.svg", height = 5, width = 15, device = svg)




###### Petrous vs Unidentifed Ossicle --------------------------------------------------
test <- data %>%
  filter(Bone %in% c("Unidentified Ossicle", "Petrous")) %>%
  group_by(Sample)%>%
  mutate(n=n())%>%
  filter(n>1)%>%
  select(Sample, Bone, `Library ID`, `Human Proportion (%)`, `Clonality (%)`, 
         `Avg Read Length (bp)`, `PMD (%)`, `Cont. Est.`)%>%
  ungroup()

test_m <- melt(test, id.vars = c("Sample", "Bone", "Library ID"))
test_m$Bone <- factor(test_m$Bone, levels = c("Petrous", "Unidentified Ossicle"))
test_m$value <- as.numeric(test_m$value)

p1 <- test_m %>%
  filter(variable %nin% c("mt_coverage", "Cont. Est.", "Avg Read Length (bp)")) %>%
  ggplot(aes(x = Sample, y = value, fill = Bone))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_grid(variable~., scales = "free", switch = "y")+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c(col_p, col_ui))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 10,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        panel.grid.major.x = element_blank());p1

p2_1 <-  test_m %>%
  filter(variable %in% c("Avg Read Length (bp)")) %>%
  ggplot(aes(x = Sample, y = value, fill = Bone))+
  geom_point(shape = 21, size = 5, stroke = .25, position = position_dodge(width = .2))+
  facet_grid(variable~., scales = "free", switch = "y")+
  ylim(40,80)+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c(col_p, col_ui))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 10,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside");p2_1

contamMix <- test_m %>%
  filter(variable %in% c("Cont. Est."))

p2_2 <- contamMix %>%
  #filter(!is.na(value)) %>%
  ggplot(aes(x = Sample, y = value, fill = Bone))+
  geom_point(shape = 21, size = 5, stroke = .25, position = position_dodge(width = .2))+
  facet_grid(variable~., scales = "free", switch = "y")+
  ylim(.5,1)+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c(col_p, col_ui))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 10,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside");p2_2

p2 <- ggarrange(p2_1, p2_2, nrow = 2, ncol = 1, align = "hv", heights = c(1,1), labels = c("B", "C"));p2
ggarrange(p1, p2, labels = c("A"), nrow = 1)

ggsave(filename = "figures/png/other/p_vs_un_1.png", height = 6, width = 11, device = "png", dpi = 300)
ggsave(filename = "figures/svg/other/p_vs_un_1.svg", height = 6, width = 11, device = svg)

# Testing changing percentage
change <- test_m %>%
  filter(variable %in% "Human Proportion (%)")

temp <- change %>%
  select(Sample, Bone, value) %>%
  spread(Bone, value)

test <- wilcox.test(temp$Petrous, temp$`Unidentified Ossicle`, paired = T)$p.value

star <- paste0(sign(test), " (p=", round(test, digits = 3), ")")

p1 <- change %>%
  ggplot(aes(x = Bone, y = value))+
  ylim(0,80)+
  labs(y = "Human Proportion (%)",
       x = "")+
  geom_line(aes(group = Sample, color = Sample), size = 2, show.legend = F)+
  geom_point(aes(fill = Sample), shape = 21, size = 5, show.legend = T, stroke = .75)+
  scale_color_viridis_d(option = "rocket")+
  scale_fill_viridis_d(option = "rocket")+
  showSignificance(c(1,2), 70, -0.1, star,
                   textParams = list(size = 4, fontface = "bold.italic"), segmentParams = list(size = 1))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 12, face = "bold"),
        plot.caption = element_text(face = "italic"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside");p1


un_ossicle <- test_m %>%
  group_by(Sample, variable) %>%
  mutate(diff = ((value[Bone %in% "Unidentified Ossicle"]-value[Bone %in% "Petrous"])/value[Bone %in% "Petrous"])*100,
         type = "Unidentified Ossicle",
         comb = "Petrous|Unidentified Ossicle") %>%
  filter(variable %in% c("Human Proportion (%)", "mt Proportion (%)",
                         "Avg Read Length (bp)", "PMD (%)", "Cont. Est.", "mt_coverage", "Clonality (%)")) %>%
  select(type, Bone, Sample, value, diff, variable, comb) %>%
  distinct(Sample, Bone, .keep_all = T)

p2_1 <- un_ossicle %>%
  filter(Bone %in% "Unidentified Ossicle") %>%
  filter(variable %in% "Human Proportion (%)") %>%
  ggplot(aes(x = variable, y = diff))+
  ylim(0,300)+
  #facet_zoom(y = diff <= 350, ylim = c(0,310), zoom.size = 1.5, horizontal = F)+
  geom_boxplot(lwd = .50, fill = "#fbfbf9", width = .5)+
  geom_hline(yintercept = 0, linetype = "dashed", size = .75, color = "red")+
  labs(x = "",
       y = "Change (%)",
       title = "Change in sequence statistics after preparing library from incus")+
  geom_point(aes(fill = Sample), shape = 21, size = 5, show.legend = T, stroke = .75)+
  scale_fill_viridis_d(option = "rocket")+
  scale_x_discrete(labels = "Human Proportion")+
  theme_bw()+
  theme(
    zoom.y = element_rect(fill = "grey90", color = "grey70"),
    axis.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold"),
    legend.background = element_rect(fill = "white"));p2_1


p2_2 <- un_ossicle %>%
  filter(Bone %in% "Unidentified Ossicle") %>%
  filter(variable %nin% c("Human Proportion (%)", "mt_coverage", "mt Proportion (%)")) %>%
  ggplot(aes(x = variable, y = diff))+
  #ylim(-10,20)+
  geom_boxplot(lwd = .50, fill = "#fbfbf9", width = .75)+
  geom_hline(yintercept = 0, linetype = "dashed", size = .75, color = "red")+
  #geom_line(aes(group = Sample, color = Sample), size = 2)+
  labs(x = "",
       y = "Change (%)",
       title = "Change in sequence statistics after preparing library from stapes")+
  geom_point(aes(fill = Sample), shape = 21, size = 5, show.legend = T, stroke = .75)+
  scale_fill_viridis_d(option = "rocket")+
  scale_x_discrete(labels = c("Clonality", "Avg Read Length", "PMD", "Cont. Est."))+
  #facet_grid(.~variable, scales = "free", space = "free")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "#121212"),
        strip.text = element_text(face = "bold", size = 11,
                                  color = "white"),
        axis.text = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        legend.background = element_rect(fill = "white"));p2_2


summary(un_ossicle$diff[un_ossicle$variable %in% "Human Proportion (%)"])
summary(un_ossicle$diff[un_ossicle$variable %in% "Avg Read Length (bp)"])
summary(un_ossicle$diff[un_ossicle$variable %in% "PMD (%)"])
summary(un_ossicle$diff[un_ossicle$variable %in% "Cont. Est."])

p2 <- ggarrange(p2_1, p2_2, ncol = 2, nrow = 1, legend = F, labels = c("B", "C"));p2

ggarrange(p1, p2_1, p2_2, ncol = 3, nrow = 1, common.legend = T, labels = c("A", "B", "C"), legend = F)

ggsave(filename = "figures/png/other/p_vs_un_2.png", height = 5, width = 15, device = "png", dpi = 300)
ggsave(filename = "figures/svg/other/p_vs_un_2.svg", height = 5, width = 15, device = svg)



# Combining Statistics ----------------------------------------------------
seqstats <- rbind(stapes, malleus, incus, un_ossicle)
write.table(seqstats, "data/seqstats.txt", sep = "\t", row.names = F, quote = F)






# Microbial Analysis -------------------------------------------------------------

# Reading krakenuniq results
raw <- read_xlsx("data/Supplementary_Tables_rev01.xlsx", na = "S1")

# ______ PCoA on Anatolian Samples ----
krakenuniq <- read_xlsx("data/Supplementary_Tables_rev01.xlsx", sheet = "S2")
krakenuniq <- krakenuniq %>%
  filter(Species %nin% "Homo sapiens") %>%
  select(Library, Species, Reads) %>%
  pivot_wider(names_from = Library, 
              values_from = Reads, 
              values_fill = 0) %>%
  column_to_rownames("Species") %>%
  as.matrix()


# Create metadata
metadata <- raw %>%
  select(`Library ID`, Bone, `Library Type`, Region)

# Edit column names
colnames(metadata) <- c("library", "bone", "library_type", "region")

# Filter metadata & Krakenuniq matrix
metadata_matched <- metadata %>%
  filter(library_type %nin% "petrous&tooth") %>%
  filter(bone %nin% "Tooth") %>%
  #filter(!grepl(library_type, pattern = "only")) %>%
  filter(region %nin% "Iberia") %>%
  ungroup()

#metadata_matched <- metadata
cols <- which(colnames(krakenuniq) %in% metadata_matched$library)
krakenuniq <- krakenuniq[,c(cols)]


# Transposing
kraken_t <- t(krakenuniq)

# Filter if number of reads equal to 0
kraken_t <- kraken_t[rowSums(kraken_t) > 0, ]

# Normalization
kraken_rel <- apply(kraken_t, 1, function(x) (x / sum(x)) * 100)
kraken_rel <- t(kraken_rel)

# Check if there is NA
any(is.na(kraken_rel))

# Match metadata and metaphlan matrix
common_samples <- intersect(rownames(kraken_rel), metadata_matched$library)
kraken_rel <- kraken_rel[common_samples, ]
metadata_matched <- metadata_matched[match(common_samples, metadata_matched$library), ]

# Check if row numbers equal or not
if(nrow(kraken_rel) != nrow(metadata_matched)) stop("ERROR")


# Calculate Bray-Curtis Distance
dist_matrix <- vegdist(kraken_rel, method = "bray")

# PCoA (Multidimentional Scaling)
pcoa_res <- cmdscale(dist_matrix, k = 2, eig = TRUE)

# Get coordinates and variance levels
pcoa_df <- as.data.frame(pcoa_res$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
var <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 1)
pcoa_df$library <- rownames(pcoa_df)

# Get bone information
pcoa_df <- merge(pcoa_df, metadata_matched, by = "library", all.x = T)

pcoa_df$bone <- factor(pcoa_df$bone, levels = c("Petrous", "Stapes", "Malleus", 
                                                "Incus", "Unidentified Ossicle", "Tooth"))

# Run PERMANOVA
metadata_matched$library_short <- substr(metadata_matched$library, start = 1, stop = 6)
permanova_res <- adonis2(dist_matrix ~ bone, data = metadata_matched, 
                         strata = metadata_matched$library_short)


# Get PERMANOVA results
r2 <- round(permanova_res$R2[1], 3); r2
f <- round(permanova_res$F[1],3);f
p_value <- round(permanova_res$`Pr(>F)`[1],3);p_value

permanova_results <- data.frame(parameter = c("r2", "f", "p", "PCoA1", "PCoA2"), 
                                value = c(r2, f, p_value, var[1], var[2]), region = "Anatolia")

# Calculate correlation between PCoA and krakenuniq relative abundance matrix
fit <- envfit(pcoa_res, kraken_rel, permutations = 999)

# Get most correlated longest vectors 
spp_scores <- as.data.frame(scores(fit, "vectors"))
spp_scores$species <- rownames(spp_scores)
spp_scores$pvals <- fit$vectors$pvals

# Filter only significant top 15 species
top_species <- spp_scores %>%
  filter(pvals < 0.05) %>%
  arrange(desc(Dim1^2 + Dim2^2)) %>%
  head(15) %>%
  mutate(region = "Anatolia")


write.table(x = top_species, file = "data/top_species", quote = F, sep = "\t", row.names = F)
write.table(x = pcoa_df, file = "data/pcoa_df", quote = F, sep = "\t", row.names = F)
write.table(x = permanova_results, file = "data/permanova_results", quote = F, sep = "\t", row.names = F)




# ______ PCoA on Iberian Samples ----
krakenuniq <- read_xlsx("data/Supplementary_Tables_rev01.xlsx", sheet = "S2")
krakenuniq <- krakenuniq %>%
  filter(Species %nin% "Homo sapiens") %>%
  select(Library, Species, Reads) %>%
  pivot_wider(names_from = Library, 
              values_from = Reads, 
              values_fill = 0) %>%
  column_to_rownames("Species") %>%
  as.matrix()


# Create metadata
metadata <- raw %>%
  select(`Library ID`, Bone, `Library Type`, Region)

# Edit column names
colnames(metadata) <- c("library", "bone", "library_type", "region")

# Filter metadata & Krakenuniq matrix
metadata_matched <- metadata %>%
  filter(library_type %nin% "petrous&tooth") %>%
  filter(bone %nin% "Tooth") %>%
  filter(!grepl(library_type, pattern = "only")) %>%
  filter(region %in% "Iberia") %>%
  ungroup()

#metadata_matched <- metadata
cols <- which(colnames(krakenuniq) %in% metadata_matched$library)
krakenuniq <- krakenuniq[,c(cols)]


# Transposing
kraken_t <- t(krakenuniq)

# Filter if number of reads equal to 0
kraken_t <- kraken_t[rowSums(kraken_t) > 0, ]

# Normalization
kraken_rel <- apply(kraken_t, 1, function(x) (x / sum(x)) * 100)
kraken_rel <- t(kraken_rel)

# Check if there is NA
any(is.na(kraken_rel))

# Match metadata and metaphlan matrix
common_samples <- intersect(rownames(kraken_rel), metadata_matched$library)
kraken_rel <- kraken_rel[common_samples, ]
metadata_matched <- metadata_matched[match(common_samples, metadata_matched$library), ]

# Check if row numbers equal or not
if(nrow(kraken_rel) != nrow(metadata_matched)) stop("ERROR")


# Calculate Bray-Curtis Distance
dist_matrix <- vegdist(kraken_rel, method = "bray")

# PCoA (Multidimentional Scaling)
pcoa_res <- cmdscale(dist_matrix, k = 2, eig = TRUE)

# Get coordinates and variance levels
pcoa_df <- as.data.frame(pcoa_res$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
var <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 1)
pcoa_df$library <- rownames(pcoa_df)

# Get bone information
pcoa_df <- merge(pcoa_df, metadata_matched, by = "library", all.x = T)

pcoa_df$bone <- factor(pcoa_df$bone, levels = c("Petrous", "Stapes", "Malleus", 
                                                "Incus", "Unidentified Ossicle", "Tooth"))

# Run PERMANOVA
metadata_matched$library_short <- substr(metadata_matched$library, start = 1, stop = 6)
permanova_res <- adonis2(dist_matrix ~ bone, data = metadata_matched, 
                         strata = metadata_matched$library_short)

# Get PERMANOVA results
r2 <- round(permanova_res$R2[1], 3); r2
f <- round(permanova_res$F[1],3);f
p_value <- round(permanova_res$`Pr(>F)`[1],3);p_value

permanova_results <- data.frame(parameter = c("r2", "f", "p", "PCoA1", "PCoA2"), 
                                value = c(r2, f, p_value, var[1], var[2]), region = "Iberia")


# Calculate correlation between PCoA and krakenuniq relative abundance matrix
fit <- envfit(pcoa_res, kraken_rel, permutations = 999)

# Get most correlated longest vectors 
spp_scores <- as.data.frame(scores(fit, "vectors"))
spp_scores$species <- rownames(spp_scores)
spp_scores$pvals <- fit$vectors$pvals

# Filter only significant top 15 species
top_species <- spp_scores %>%
  filter(pvals < 0.05) %>%
  arrange(desc(Dim1^2 + Dim2^2)) %>%
  head(15) %>%
  mutate(region = "Iberia")


write.table(x = top_species, file = "data/top_species", quote = F, sep = "\t", row.names = F, append = T, col.names = F)
write.table(x = pcoa_df, file = "data/pcoa_df", quote = F, sep = "\t", row.names = F, append = T, col.names = F)
write.table(x = permanova_results, file = "data/permanova_results", quote = F, sep = "\t", row.names = F, append = T, col.names = F)



# ______ Comparison of Non-oral and Oral reads ----
raw <- read_xlsx("data/Supplementary_Tables_rev01.xlsx", sheet = "S2")

raw <- raw %>%
  filter(Species != "bacterium") %>%
  filter(Library %nin% c("Ind039_b"))

table(raw$Source)

# Calculate relative abundances for species sources
data <- raw %>%
  group_by(Library) %>%
  mutate(sum_reads = sum(Reads)) %>%
  group_by(Library, Source) %>%
  mutate(sum_source = sum(Reads)) %>%
  distinct(Library, Source, .keep_all = T) %>%
  select(Library, Source, sum_reads, sum_source, `Tissue Type`) %>%
  rowwise() %>%
  mutate(rel_source = (sum_source/sum_reads)*100)

# Modify dataset
#data$`Tissue Type`[data$`Tissue Type` %in% "Unidentified Ossicle"] <- "UI Ossicle"
data$`Tissue Type` <- factor(data$`Tissue Type`, 
                             levels = c("Petrous", "Stapes", "Malleus", "Incus", "Unidentified Ossicle", "Tooth"))

biplot <- data %>%
  ungroup() %>%
  select(Library, `Tissue Type`, rel_source, Source) %>%
  pivot_wider(names_from = Source,
              values_from = rel_source)

# Convert NA values to 0
biplot[,c(3:5)][is.na(biplot[,c(3:5)])] <- 0


# Writing results file
write.table(biplot, "data/source_comparison", sep = "\t", row.names = F)




