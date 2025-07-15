# This script contains all codes for Sağlıcan & Sevkar et. al, 2025 "The Mini Yet Mighty Stapes: A Superior DNA Source Than The Petrous Powder" manuscript.



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
library(ggpubr)
library(ggforce)
library(rnaturalearth)
library(rnaturalearthdata)
library(superb)
library(readxl)
library(reshape2)
library(tidyverse)


# Loading Dataset
raw <- read_xlsx("data/Supplementary Table 1.xlsx", na = "..")


# Statistics on Petrous and Ossicle Libraries Prepared From the Same Indivudals-------------------------

#Filtering dataset
data <- raw %>%
  filter(`Library Type` %in% c("ossicle&petrous"))

# Changing column names
x <- match(x = c("Mean Read Length (bp)", 
                 "PMD 5'", "Contamination Estimate (contamMix)"), table = colnames(data))

colnames(data)[x] <- c("Avg Read Length (bp)", "PMD (%)", "Cont. Est.")

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
  select(Sample, Bone, `Library ID`, `Human Proportion (%)`, `Downsampled Clonality (%)`, 
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
  scale_fill_manual(values = c("#e63946", "#0A2342"))+
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
  ylim(50,95)+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c("#e63946", "#3d405b"))+
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
  ggplot(aes(x = Sample, y = value, fill = Bone))+
  geom_point(shape = 21, size = 5, stroke = .25, position = position_dodge(width = .2))+
  facet_grid(variable~., scales = "free", switch = "y")+
  ylim(.9,1)+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c("#e63946", "#3d405b"))+
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

ggsave(filename = "figures/png/p_vs_s_1.png", height = 8, width = 11, device = png, dpi = 300, create.dir = T)
ggsave(filename = "figures/svg/p_vs_s_1.svg", height = 8, width = 11, device = svg, create.dir = T)


# Calculating changing percentages
change <- test_m %>%
  filter(variable %in% "Human Proportion (%)")

temp <- change %>%
  select(Sample, Bone, value) %>%
  spread(Bone, value)

# Running Wilcoxon signed-rank test
test <- wilcox.test(temp$Petrous, temp$Stapes, paired = T)$p.value
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


stapes <- test_m %>%
  #filter(Sample %nin% c("bah032", "bah037")) %>%
  group_by(Sample, variable) %>%
  mutate(diff = ((value[Bone %in% "Stapes"]-value[Bone %in% "Petrous"])/value[Bone %in% "Petrous"])*100,
         type = "Stapes",
         comb = "Petrous|Stapes") %>%
  filter(variable %in% c("Human Proportion (%)", "mt Proportion (%)", "Avg Read Length (bp)", 
                         "PMD (%)", "Cont. Est.", "mt_coverage", "Downsampled Clonality (%)")) %>%
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

ggarrange(p1, p2_1, p2_2, ncol = 3, nrow = 1, common.legend = T, labels = c("A", "B", "C"))

ggsave(filename = "figures/png/p_vs_s_2.png", height = 5, width = 15, device = png, dpi = 300)
ggsave(filename = "figures/svg/p_vs_s_2.svg", height = 5, width = 15, device = svg)


###### Petrous vs Malleus ------------------------------------------------------
test <- data %>%
  filter(Bone %in% c("Malleus", "Petrous")) %>%
  group_by(Sample)%>%
  mutate(n=n())%>%
  filter(n>1)%>%
  select(Sample, Bone, `Library ID`, `Human Proportion (%)`, `Downsampled Clonality (%)`, 
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
  scale_fill_manual(values = c("#e63946", "#0a9dae"))+
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
  ylim(50,95)+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c("#e63946", "#0a9dae"))+
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
  ggplot(aes(x = Sample, y = value, fill = Bone))+
  geom_point(shape = 21, size = 5, stroke = .25, position = position_dodge(width = .2))+
  facet_grid(variable~., scales = "free", switch = "y")+
  ylim(.9,1)+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c("#e63946", "#3d405b"))+
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

ggsave(filename = "figures/png/p_vs_m_1.png", height = 8, width = 11, device = png, dpi = 300, create.dir = T)
ggsave(filename = "figures/svg/p_vs_m_1.svg", height = 8, width = 11, device = svg, create.dir = T)


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
                         "PMD (%)", "Cont. Est.", "mt_coverage", "Downsampled Clonality (%)")) %>%
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
  filter(Sample %nin% "bah032") %>%
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

ggarrange(p1, p2_1, p2_2, ncol = 3, nrow = 1, common.legend = T, labels = c("A", "B", "C"))

ggsave(filename = "figures/png/p_vs_m_2.png", height = 5, width = 15, device = png, dpi = 300)
ggsave(filename = "figures/svg/p_vs_m_2.svg", height = 5, width = 15, device = svg)




###### Petrous vs Incus --------------------------------------------------------

test <- data %>%
  filter(Bone %in% c("Incus", "Petrous")) %>%
  group_by(Sample)%>%
  mutate(n=n())%>%
  filter(n>1)%>%
  select(Sample, Bone, `Library ID`, `Human Proportion (%)`, `Downsampled Clonality (%)`, 
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
  scale_fill_manual(values = c("#e63946", "#c7ffda"))+
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


contamMix <- test_m %>%
  filter(variable %in% c("Cont. Est."))

p2_2 <- contamMix %>%
  ggplot(aes(x = Sample, y = value, fill = Bone))+
  geom_point(shape = 21, size = 5, stroke = .25, position = position_dodge(width = .2))+
  facet_grid(variable~., scales = "free", switch = "y")+
  ylim(.9,1)+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c("#e63946", "#c7ffda"))+
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

ggsave(filename = "figures/png/p_vs_i_1.png", height = 6, width = 11, device = png, dpi = 300)
ggsave(filename = "figures/svg/p_vs_i_1.svg", height = 6, width = 11, device = svg)


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
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside");p1


incus <- test_m %>%
  group_by(Sample, variable) %>%
  mutate(diff = ((value[Bone %in% "Incus"]-value[Bone %in% "Petrous"])/value[Bone %in% "Petrous"])*100,
         type = "Incus",
         comb = "Petrous|Incus") %>%
  filter(variable %in% c("Human Proportion (%)", "mt Proportion (%)",
                         "Avg Read Length (bp)", "PMD (%)", "Cont. Est.", "mt_coverage", "Downsampled Clonality (%)")) %>%
  select(type, Bone, Sample, value, diff, variable, comb) %>%
  distinct(Sample, Bone, .keep_all = T)

# Drawing Second Panel
p2_1 <- incus %>%
  filter(Bone %in% "Incus") %>%
  filter(variable %in% "Human Proportion (%)") %>%
  ggplot(aes(x = variable, y = diff))+
  ylim(0,80)+
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
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold"),
    legend.background = element_rect(fill = "white"));p2_1


# Drawing Third Panel
p2_2 <- incus %>%
  filter(Bone %in% "Incus") %>%
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
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        legend.background = element_rect(fill = "white"));p2_2


p2 <- ggarrange(p2_1, p2_2, ncol = 2, nrow = 1, legend = F, labels = c("B", "C"));p2

ggarrange(p1, p2_1, p2_2, ncol = 3, nrow = 1, common.legend = T, labels = c("A", "B", "C"))

ggsave(filename = "figures/png/p_vs_i_2.png", height = 5, width = 15, device = png, dpi = 300)
ggsave(filename = "figures/svg/p_vs_i_2.svg", height = 5, width = 15, device = svg)





###### Petrous vs Unidentifiable Ossicle --------------------------------------------------
test <- data %>%
  filter(Bone %in% c("Unidentifiable Ossicle", "Petrous")) %>%
  group_by(Sample)%>%
  mutate(n=n())%>%
  filter(n>1)%>%
  select(Sample, Bone, `Library ID`, `Human Proportion (%)`, `Downsampled Clonality (%)`, 
         `Avg Read Length (bp)`, `PMD (%)`, `Cont. Est.`)%>%
  ungroup()

test_m <- melt(test, id.vars = c("Sample", "Bone", "Library ID"))
test_m$Bone <- factor(test_m$Bone, levels = c("Petrous", "Unidentifiable Ossicle"))
test_m$value <- as.numeric(test_m$value)

p1 <- test_m %>%
  filter(variable %nin% c("mt_coverage", "Cont. Est.", "Avg Read Length (bp)")) %>%
  ggplot(aes(x = Sample, y = value, fill = Bone))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_grid(variable~., scales = "free", switch = "y")+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c("#e63946", "#f58f29"))+
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
  ylim(50,95)+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c("#e63946", "#f58f29"))+
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
  ggplot(aes(x = Sample, y = value, fill = Bone))+
  geom_point(shape = 21, size = 5, stroke = .25, position = position_dodge(width = .2))+
  facet_grid(variable~., scales = "free", switch = "y")+
  ylim(.9,1)+
  labs(x = "",
       y = "",
       title = "")+
  scale_fill_manual(values = c("#e63946", "#f58f29"))+
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

ggsave(filename = "figures/png/p_vs_un_1.png", height = 6, width = 11, device = png, dpi = 300)
ggsave(filename = "figures/svg/p_vs_un_1.svg", height = 6, width = 11, device = svg)

# Testing changing percentage
change <- test_m %>%
  filter(variable %in% "Human Proportion (%)")

temp <- change %>%
  select(Sample, Bone, value) %>%
  spread(Bone, value)

test <- wilcox.test(temp$Petrous, temp$`Unidentifiable Ossicle`, paired = T)$p.value

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
  mutate(diff = ((value[Bone %in% "Unidentifiable Ossicle"]-value[Bone %in% "Petrous"])/value[Bone %in% "Petrous"])*100,
         type = "Unidentifiable Ossicle",
         comb = "Petrous|Unidentifiable Ossicle") %>%
  filter(variable %in% c("Human Proportion (%)", "mt Proportion (%)",
                         "Avg Read Length (bp)", "PMD (%)", "Cont. Est.", "mt_coverage", "Downsampled Clonality (%)")) %>%
  select(type, Bone, Sample, value, diff, variable, comb) %>%
  distinct(Sample, Bone, .keep_all = T)

p2_1 <- un_ossicle %>%
  filter(Bone %in% "Unidentifiable Ossicle") %>%
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
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold"),
    legend.background = element_rect(fill = "white"));p2_1


p2_2 <- un_ossicle %>%
  filter(Bone %in% "Unidentifiable Ossicle") %>%
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
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        strip.placement = "outside",
        legend.background = element_rect(fill = "white"));p2_2


summary(un_ossicle$diff[un_ossicle$variable %in% "Human Proportion (%)"])
summary(un_ossicle$diff[un_ossicle$variable %in% "Avg Read Length (bp)"])
summary(un_ossicle$diff[un_ossicle$variable %in% "PMD (%)"])
summary(un_ossicle$diff[un_ossicle$variable %in% "Cont. Est."])

p2 <- ggarrange(p2_1, p2_2, ncol = 2, nrow = 1, legend = F, labels = c("B", "C"));p2

ggarrange(p1, p2_1, p2_2, ncol = 3, nrow = 1, common.legend = T, labels = c("A", "B", "C"))

ggsave(filename = "figures/png/p_vs_un_2.png", height = 5, width = 15, device = png, dpi = 300)
ggsave(filename = "figures/svg/p_vs_un_2.svg", height = 5, width = 15, device = svg)


# Combining Statistics ----------------------------------------------------
seqstats <- rbind(stapes, malleus, incus, un_ossicle)
write.table(seqstats, "data/seqstats.txt", sep = "\t", row.names = F, quote = F)





