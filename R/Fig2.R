## Code for Fig. 2
## Author: V. Kim

library(tidyverse)
library(ggpubr)
library(ggsignif)
setwd(here::here())

datadir = '../data/'
figdir = '../figures/'

viab_df = read.delim(paste0(datadir, 'drugresponse.tsv'))
# extract diagnosis
viab_df = mutate(viab_df, Diagnosis = gsub("(.+)(_[0-9]+)", "\\1", patientID))


ctrl_df = filter(viab_df, Drug == 'DMSO') %>%
  group_by(plate, Drug, Culture, patientID, Diagnosis) %>%
  summarise(viab = median(viab))

mono = filter(ctrl_df, Culture == "Monoculture")
co = filter(ctrl_df, Culture == "Coculture")

all_ctrl = inner_join(mono, co, by = c("plate", "patientID", "Diagnosis"))

# aggregate across the replicates
all_ctrl = group_by(all_ctrl, patientID, Diagnosis) %>%
  summarise(viab.x = mean(viab.x),
            viab.y = mean(viab.y))

# make Fig. 2A
disease_pal = c("#4daf4a", "#377eb8", "#9e4ea3", "#ffc633", "#ff5f00")
p = ggplot(all_ctrl, 
                 aes(x=viab.x, y = viab.y,
                     color = Diagnosis))+
  geom_abline(slope = 1, 
              size=0.7, 
              color="#a6a6a6") + 
  geom_point(size=1) +
  scale_color_manual(values = disease_pal)+
  labs(x = "Baseline viability in monoculture",
       y = "Baseline viability in coculture",
       title = "",
       color = "") +  
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,1),
                     breaks = seq(0.25,1,by = 0.25)) + 
  scale_x_continuous(expand = c(0,0), limits = c(0,1),
                     labels = c(0, format(seq(0.25,1,by=0.25),
                                          digits = 2))) + 
  theme_classic(base_size = 9) +
  theme(aspect.ratio = 1,
        legend.key.size = unit(3, 'mm'),
        plot.title = element_text(hjust = 0.5))
p

# positional effects
pos_df = filter(viab_df, Drug == 'DMSO') %>%
  mutate(Edge = factor(Edge, levels = 2:0,
                       labels = c('Edge', 
                                  'Next to edge', 
                                  'Center')),
         Culture = factor(Culture, levels = c("Monoculture", "Coculture")))
# make Fig. 2B
ggplot(pos_df,
       aes(x = Edge, y = viab_norm,
           fill = Edge, group = Edge)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.25) + 
  facet_wrap(~ Culture) + 
  theme_classic(base_size = 9) +
  scale_fill_manual(values = c("#448ac1", "#90d1c2", "#eef8b8")) + 
  ylab("Control viability") + 
  ylim(c(0,1.2)) + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust=0.5),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank()) 

# mean drug response
mean_drugresp = filter(viab_df, Drug != "DMSO") %>%
  group_by(plate, Drug, Culture) %>%
  summarise(viab_mean = mean(viab_norm))

spapop_df = select(ungroup(ctrl_df), plate, Culture, viab) %>%
  mutate(spapop = 1 - viab)
mean_drugresp = inner_join(mean_drugresp, spapop_df)
# correlations between SA and drug response
cor_df = group_by(mean_drugresp, Culture, Drug) %>%
  summarise(cor = cor(viab_mean, spapop, method='spearman')) %>%
  mutate(Culture = factor(Culture,
                          levels = c("Monoculture", 
                                     "Coculture"),
                          labels = c("Mono.", "Co.")))

# make Fig. 2C
pal = c('#4878d0', '#ee854a')
ggplot(cor_df, aes(x=Culture, y = cor, fill = Culture))+
  geom_boxplot(outlier.size = 0.25, lwd = 0.25) + 
  scale_fill_manual(values =  pal)+
  ylim(c(-1,1)) +
  ylab("Spearman correlation") + 
  theme_classic(base_size = 9) + 
  theme(legend.position = 'none')

# venetoclax
ggplot(filter(mean_drugresp, Drug =='Venetoclax') %>%
         mutate(Culture = factor(Culture,
                                 levels = c("Monoculture", 
                                            "Coculture"))),
       aes(x = spapop, y = viab_mean,
           color = Culture)) +
  geom_point(size = 0.5) + 
  facet_wrap(~ Culture) + 
  labs(x = "Spontaneous apoptosis",
       y = "Normalized viability",
       title = 'Venetoclax') + 
  theme_classic(base_size = 9) + 
  scale_color_manual(values = pal) + 
  stat_cor(aes(label = ..p.label..),
           method = 'spearman',
           color = 'black',
           size = 2.5, label.y = 0.1, label.x = 0.1) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.2)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1.2)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        legend.position = 'none',
        strip.background = element_blank())

# ibrutinib
ggplot(filter(mean_drugresp, Drug =='Ibrutinib') %>%
         mutate(Culture = factor(Culture,
                                 levels = c("Monoculture", 
                                            "Coculture"))),
       aes(x = spapop, y = viab_mean,
           color = Culture)) +
  geom_point(size = 0.5) + 
  facet_wrap(~ Culture) + 
  labs(x = "Spontaneous apoptosis",
       y = "Normalized viability",
       title = 'Ibrutinib') + 
  theme_classic(base_size = 9) + 
  scale_color_manual(values = pal) + 
  stat_cor(aes(label = ..p.label..),
           method = 'spearman',
           color = 'black',
           size = 2.5, label.y = 0.1, label.x = 0.1) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.2)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1.2)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        legend.position = 'none',
        strip.background = element_blank())

# fludarabine
ggplot(filter(mean_drugresp, Drug =='Fludarabine') %>%
         mutate(Culture = factor(Culture,
                                 levels = c("Monoculture", 
                                            "Coculture"))),
       aes(x = spapop, y = viab_mean,
           color = Culture)) +
  geom_point(size = 0.5) + 
  facet_wrap(~ Culture) + 
  labs(x = "Spontaneous apoptosis",
       y = "Normalized viability",
       title = 'Fludarabine') + 
  theme_classic(base_size = 9) + 
  scale_color_manual(values = pal) + 
  stat_cor(aes(label = ..p.label..),
           method = 'spearman',
           color = 'black',
           size = 2.5, label.y = 0.1, label.x = 0.1) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.2)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1.2)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        legend.position = 'none',
        strip.background = element_blank())

## now only CLL samples
cll_samples = unique(filter(ctrl_df, Diagnosis == 'CLL')$plate)
mean_drugresp = filter(mean_drugresp, plate %in% cll_samples)

# venetoclax
ggplot(filter(mean_drugresp, Drug =='Venetoclax') %>%
         mutate(Culture = factor(Culture,
                                 levels = c("Monoculture", 
                                            "Coculture"))),
       aes(x = spapop, y = viab_mean,
           color = Culture)) +
  geom_point(size = 0.5) + 
  facet_wrap(~ Culture) + 
  labs(x = "Spontaneous apoptosis",
       y = "Normalized viability",
       title = 'Venetoclax') + 
  theme_classic(base_size = 9) + 
  scale_color_manual(values = pal) + 
  stat_cor(aes(label = ..p.label..),
           method = 'spearman',
           color = 'black',
           size = 2.5, label.y = 0.1, label.x = 0.1) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.2)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1.2)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        legend.position = 'none',
        strip.background = element_blank())

# ibrutinib
ggplot(filter(mean_drugresp, Drug =='Ibrutinib') %>%
         mutate(Culture = factor(Culture,
                                 levels = c("Monoculture", 
                                            "Coculture"))),
       aes(x = spapop, y = viab_mean,
           color = Culture)) +
  geom_point(size = 0.5) + 
  facet_wrap(~ Culture) + 
  labs(x = "Spontaneous apoptosis",
       y = "Normalized viability",
       title = 'Ibrutinib') + 
  theme_classic(base_size = 9) + 
  scale_color_manual(values = pal) + 
  stat_cor(aes(label = ..p.label..),
           method = 'spearman',
           color = 'black',
           size = 2.5, label.y = 0.1, label.x = 0.1) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.2)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1.2)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        legend.position = 'none',
        strip.background = element_blank())

# fludarabine
ggplot(filter(mean_drugresp, Drug =='Fludarabine') %>%
         mutate(Culture = factor(Culture,
                                 levels = c("Monoculture", 
                                            "Coculture"))),
       aes(x = spapop, y = viab_mean,
           color = Culture)) +
  geom_point(size = 0.5) + 
  facet_wrap(~ Culture) + 
  labs(x = "Spontaneous apoptosis",
       y = "Normalized viability",
       title = 'Fludarabine') + 
  theme_classic(base_size = 9) + 
  scale_color_manual(values = pal) + 
  stat_cor(aes(label = ..p.label..),
           method = 'spearman',
           color = 'black',
           size = 2.5, label.y = 0.1, label.x = 0.1) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.2)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1.2)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        legend.position = 'none',
        strip.background = element_blank())

