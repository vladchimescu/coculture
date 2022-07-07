## Code for Fig. 4
## Author: V. Kim

library(tidyverse)
library(ggpubr)
library(ggsignif)
setwd(here::here())

datadir = '../data/'
figdir = '../figures/'

# function to find p-value threshold at FDR level alpha
get_bh_threshold <- function (pvals, alpha, mtests = length(pvals)) 
{
  m <- length(pvals)
  pvals <- sort(pvals)
  prejected <- which(pvals <= (1:m)/mtests * alpha)
  ifelse(length(prejected) == 0, 0, pvals[prejected[which.max(prejected)]])
}

viab_df = read.delim(paste0(datadir, 'drugresponse.tsv'))
# extract diagnosis
viab_df = mutate(viab_df, Diagnosis = gsub("(.+)(_[0-9]+)", "\\1", patientID))

# for this analysis subset only to drug response data (no control wells)
viab_df = filter(viab_df, Drug != "DMSO")
# subset to CLL samples
viab_df = filter(viab_df, Diagnosis == 'CLL')

# add concentration step
drugconc = distinct(viab_df, Drug, Concentration) %>%
  group_by(Drug) %>%
  mutate(Concentration_step = ifelse(Drug != "DMSO", rank(Concentration), 0))
# add concentration step
viab_df = inner_join(viab_df, drugconc)

# load compound class annotation
drugclass = readxl::read_xlsx(paste0(datadir, "Drugs_Pathways.xlsx"))

# load sample genomic information
mutdata = read.delim(paste0(datadir, 'CLL_genomic.tsv'))

viab_mut = inner_join(viab_df,  mutdata, by = "patientID")

# here pick maximum variance concentration
conc_df = group_by(viab_df, Drug, Concentration) %>%
  summarise(var = var(viab_norm)) %>%
  ungroup() %>%
  group_by(Drug) %>%
  filter(var == max(var))
viab_varmax = inner_join(viab_df, conc_df)

# only paired genetic features (mutated - wildtype)
# and at least 5 mutated cases
feats_paired = group_by(mutdata, mutation, value) %>%
  filter(!is.na(value)) %>%
  summarise(n_val = n()) %>%
  filter(n_val > 4) %>%
  group_by(mutation) %>%
  summarise(n=n()) %>%
  filter(n == 2)

viab_mut = filter(viab_mut, mutation %in% unique(feats_paired$mutation))
viab_mut = filter(viab_mut, !is.na(value))

res_ttest = mutate(viab_mut, value = plyr::mapvalues(value, from=c("U", "M", "wt", "Mut"),
                                                     to=c(0,1,0,1))) %>%
  group_by(Drug, Concentration, Culture, mutation) %>%
  do(broom::tidy(t.test(viab_norm ~ value, data=., paired=F))) %>%
  select(mutation, p.value, statistic)
res_ttest = mutate(ungroup(res_ttest),
                   padj = p.adjust(p.value, method='BH'))

res_out = filter(res_ttest, padj <= 0.1) %>%
  select(Drug, Concentration, mutation) %>%
  inner_join(res_ttest)

res_out = inner_join(filter(res_out, Culture == 'Monoculture'),
                     filter(res_out, Culture == 'Coculture'),
                     by = c("Drug", "Concentration", "mutation")) %>%
  select(-c(Culture.x, Culture.y, p.value.x, p.value.y)) %>%
  mutate(Significant =ifelse(padj.x <= 0.1, ifelse(padj.y <= 0.1, "Both", "Monoculture"), 'Coculture')) %>%
  select(-c(padj.x, padj.y)) %>%
  arrange(mutation, Drug, Concentration, Significant)
res_out = arrange(res_out, Drug, mutation)
xlsx::write.xlsx(res_out, file=paste0(datadir, "drug-gene-associations.xlsx"))

pth = get_bh_threshold(res_ttest$p.value, alpha = 0.1)
t_thresh = abs(filter(res_ttest, p.value == pth)$statistic)

drugmut_class = inner_join(filter(res_ttest, Culture == 'Monoculture'),
                           filter(res_ttest, Culture == 'Coculture'),
                           by=c("Drug", "mutation", "Concentration")) %>%
  inner_join(drugclass)
# remove drug class 'Other' and those that don't have any significant hits
drugmut_class = filter(drugmut_class, !superclass %in% c("Other"))

# color drug-gene associations below 10% FDR
drugmut_class = dplyr::mutate(drugmut_class, 
                              mutation = ifelse(padj.x > 0.1 & padj.y > 0.1, "FDR > 0.1", mutation))
# drugmut_class = dplyr::group_by(drugmut_class, Drug, mutation, superclass) %>%
#   summarise(statistic.x = statistic.x[which.max(abs(statistic.x))],
#             statistic.y = statistic.y[which.max(abs(statistic.y))])
levs = c("ATM", "del11q", "IGHV", "TP53", "trisomy12", "FDR > 0.1")
df_out = filter(drugmut_class, mutation != "FDR > 0.1") %>%
  dplyr::group_by(Drug, mutation, superclass) %>%
  summarise(statistic.x = statistic.x[which.max(abs(statistic.x))],
            statistic.y = statistic.y[which.max(abs(statistic.y))])
p = ggplot(df_out,
           aes(x = statistic.x,
               y = statistic.y,
               color = factor(mutation, levels = levs))) +
  scale_color_manual(name = "Mutations",
                     values = c("#75cfb8", "#eb596e", 
                                "#9a181a", "#3d7ea6", 
                                "#fce38a")) + 
  geom_point(data=filter(drugmut_class, mutation == 'FDR > 0.1'),
             color = "#C0C0C0", shape = 1, alpha=0.8) +
  geom_point() +
  ggpubr::stat_cor(data = drugmut_class,
                   aes(x = statistic.x,
                       y = statistic.y,
                       label = ..r.label..), 
                   inherit.aes = F, size = 5*8 / 14)+
  facet_wrap(~ superclass, scales = 'fixed')+
  xlim(c(-5.5,5.5)) + ylim(c(-5.5,5.5))+
  labs(x = "t-value in monoculture",
       y = "t-value in coculture",
       color = "") + 
  geom_abline(slope = 1, linetype = 'dotted') + 
  theme_bw(base_size = 9) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        #strip.background = element_blank(),
        legend.position = 'none')
ggsave(p, filename = paste0(figdir, "tval-mono-coculture-by-drugclass.pdf"),
       width = 10, height = 10, units = 'cm')

# supplementary Fig. 6
co_signs = filter(res_ttest, Culture == 'Coculture') %>%
  mutate(sgn = sign(statistic))

mono_signs = filter(res_ttest, Culture == 'Monoculture') %>%
  mutate(sgn = sign(statistic))

signif_signs = inner_join(mono_signs, 
                          co_signs,
                          by = c("Drug",
                                 "Concentration",
                                 "mutation")) %>%
  filter(padj.x <= 0.1 | padj.y <= 0.1) %>%
  mutate(sgn.x = ifelse(sgn.x == -1, "Negative", "Positive")) %>%
  mutate(sgn.y = ifelse(sgn.y == -1, "Negative", "Positive"))

p = as.data.frame(table(signif_signs$sgn.x, signif_signs$sgn.y)) %>%
  ggplot(aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), size = 9*5/14)+
  labs(x = "Monoculture", 
       y = "Coculture")+
  scale_fill_gradient(low='white', high='#B91646') +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme_classic(base_size = 9)+
  theme(legend.position = 'none',
        aspect.ratio = 1)
ggsave(p, filename = paste0(figdir, 'druggene-sign-concordance.pdf'),
       width = 6, height = 5, units = 'cm')


df = filter(viab_mut, Drug == 'Ibrutinib' & mutation == 'IGHV' & Concentration_step==1) 
df = mutate(df, Culture = factor(Culture, levels=c("Monoculture", "Coculture")))
p = ggplot(df, 
           aes(x = factor(value, levels = c("U", "M")), 
               y = viab_norm,
               fill = Culture))+
  geom_dotplot(aes(color = Culture),
               stackdir = 'center',
               binaxis = 'y',
               dotsize=0.75)+
  geom_boxplot(alpha=0.5,lwd = 0.25, outlier.colour = NA) +
  scale_fill_manual(values = c('#4878d0', '#ee854a')) +
  scale_color_manual(values = c('#4878d0', '#ee854a')) + 
  facet_wrap(~ Culture) + 
  ggpubr::stat_compare_means(comparisons = list(c("M", "U"),
                                                c("U", "M")),
                             method = "t.test",
                             paired = F, size=1.5)+
  labs(x = "IGHV", y = "Normalized viability",
       title = expression(paste("Ibrutinib (0.04 ", mu, "M)")))+ 
  theme_classic(base_size = 9)+
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = rel(0.9)),
        legend.position = 'none')
ggsave(p, filename = paste0(figdir, 'ibrutinib-IGHV.pdf'),
       width = 6.5, height = 4.5, units = 'cm')

df = filter(viab_mut, Drug %in% c('Nutlin 3a') & mutation == 'TP53' & Concentration_step == 3)
df = mutate(df, Culture = factor(Culture, levels=c("Monoculture", "Coculture")))
p = ggplot(df, 
           aes(x = factor(value, levels = c("wt", "Mut")), 
               y = viab_norm,
               fill = Culture))+
  geom_dotplot(aes(color = Culture),
               stackdir = 'center',
               binaxis = 'y',
               dotsize=0.75)+
  geom_boxplot(alpha=0.5,lwd = 0.25, outlier.colour = NA) +
  scale_fill_manual(values = c('#4878d0', '#ee854a')) +
  scale_color_manual(values = c('#4878d0', '#ee854a')) + 
  facet_wrap(~ Culture) + 
  ggpubr::stat_compare_means(comparisons = list(c("Mut", "wt"),
                                                c("wt", "Mut")),
                             method = "t.test",
                             paired = F,
                             size =1.5)+
  labs(x = "TP53", y = "Normalized viability",
       title = expression(paste("Nutlin 3a (9 ", mu, "M)"))) + 
  theme_classic(base_size = 9)+
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = rel(0.9)),
        legend.position = 'none')
ggsave(p, filename = paste0(figdir, 'nutlin3a-TP53.pdf'),
       width = 6.5, height = 4.5, units = 'cm')

# compute the absolute value of the effect size of drug gene associations
res_efsize = mutate(viab_mut, value = plyr::mapvalues(value, from=c("U", "M", "wt", "Mut"),
                                                      to=c(0,1,0,1))) %>%
  group_by(Drug, Concentration, Culture, mutation) %>%
  summarise(mu0 = mean(viab_norm[value == 0]),
            mu1 = mean(viab_norm[value == 1])) %>%
  mutate(efsize_abs = 100 * abs(mu0 - mu1) / mu0) %>%
  select(-c(mu0, mu1))

res_efsize = inner_join(res_efsize, res_ttest)
drugs_sel = c("Ibrutinib",
              "Idelalisib", 
              "Duvelisib",
              "PRT062607",
              "Dasatinib", 
              "Selumetinib")
df_efsize = filter(res_efsize, Drug %in% drugs_sel) %>%
  filter(mutation %in% c("IGHV", "trisomy12")) %>%
  filter(efsize_abs < 50)

df_efsize = inner_join(df_efsize,
                       distinct(viab_mut, Drug,
                                Concentration, Concentration_step)) %>%
  mutate(conc = plyr::mapvalues(Concentration_step, from = 1:3,
                                to=c("Lowest", "Medium", "Highest")))

p = ggplot(df_efsize, aes(x = factor(Drug, levels = rev(drugs_sel)), 
                      y = efsize_abs,
                      color = factor(Culture, levels = c("Monoculture", "Coculture")),
                      fill = factor(Culture, levels = c("Monoculture", "Coculture")))) + 
  geom_point(size = 4, shape = 108)+
  scale_color_manual(values =  c('#4878d0', '#ee854a'))+
  scale_fill_manual(values =  c('#4878d0', '#ee854a'))+
  scale_alpha_manual(values = c(0.4, 0.75, 1))+
  facet_wrap(~ mutation, scales = "free_x") +
  coord_flip() + 
  labs(fill = "",
       color = "",
       y = "Effect size %",
       x = "",
       title = "") +
  # scale_y_continuous(expand = c(0,0))+
  # scale_x_discrete(expand=c(0,0))+
  theme_classic(base_size = 9) + 
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(size = 0.25, colour = "#a6a6a6"),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text())
ggsave(p+theme(legend.position = 'none'),
       filename = paste0(figdir, "kinase-inhibitor-IGHV-trisomy12-tickplot.pdf"),
       width = 7, height = 7, units = 'cm')


plot_df = inner_join(viab_mut, drugclass) %>%
  filter(class %in% c("Syk inhibitor", "BTK inhibitor", "PI3K inhibitor")) %>%
  filter(mutation %in% c("IGHV", "trisomy12"))

plot_df = mutate(plot_df, value = plyr::mapvalues(value, from = c("U", "M", "wt", "Mut"),
                                                  to = c("Wt (U)", "Mut (M)", "Wt (U)", "Mut (M)")))

plot_df = mutate(plot_df, mutgroup = ifelse(mutation == 'IGHV', ifelse(value == 'Wt (U)', 'U-CLL', 'M-CLL' ), 
                                            ifelse(value == 'Wt (U)', 'Tri12(-)', 'Tri12(+)')))
p = ggplot(plot_df, aes(x = factor(mutgroup, levels = c("U-CLL", "M-CLL", "Tri12(-)", "Tri12(+)")), 
                        y = viab_norm, 
                        fill = factor(Culture, levels =c ("Monoculture", "Coculture")))) +
  geom_boxplot(outlier.shape=NA,
               lwd = 0.25,
               width=0.5,
               position = position_dodge(width = 0.8)) +
  facet_wrap(~ mutation, scales = 'free_x') + 
  scale_fill_manual(values = c('#4878d0', '#ee854a')) + 
  labs(x = "", y = "Normalized viability",
       title = "",
       fill = "")+ 
  theme_classic(base_size = 9) + 
  theme(plot.title = element_text(size = 9, hjust = 0.5),
        strip.background = element_blank())
ggsave(p+theme(legend.position = 'none'),
       filename = paste0(figdir, "viabnorm-vs-mutation-BCR-inhibitors.pdf"),
       width = 5, height = 7, units = 'cm')

drugs_sel = c("Ibrutinib",
              "Idelalisib", 
              "Duvelisib",
              "PRT062607",
              "Dasatinib", 
              "Selumetinib",
              "Nutlin 3a",
              "Carfilzomib", 
              "Ixazomib")
df_var = filter(viab_df, Drug %in% drugs_sel)
p = ggplot(df_var, aes(y = factor(Drug, levels = rev(drugs_sel)),
                       x= viab_norm,
                       fill = factor(Culture, levels = c("Monoculture", "Coculture"))))+
  geom_boxplot(width=0.7, coef=0,
               lwd=0.25,
               outlier.shape = NA) +
  scale_fill_manual(values =  c('#4878d0', '#ee854a'))+
  labs(fill = "",
       x = "Normalized viability",
       y = "",
       title = "") +
  theme_classic(base_size = 9) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1.1))+
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size =9),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text())
ggsave(p+theme(legend.position = 'none'),
       filename = paste0(figdir, "drugcult-variability.pdf"),
       width = 5.75, height = 7, units = 'cm')
