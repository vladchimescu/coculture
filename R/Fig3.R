## Code for Fig. 3
## Author: V. Kim

library(tidyverse)
library(ggpubr)
library(ggsignif)
setwd(here::here())

datadir = '../data/'
figdir = '../figures/'

# FDR and effect size thresholds
fdr = 0.01
thresh = 5

# effect size based on AUC's of median dose-response curves
compute_efsize <- function(viab_df) {
  efsize_df = group_by(viab_df, Drug, Concentration,
                       Concentration_step, Culture) %>%
    summarise(viab_median = median(viab_norm)) %>%
    mutate(Concentration = as.numeric(Concentration)) %>%
    group_by(Drug, Culture) %>%
    summarise(AUC = DescTools::AUC(log(Concentration), viab_median)) %>%
    group_by(Drug) %>%
    summarise(efsize_AUC = 100*(AUC[Culture == 'Coculture'] - AUC[Culture == 'Monoculture'])/AUC[Culture=='Monoculture'])

  efsize_df
}

compute_pval <- function(viab_df) {
  # only paired coculture - monoculture observations
  viab_df = group_by(viab_df, plate, Drug, Concentration) %>%
    mutate(n=n()) %>%
    filter(n>1) %>%
    ungroup()
  
  conc_df = group_by(viab_df, Drug, Concentration) %>%
    summarise(var = var(viab_norm)) %>%
    ungroup() %>%
    group_by(Drug) %>%
    filter(var == max(var))
  
  viab_varmax = inner_join(viab_df, conc_df)
  
  # the sample is paired (same sample, same compound tested in mono- and co-culture)
  res_ttest = group_by(viab_varmax, Drug, Concentration_step) %>%
    summarise(pval = t.test(x=viab_norm[Culture == 'Monoculture'],
                            y= viab_norm[Culture=='Coculture'],
                            paired = T)$p.value,
              tval = t.test(x=viab_norm[Culture == 'Coculture'],
                            y= viab_norm[Culture=='Monoculture'],
                            paired = T)$statistic,
              viab_median_mono = median(viab_norm[Culture == 'Monoculture']),
              viab_median_co =  median(viab_norm[Culture=='Coculture'])) %>%
    mutate(efsize_point = 100*(viab_median_co - viab_median_mono) / viab_median_mono )
  
  res_ttest = rename(res_ttest, concstep = Concentration_step)
  res_ttest = mutate(ungroup(res_ttest), 
                     padj = p.adjust(pval, method = 'BH'))
  res_ttest
}

viab_df = read.delim(paste0(datadir, 'drugresponse.tsv'))
# extract diagnosis
viab_df = mutate(viab_df, Diagnosis = gsub("(.+)(_[0-9]+)", "\\1", patientID))

drugconc = distinct(viab_df, Drug, Concentration) %>%
  group_by(Drug) %>%
  mutate(Concentration_step = ifelse(Drug != "DMSO", rank(Concentration), 0))
# add concentration step
viab_df = inner_join(viab_df, drugconc)

# for this analysis subset only to drug response data (no control wells)
viab_df = filter(viab_df, Drug != "DMSO")

# subset to CLL samples
viab_cll = filter(viab_df, Diagnosis == 'CLL')
# compute the effect size of stromal effects
efsize_cll = compute_efsize(viab_cll)

# exclude toxic conditions
toxic_df = read.delim(paste0(datadir, 'toxic_conditions.tsv'))
# remove concentrations toxic to stroma
viab_cll = anti_join(viab_cll, toxic_df)

# compute p-values for CLL
pval_cll = compute_pval(viab_cll)
res_cll = inner_join(pval_cll, efsize_cll)

res_cll = mutate(res_cll, type_inter = ifelse(padj < fdr & efsize_AUC > thresh, "Decreased",
                                                  ifelse(padj < fdr & efsize_AUC < -thresh, "Increased", "Unchanged")))

# load compound class annotation
drugclass = readxl::read_xlsx(paste0(datadir, "Drugs_Pathways.xlsx"))
drugclass = mutate(drugclass, Drug = plyr::mapvalues(Drug,
                                                     from = c("IRAK4 Inhibitor, Compound 26"),
                                                     to = c("IRAK4 inhibitor")))
res_cll = inner_join(res_cll, drugclass)

levs = rev(c("Decreased",
             "Unchanged",
             "Increased"))

class_ord = group_by(res_cll, class) %>%
  summarise(mean_efsize = mean(efsize_AUC)) %>%
  arrange(mean_efsize)

# plot by compound subclass and color by compound efficacy
res_cll = mutate(res_cll, 
                     efsize_AUC=ifelse(Drug %in% c('Carfilzomib', 'JQ1'),
                                       efsize_AUC-25, efsize_AUC))

viab_aml = filter(viab_df, Diagnosis == 'AML')

# compute the effect size of stromal effects
efsize_aml = compute_efsize(viab_aml)

# remove concentrations toxic to stroma
viab_aml = anti_join(viab_aml, toxic_df)

# compute p-values for aml
pval_aml = compute_pval(viab_aml)
res_aml = inner_join(pval_aml, efsize_aml)

res_aml = mutate(res_aml, type_inter = ifelse(padj < fdr & efsize_AUC > thresh, "Decreased",
                                              ifelse(padj < fdr & efsize_AUC < -thresh, "Increased", "Unchanged")))
aml_efsize = inner_join(res_aml, select(drugclass, -Comment))

aml_vs_cll = inner_join(aml_efsize, res_cll, by=c("Drug", "superclass", "class")) %>%
  mutate(type_inter = ifelse(type_inter.x == 'Unchanged' & type_inter.y != 'Unchanged',
                             ifelse((sign(efsize_AUC.x) == sign(efsize_AUC.y)) & abs(efsize_AUC.x) > 5, type_inter.y, type_inter.x),
                             type_inter.x))

aml_res = select(aml_vs_cll, Drug, type_inter, efsize_AUC.x, pval.x, superclass, class) %>%
  rename(efsize_AUC = efsize_AUC.x,
         pval = pval.x)
res_aml = mutate(aml_res, 
                     efsize_AUC=ifelse(Drug %in% c('Carfilzomib'),
                                       efsize_AUC-75, efsize_AUC))

res_cll$entity = 'CLL + stroma coculture'
res_aml$entity = 'AML + stroma coculture'

pad = 3
p = ggplot(bind_rows(res_cll, res_aml) %>%
         mutate(entity = factor(entity,
                                levels = c("CLL + stroma coculture",
                                           "AML + stroma coculture"))),
       aes(y = factor(class, levels = class_ord$class),
           x = efsize_AUC,
           color=factor(type_inter, levels=levs))) +
  geom_vline(xintercept = -5, linetype='dotted')+
  geom_vline(xintercept = 5, linetype='dotted')+
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = rev(c("#d8345f", "#ccafaf", "#588da8"))) +
  #scale_size_continuous(range = c(0.5,4))+
  scale_x_continuous(expand = c(0,0),
                     limits = c(min(res_aml$efsize_AUC)-pad, max(res_aml$efsize_AUC)+pad),
                     breaks = c(0, 25, 50, 75),
                     labels = c(0, 25, 50, 100))+
  facet_wrap(~ entity) + 
  labs(y = "",
       x = "Effect size %",
       title = "",
       color = "Efficacy")+
  theme_bw() +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_text(hjust = 1,
                                   size = 9),
        axis.text.x = element_text(size = 9),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 9),
        #legend.position = c(0.8,0.3),
        legend.box.background = element_rect(colour = "black"))
ggsave(p+theme(legend.position = 'none'), 
       filename = paste0(figdir, "efsize-by-drugclass-CLL-and-AML-without-pointsize.pdf"),
       height = 11.5, width = 14, units = 'cm')
