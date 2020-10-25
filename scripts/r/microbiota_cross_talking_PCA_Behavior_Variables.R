#!/usr/bin/env Rscript

#############################################################################
### Jose Espinosa-Carrasco CB-CRG Group. Sep 2020                         ###
#############################################################################
### Cross-talking between behavior and microbiota                         ###
#############################################################################

library(ggplot2)
library(tidyverse)
library(dplyr)

## Functions
transpose_df <- function(df) {
  # keep the first column 
  names <-  df[,1]
  
  # Transpose everything other than the first column
  df.T <- as.data.frame(as.matrix(t(df[,-1])))
  
  # Assign first column as the column names of the transposed dataframe
  colnames(df.T) <- names
  return(df.T)
}

first_up <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

### microbiota data
# axis_text_size<-24;
# taxon <- "phylum"; sep_f=";";
# first_taxon <- 'Verrucomicrobia';last_taxon <- 'Actinobacteria';

axis_text_size<-24;
taxon <- "family"; sep_f=";"
first_taxon <- 'Alcaligenaceae';last_taxon <- 'Others';

# axis_text_size<-14;
# taxon <- "genus"; sep_f=";"
# first_taxon <- 'Acetatifactor';last_taxon <- 'Tyzzerella';

home_dir <- Sys.getenv("HOME")
rel_abundance_by_taxon <- paste0(home_dir, "/git/food_addiction_analysis/data/microbiota/relative_abundances_by_", taxon, ".csv")

microbiota_by_taxon_ori <- read.csv(rel_abundance_by_taxon,
                                     dec=",",
                                     sep=sep_f,
                                     # sep="\t",
                                     check.names = F,
                                     stringsAsFactors = F)

microbiota_by_taxon_tmp <- transpose_df(microbiota_by_taxon_ori)
write.csv(microbiota_by_taxon_tmp, paste0(home_dir, "/tmp.csv"))
microbiota_by_taxon <- read.csv(paste0(home_dir, "/tmp.csv"),
                                 dec=",",
                                 check.names = F,
                                 stringsAsFactors = F)
head(microbiota_by_taxon)

## Only addicts
microbiota_by_taxon <-subset(microbiota_by_taxon, Grouping=="Addict")

# Only behavioral variables used in PCA
behavioral_data_path <- paste0(home_dir,
                               "/git/food_addiction_analysis/data/microbiota/behavioral_data_from_results_microbiota_16_04_20_PCA_variables.csv")
behavioral_data <- read.csv(behavioral_data_path,
                            dec=",",
                            sep=";",
                            check.names = F,
                            stringsAsFactors = F)
# head(behavioral_data)

####################################
## Only dataframe selected columns
# head(microbiota_by_phylum)
# head(behavioral_data)
microbiota_relAbund <- subset(microbiota_by_taxon, select=-c(Grouping))
# microbiota_relAbund <- subset(microbiota_by_phylum)
behavioral_cont_data <- behavioral_data
microbio_behavioral_merged <- merge (microbiota_relAbund, behavioral_cont_data, by= "mouse_id")
head(microbio_behavioral_merged)
head(behavioral_cont_data)

## taxa
data <- gather(microbio_behavioral_merged, taxon, microbio_rel_ab, first_taxon:last_taxon)%>%
        gather(behavior_idx, index, Persistence_LP:Aversive_LP)

# head(data)
data_nest <- group_by(data, taxon, behavior_idx) %>% nest()
# head(data_nest)

# data_nest <- group_by(data, city, telecon) %>% nest()
# cor_method <- "pearson"
cor_method <- "spearman"
cor_fun <- function(df) cor.test(df$microbio_rel_ab, df$index, method = cor_method) %>% tidy()

library(broom) # Convert results of statistical functions (lm, t.test, cor.test, etc.) into tidy tables

# library(fs)
# library(lubridate)

data_nest <- mutate(data_nest, model = map(data, cor_fun))
# data_nest

# str(slice(data_nest, 1))

corr_pr <- select(data_nest, -data) %>% unnest()
corr_pr <- mutate(corr_pr, sig = ifelse(p.value <0.05, "Sig.", "Non Sig."))

######
# FDR
# test_p <- corr_pr$p.value
# class (corr_pr$p.value)
# 
# test_p<-c(test_p,0.001459265)
# p.adjust(test_p, method = 'BH', n = length(test_p))
# # corr_pr$fdr <- p.adjust(corr_pr$p.value, method = 'BH', n = length(corr_pr$p.value))
# # corr_pr <- mutate(corr_pr, sig = ifelse(fdr <0.05, "Sig.", "Non Sig."))
# subset (corr_pr,  p.value > 0.9)

hm <- ggplot() + geom_tile(data = corr_pr,
                           # aes(taxon, behavior_idx, fill = estimate),
                           aes( behavior_idx, taxon, fill = estimate),
                           size = 1,
                           colour = "white")+
  geom_tile(data = filter(corr_pr, sig == "Sig."),
            # aes(taxon, behavior_idx),
            aes(behavior_idx, taxon),
            size = 1,
            colour = "black",
            fill = "transparent") +
  geom_text(data = corr_pr, 
            # angle = 270,
            # aes(taxon, behavior_idx, #label = round(estimate, 2),
            aes(behavior_idx, taxon,    
                ## print cor estimate
                label = ifelse(sig == "Sig.", round(estimate, 2),""),
                ## print p-value
                # label = ifelse(sig == "Sig.", round(p.value, 4),""),
                fontface = ifelse(sig == "Sig.", "bold", "plain")))+
  # scale_fill_gradient2(breaks = seq(-1, 1, 0.2))+
  scale_fill_gradient2(breaks = seq(-1, 1, 0.2), #midpoint = mid, 
                       low = "#d53e4f", mid = "white",
                       high = "#abdda4" ) + 
  labs(y = first_up(taxon), x = "", fill = "", p.value = "") +
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size=axis_text_size),
        axis.title= element_text(size=axis_text_size),
        axis.text.x = element_text(angle = 90))

hm

out_dir <- paste0(home_dir, "/git/food_addiction_analysis/figures/cross_talking_microbio/")
dpi_q <- 200
extension_img <- ".png"
ggsave (hm, file=paste0(out_dir, "heatmap_fewVar_", taxon, extension_img), 
        width = 20, height = 12, dpi=dpi_q)







cor.test(microbio_behavioral_merged$Aversive_MP, microbio_behavioral_merged$Deferribacteres,
         method = "spearman")

ggplot (microbio_behavioral_merged, aes(x=Aversive_MP, y=Deferribacteres)) + 
  geom_point() +
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank())
