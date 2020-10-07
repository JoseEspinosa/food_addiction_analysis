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
# taxon <- "phylum"; sep_f="\t"; first_taxon <- 'Verrucomicrobia'; last_taxon <- 'Actinobacteria'
taxon <- "family"; sep_f=";"; first_taxon <- 'Alcaligenaceae';last_taxon <- 'Others'
# taxon <- "genus"; sep_f=";"; first_taxon<-'Acetatifactor'; last_taxon<-'Tyzzerella';

home_dir <- Sys.getenv("HOME")
rel_abundance_by_taxon <- paste0(home_dir, "/git/food_addiction_analysis/data/microbiota/relative_abundances_by_", taxon, ".csv")

microbiota_by_taxon_ori <- read.csv(rel_abundance_by_taxon,
                                    dec=",",
                                    sep=sep_f,
                                    check.names = F,
                                    stringsAsFactors = F)

microbiota_by_taxon_tmp <- transpose_df(microbiota_by_taxon_ori)
write.csv(microbiota_by_taxon_tmp, paste0(home_dir, "/tmp.csv"))
microbiota_by_taxon <- read.csv(paste0(home_dir, "/tmp.csv"),
                                dec=",",
                                check.names = F,
                                stringsAsFactors = F)
# head(microbiota_by_taxon)

# microRNAs data
# discovery
miRNAs_data_path_discovery <- paste0(home_dir,
                                     "/git/food_addiction_analysis/data/microbiota/miRNAs/resultats_miRNAs_FA_14.04.20.csv")

miRNAs_data_discovery <- read.csv(miRNAs_data_path_discovery,
                                  dec=",",
                                  sep=";",
                                  check.names = F,
                                  stringsAsFactors = F)

miRNAs_data_path_replica <- paste0(home_dir,
                                   "/git/food_addiction_analysis/data/microbiota/miRNAs/resultats_miRNAs_FA_14.04.20_replica.csv")

miRNAs_data_replica <- read.csv(miRNAs_data_path_replica,
                                dec=",",
                                sep=";",
                                check.names = F,
                                stringsAsFactors = F)
miRNAs_data_to_transp <- merge(miRNAs_data_replica, 
                               miRNAs_data_discovery, 
                               by="miRNA")

miRNAs_data <- transpose_df(miRNAs_data_to_transp)
miRNAs_data$mouse_id <- row.names(miRNAs_data)

## Select miRNAs from Elena's table
miRNAs_data_selected <- subset(miRNAs_data, select=c('mouse_id',
                                                     'mmu-miR-876-5p',
                                                     'mmu-miR-211-5p',
                                                     'mmu-miR-3085-3p',
                                                     'mmu-miR-665-3p',
                                                     'mmu-miR-3072-3p',
                                                     'mmu-miR-124-3p',
                                                     'mmu-miR-29c-3p',
                                                     'mmu-miR-544-3p',
                                                     'mmu-miR-137-3p',
                                                     'mmu-miR-100-5p',
                                                     'mmu-miR-192-5p'))

####################################
## Only dataframe selected columns
microbiota_relAbund <- subset(microbiota_by_taxon, select=-c(Grouping))
# miRNAs_data <- miRNAs_data_selected
microbio_behavioral_merged <- merge (microbiota_relAbund, miRNAs_data, by= "mouse_id")
# head(microbio_behavioral_merged)

## All miRNAs
first_miRNA <- 'mmu-miR-34c-3p'; last_miRNA <- 'mmu-miR-7666-3p';
axis_text_size<-24; size_p_values <- 4; angle_reg<-270; microbio_set <- "all"
width_p <- 45; height_p <- 14

## Selected miRNAs
# first_miRNA <- 'mmu-miR-876-5p'; last_miRNA <- 'mmu-miR-192-5p';
# axis_text_size<-16; size_p_values <- 5; angle_reg <- 0; microbio_set <- "selected"
# width_p <- 20; height_p <- 12

## Genus needs smaller labels
# axis_text_size<-8

## Variables
data <- gather(microbio_behavioral_merged, taxon, microbio_rel_ab, first_taxon:last_taxon)%>%
        gather(miRNA, value, first_miRNA:last_miRNA)

data_nest <- group_by(data, taxon, miRNA) %>% nest()


# data_nest <- group_by(data, city, telecon) %>% nest()
# cor_method <- "pearson"
cor_method <- "spearman"
cor_fun <- function(df) cor.test(df$microbio_rel_ab, df$value, method = cor_method) %>% tidy()

library(broom) # Convert results of statistical functions (lm, t.test, cor.test, etc.) into tidy tables

# library(fs)
# library(lubridate)

data_nest <- mutate(data_nest, model = map(data, cor_fun))
# data_nest

# str(slice(data_nest, 1))

corr_pr <- select(data_nest, -data) %>% unnest()
corr_pr <- mutate(corr_pr, sig = ifelse(p.value <0.05, "Sig.", "Non Sig."))
title_p <- paste("Correlations between miRNA expression and", first_up(taxon), "relative abundances\n")
axis_text_size<-18
hm <- ggplot() + geom_tile(data = corr_pr,
                           # aes(taxon, miRNA, fill = estimate),
                           aes(miRNA, taxon, fill = estimate),
                           size = 1,
                           colour = "white") +
      ## black lines around tiles
      geom_tile(data = filter(corr_pr, sig == "Sig."),
                # aes(taxon, miRNA),
                aes(miRNA, taxon),
                size = 1,
                colour = "black",
                fill = "transparent") +
      # geom_text(data = corr_pr, 
      #           angle = angle_reg,
      #           size = size_p_values,
      #           aes(miRNA, taxon,    
      #               label = ifelse(sig == "Sig.", round(estimate, 2),""))) +
                    # fontface = ifelse(sig == "Sig.", "bold", "plain"))) +
      scale_fill_gradient2(breaks = seq(-1, 1, 0.2), #midpoint = mid, 
                           low = "#d53e4f", mid = "white",
                           high = "#abdda4") + 
      labs(y = "", x = "", fill = "", p.value = "", title = title_p) +
      theme_minimal()+
      theme(panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_text(size=axis_text_size),
            plot.title = element_text(size=24, hjust = 0.5),
            axis.text.x = element_text(angle=90),
            legend.text = element_text( size=14))

out_dir <- paste0(home_dir, "/git/food_addiction_analysis/figures/cross_talking_microbio_miRNAs/")
dpi_q <- 200
extension_img <- ".png"

ggsave (hm, file=paste0(out_dir, "heatmap_", microbio_set ,"_microbio_", taxon, extension_img), 
        width = width_p, height = height_p, dpi=dpi_q)
