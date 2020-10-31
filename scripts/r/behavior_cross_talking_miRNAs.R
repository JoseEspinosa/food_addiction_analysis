#!/usr/bin/env Rscript

#############################################################################
### Jose Espinosa-Carrasco CB-CRG Group. Sep 2020                         ###
#############################################################################
### Cross-talking between behavior and miRNAs                             ###
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

home_dir <- Sys.getenv("HOME")

#########################################
## Only behavioral variables used in PCA
# behavioral_data_path <- paste0(home_dir,
#                                "/git/food_addiction_analysis/data/microbiota/behavioral_data_from_results_microbiota_16_04_20_PCA_variables.csv")
# behavioral_data <- read.csv(behavioral_data_path,
#                             dec=",",
#                             sep=";",
#                             check.names = F,
#                             stringsAsFactors = F)
# first_var <-"Persistence_LP"; last_var <- "Aversive_LP"

## Only PCA variables used by Alejandra
behavioral_data_path <- paste0(home_dir, "/git/food_addiction_analysis/data/behavior/PCA_behavioral_data_fromMatrixAlejandra.csv")
behavioral_data <- read.csv(behavioral_data_path,
                            dec=",",
                            sep=";",
                            check.names = F,
                            stringsAsFactors = F)
# head(behavioral_data)
first_var <-"Persistence_Index"; last_var <- "PFPeriod_First_5_min"

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

miRNAs_data_to_transp <- merge (miRNAs_data_replica, miRNAs_data_discovery, by="miRNA")

miRNAs_data <- transpose_df(miRNAs_data_to_transp)
miRNAs_data$mouse_id <- row.names(miRNAs_data)

## Select miRNAs from Elena's table
# miRNAs_data_selected <- subset(miRNAs_data, select=c('mouse_id',
#                                                      'mmu-miR-876-5p',
#                                                      'mmu-miR-211-5p',
#                                                      'mmu-miR-3085-3p',
#                                                      'mmu-miR-665-3p',
#                                                      'mmu-miR-3072-3p',
#                                                      'mmu-miR-124-3p',
#                                                      'mmu-miR-29c-3p',
#                                                      'mmu-miR-544-3p',
#                                                      'mmu-miR-137-3p',
#                                                      'mmu-miR-100-5p',
#                                                      'mmu-miR-192-5p'))
###################
## All miRNAs
# miRNAs_data_selected <-miRNAs_data

###################
## Merging of behavioral data + miRNAs 
miRNAs_behavioral_merged <- merge (miRNAs_data_selected, behavioral_data, by= "mouse_id")
# head(miRNAs_behavioral_merged)

## All miRNAs
first_miRNA <- 'bta-miR-2478'; last_miRNA <- 'xtr-miR-9b-5p';
axis_text_size<-24; size_p_values <- 4; angle_reg<-270; microbio_set <- "all"
width_p <- 45; height_p <- 14

## Selected miRNAs
# first_miRNA <- 'mmu-miR-876-5p'; last_miRNA <- 'mmu-miR-192-5p'
# axis_text_size<-16; size_p_values <- 5; angle_reg <- 0; microbio_set <- "selected"
# width_p <- 20; height_p <- 12

## My PCA var
# data <- gather(miRNAs_behavioral_merged, miRNA, value, first_miRNA:last_miRNA)%>%
#         gather(behavior_idx, index, Persistence_LP:Aversive_LP)
## PCA var Alejandra
data <- gather(miRNAs_behavioral_merged, miRNA, value, first_miRNA:last_miRNA)%>%
        gather(behavior_idx, index, first_var:last_var)
data_nest <- group_by(data, miRNA, behavior_idx) %>% nest()

# data_nest <- group_by(data, city, telecon) %>% nest()
# cor_method <- "pearson"
cor_method <- "spearman"
cor_fun <- function(df) cor.test(df$value, df$index, method = cor_method) %>% tidy()

library(broom) # Convert results of statistical functions (lm, t.test, cor.test, etc.) into tidy tables

# library(fs)
# library(lubridate)

data_nest <- mutate(data_nest, model = map(data, cor_fun))

data_nest

# str(slice(data_nest, 1))

corr_pr <- select(data_nest, -data) %>% unnest()
corr_pr <- mutate(corr_pr, sig = ifelse(p.value <0.05, "Sig.", "Non Sig."))

###########
## FDR
fdr_cutoff <- 0.2
test_p <- corr_pr$p.value
# class (corr_pr$p.value)

# test_p<-c(test_p,0.001459265)
p.adjust(test_p, method = 'BH', n = length(test_p))
corr_pr$fdr <- p.adjust(corr_pr$p.value, method = 'BH', n = length(corr_pr$p.value))
corr_pr <- mutate(corr_pr, sig = ifelse(fdr < fdr_cutoff, "Sig.", "Non Sig."))
corr_pr$fdr

hm <- ggplot() + geom_tile(data = corr_pr,
                           # aes(behavior_idx, miRNA, fill = estimate),
                           aes(miRNA, behavior_idx, fill = estimate),
                               size = 1,
                               colour = "white") +
      ## black lines around tiles
      geom_tile(data = filter(corr_pr, sig == "Sig."),
                # aes(behavior_idx, miRNA),
                aes(miRNA, behavior_idx),
                size = 1,
                colour = "black",
                fill = "transparent") +
      # geom_text(data = corr_pr, 
      #           angle = angle_reg,
      #           size=5,
      #           # aes(behavior_idx, miRNA, #label = round(estimate, 2),
      #           aes(miRNA, behavior_idx,    
      #               label = ifelse(sig == "Sig.", round(estimate, 2), ""))) +
      #               # fontface = ifelse(sig == "Sig.", "bold", "plain"))) +
      #     # scale_fill_gradient2(breaks = seq(-1, 1, 0.2)) +
      scale_fill_gradient2(breaks = seq(-1, 1, 0.2), #midpoint = mid, 
                           low = "#d53e4f", mid = "white",
                           high = "#abdda4") + 
      labs(y = "", x = "", fill = "", p.value = "") +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_text(size=axis_text_size),
            axis.text.x = element_text(angle=90),
            legend.text = element_text( size=14))

hm
# out_dir <- paste0(home_dir, "/git/food_addiction_analysis/figures/cross_talking_behavior_miRNAs/")
out_dir <- paste0(home_dir, "/git/food_addiction_analysis/figures/cross_talking_behavior_miRNAs_selected/")
dpi_q <- 200
extension_img <- ".png"
ggsave (hm, file=paste0(out_dir, "heatmap_fewVar_", microbio_set ,
                        "miRNAs", extension_img), 
        width=width_p, height = height_p, dpi=dpi_q)

