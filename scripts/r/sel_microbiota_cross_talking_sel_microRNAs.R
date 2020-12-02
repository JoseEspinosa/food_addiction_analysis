
#############################################################################
### Jose Espinosa-Carrasco CB-CRG Group. Noc 2020                         ###
#############################################################################
### Cross-talking between miRNAs and microbiota                           ###
#############################################################################

### Based on
## https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-0710-2#Sec2

library(ggplot2)
library(tidyverse)
library(dplyr)
library(broom) # Convert results of statistical functions (lm, t.test, cor.test, etc.) into tidy tables
library(Hotelling) # center log ratio transformation
library(microbiome)
library(qvalue)

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
taxon <- "genus"; sep_f=";"; first_taxon<-'Acetatifactor'; last_taxon<-'Tyzzerella'

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

## Hardcoded taxons to keep
taxon_to_keep <- c("Enterorhabdus",
                   "Gastranaerophilales_uncultured bacterium",
                   "Anaeroplasma",
                   "Blautia",
                   "Lachnospiraceae_UCG001",
                   "Allobaculum")

microbiota_by_taxon_filt <- subset(microbiota_by_taxon, select=taxon_to_keep)

# microbiota_by_taxon_filt
microbiota_by_taxon_filt_transp <- as.data.frame(t(microbiota_by_taxon_filt))

## Clr transform
clr_microbiomePack_transform <- microbiome::transform(microbiota_by_taxon_filt_transp, "clr")
clr_microbiomePack_transf_toBind <- as.data.frame(t(clr_microbiomePack_transform))
microbiota_by_taxon_filt_clr_microbiome <- cbind (microbiota_by_taxon[,c(2:3)], 
                                                  clr_microbiomePack_transf_toBind)

### All animals
## microbiota_by_taxon <- microbiota_by_taxon; 
suffix <- "_all"; title_tag <- "All individuals (addicts and no-addicts)"

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

## mice (samples in columns for microbiome package), before transpose
## Important to use cbind.data.frame to avoid conversion to character!!!
# miRNAs_data_to_transp_transform <-  as.data.frame(cbind.data.frame(miRNAs_data_to_transp[,1], microbiome::transform(miRNAs_data_to_transp[,-1], "clr")))
# miRNAs_data <- transpose_df(miRNAs_data_to_transp_transform)
# miRNAs_data$mouse_id <- row.names(miRNAs_data)

#################
# Selected miRNAs from transformed data
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

# first_miRNA <- 'mmu-miR-876-5p'; last_miRNA <- 'mmu-miR-192-5p';
# axis_text_size_x <- 16; size_p_values <- 5; angle_reg <- 0; microbio_set <- "selected"
# width_p <- 20; height_p <- 12

#################
# Selected miRNAs from not transformed data
## Select miRNAs from Elena's table
miRNAs_data <- miRNAs_data_to_transp
miRNAs_data_selected_to_transform <- miRNAs_data[miRNAs_data$miRNA %in% c('mouse_id',
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
                                                                          'mmu-miR-192-5p'), ]

miRNAs_data_to_transp_transform <-  as.data.frame(cbind.data.frame(miRNAs_data_selected_to_transform[,1], 
                                                  microbiome::transform(miRNAs_data_selected_to_transform[,-1], "clr")))
miRNAs_data <- transpose_df(miRNAs_data_to_transp_transform)
miRNAs_data$mouse_id <- row.names(miRNAs_data)
miRNAs_data_selected <- miRNAs_data

first_miRNA <- 'mmu-miR-100-5p'; last_miRNA <- 'mmu-miR-876-5p';

###############
# microbiota_relAbund <- subset(microbiota_by_taxon_filt_transf,
#                               select=-c(Grouping))
microbiota_relAbund <- subset(microbiota_by_taxon_filt_clr_microbiome, 
                              select=-c(Grouping))

## selected miRNAs
microbio_behavioral_merged <- merge (microbiota_relAbund,
                                     miRNAs_data_selected,
                                     by= "mouse_id")
first_taxon <- taxon_to_keep[1]
last_taxon <- tail(taxon_to_keep, n=1)
first_miRNA <- "mmu-miR-100-5p";
last_miRNA <- "mmu-miR-876-5p"

## Variables
data <- gather(microbio_behavioral_merged, taxon, microbio_rel_ab, first_taxon:last_taxon)%>%
        gather(miRNA, value, first_miRNA:last_miRNA)

data_nest <- group_by(data, taxon, miRNA) %>% nest()

# cor_method <- "pearson"
cor_method <- "spearman"
cor_fun <- function(df) cor.test(df$microbio_rel_ab, df$value, method = cor_method) %>% tidy()

data_nest <- mutate(data_nest, model = map(data, cor_fun))

corr_pr <- select(data_nest, -data) %>% unnest()
corr_pr <- mutate(corr_pr, sig = ifelse(p.value <0.05, "Sig.", "Non Sig."))

###########
## FDR
# fdr_cutoff <- 0.3
# test_p <- corr_pr$p.value
# test_p <- test_p[order(test_p)]
# 
# corr_pr$fdr <- p.adjust(corr_pr$p.value, method = 'BH', n = length(corr_pr$p.value))
# corr_pr <- mutate(corr_pr, sig = ifelse(fdr < fdr_cutoff, "Sig.", "Non Sig."))
# 
# q_value <- qvalue(corr_pr$p.value)
# plot(q_value)
# corr_pr$qvalue <- q_value$qvalues
# corr_pr <- mutate(corr_pr, sig = ifelse(qvalue <0.2, "Sig.", "Non Sig."))

# sign <- subset(corr_pr, fdr<0.3)
# sign

#######
## Plot
title_tag <- "log transf"
title_p <- paste("Correlations between miRNA expression and", 
                 first_up(taxon), "relative abundances\n", title_tag)
title_p <- ""

axis_text_size_x <- 18; size_p_values <- 5; angle_reg <- 0; microbio_set <- "selected"
axis_text_size_y <- 18
width_p <- 20; height_p <- 12

# axis_text_size<-18
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
      scale_fill_gradient2(breaks = seq(-1, 1, 0.2), #midpoint = mid, 
                           low = "#d53e4f", mid = "white",
                           high = "#abdda4") + 
      labs(y = "", x = "", fill = "", p.value = "", title = title_p) +
      theme_minimal()+
      theme(panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(size=axis_text_size_x, angle=90),
            axis.text.y = element_text(size=axis_text_size_y),
            plot.title = element_text(size=24, hjust = 0.5),
            legend.text = element_text( size=14))
hm

out_dir <- paste0(home_dir, "/git/food_addiction_analysis/figures/cross_talking_microbio_miRNAs_transformed/")
dpi_q <- 300
width_p <- 20; height_p <- 12
extension_img <- ".png"
ggsave (hm, file=paste0(out_dir, "heatmap_", microbio_set ,"_sel_microbio_logTransform_", "sel_miRNAs_", taxon, suffix, extension_img), 
        width = width_p, height = height_p, dpi=dpi_q)
