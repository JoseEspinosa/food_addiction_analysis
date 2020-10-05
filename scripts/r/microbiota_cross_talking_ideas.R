#!/usr/bin/env Rscript

#############################################################################
### Jose Espinosa-Carrasco CB-CRG Group. Sep 2020                         ###
#############################################################################
### Cross-talking between behavior and microbiota                         ###
#############################################################################

# https://www.nature.com/articles/s41598-018-38218-7

# Make a heatmap and clustering of addicted and no addicted mouse

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

### microbiota data
taxon <- "phylum"; sep_f="\t"
# taxon <- "family"; sep_f=";"
# taxon <- "genus"; sep_f=";"

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

microbiota_by_taxon[,1]
microbiota_relAbund <- subset(microbiota_by_taxon, select=-c(Grouping))
head(microbiota_relAbund )
library ("corrplot")
com = microbiota_by_taxon[,4:length(microbiota_by_taxon[1,])]
cc = cor(com, method = "spearman")
corrplot(cc)
install.packages("heatmaply")
library(heatmaply)
heatmaply_cor(
  cor(com),
  xlab = "Features", 
  ylab = "Features",
  k_col = 2, 
  k_row = 2
)

addicts <- subset(microbiota_by_taxon, Grouping=="Addict")[,4:length(microbiota_by_taxon[1,])]
no_addicts <- subset(microbiota_by_taxon, Grouping=="Non-Addict")[,4:length(microbiota_by_taxon[1,])]

heatmaply_cor(
  cor(addicts),
  xlab = "Features", 
  ylab = "Features",
  k_col = 2, 
  k_row = 2
)

heatmaply_cor(
  cor(no_addicts),
  xlab = "Features", 
  ylab = "Features",
  k_col = 2, 
  k_row = 2
)

# How predictive are variants --> how predictive is micribiome
# Polygenic risk score
# Correlacionar PC de los comportamientos con una pCA de la microbiota

