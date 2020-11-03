
#############################################################################
### Jose Espinosa-Carrasco CB-CRG Group. Nov 2020                         ###
#############################################################################
### Cross-talking between behavior and microbiome                         ###
#############################################################################

home_dir <- Sys.getenv("HOME")
taxon <- "genus"
ra <- paste0(home_dir, "/git/food_addiction_analysis/forJose/relative_abundances_by_", taxon, ".txt")
ma <- read.table(ra, sep="\t", skip=1,
                 row.names=1, colClasses = "character")
ma
m_microbiome=read.table(ra,sep="\t",skip=3,row.names=1)
m
colnames(m)=ma[1,]
m
cor_m_microbiome <- cor(m, method="pearson")
cor_m
dist_m_microbiome <- 1-cor_m

## Behavior
bf <- paste0(home_dir, "/git/food_addiction_analysis/forJose/behavior.txt")
mb=read.table(bf,header=TRUE,sep="\t",row.names=1)
mb_t <- t(mb)

## reorder to match microbiome columns (mice)
col_names_microb <- colnames(m)
mb_t_ord <- mb_t [, col_names_microb]
cor_m_beh <- cor(mb_t_ord, method="pearson")
dist_m_beh <- 1 - cor_m_beh

colnames (dist_m_beh)
colnames (dist_m_microbiome)

str(dist_m_beh)
str(dist_m_microbiome)

cor(c(dist_m_beh), c(dist_m_beh))
cor(c(dist_m_beh), c(dist_m_microbiome), method="spearman")






MB=mb[rownames(M2),]

?cor


# install.packages("tidyverse")
library(tidyverse)

# install.packages("corrr")
library(corrr)

m %>%
corrr::correlate() %>%
  corrr::focus(mpg:hp, mirror = TRUE) %>%
  # converts the upper triangle (default) to missing values
  corrr::shave() %>%
  # converts a correlation df into clean matrix
  corrr::fashion()


M=t(m)
rownames(M)=ma[1,]

gr=as.character(ma[2,])


