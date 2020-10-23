
#############################################################################
### Jose Espinosa-Carrasco CB-CRG Group. Sep 2020                         ###
#############################################################################
### Cross-talking between behavior and microbiota                         ###
#############################################################################

### Based on
## https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-0710-2#Sec2

library(ggplot2)
library(tidyverse)
library(dplyr)
library(broom) # Convert results of statistical functions (lm, t.test, cor.test, etc.) into tidy tables
library(Hotelling) # center log ratio transformation

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
taxon <- "genus"; sep_f=";"; first_taxon<-'Acetatifactor'; last_taxon<-'Tyzzerella'; axis_text_size_y <- 10

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

########################
## Filter 
## Half of the animals has at least 0.1 relative abundance
# microbiota_by_taxon <- microbiota_by_taxon [, -1]
# microbiota_by_taxon_filt <- microbiota_by_taxon [ , colSums(microbiota_by_taxon > 0.1) >= 6 ]

# microbiota_by_taxon[,c(-1,-2)] %>%
#   group_by(Grouping) %>% 
#   filter(n()>50, colSums(valuation==0) <10)

min_n_samples <- microbiota_by_taxon[,c(-1,-2)]  %>%
                  group_by(Grouping) %>%
                  summarise_all(funs(sum(.!=0))) %>%
                  summarise_if(is.numeric, min) %>%
                  as.data.frame()

min_n_samples <- floor(min_n_samples/2)

n_0.1_rel_ab <- microbiota_by_taxon[,c(-1,-2)]  %>%
  group_by(Grouping) %>%
  summarise_all(funs(sum(.>0.1)))%>%
  summarise_if(is.numeric, min) %>%
  as.data.frame()

df_to_filter <- dplyr::bind_rows(min_n_samples, n_0.1_rel_ab)

v <- df_to_filter[2,] - df_to_filter[1,] 
genus_to_keep <- colnames(v[which(v >= 1)])

microbiota_by_taxon_filt <- subset(microbiota_by_taxon, select=genus_to_keep)

microbiota_by_taxon_filt_pseudoCts <- microbiota_by_taxon_filt

## Adding pseudo-counts from here # https://genominfo.org/journal/view.php?number=549
microbiota_by_taxon_filt_pseudoCts[microbiota_by_taxon_filt_pseudoCts == 0] <- 0.001

## Esta bien by rows porque cada row es un sample, un raton
microbiota_by_taxon_filt_transf <- cbind (microbiota_by_taxon[,c(1:3)], 
                                          Hotelling::clr(microbiota_by_taxon_filt_pseudoCts[,c(-1,-2,-3)]))

### All animals
## microbiota_by_taxon <- microbiota_by_taxon; 
suffix <- "_all"; title_tag <- "All individuals"

## Only addicts
# microbiota_by_taxon <- subset(microbiota_by_taxon, Grouping=="Addict"); suffix <- "_addict"; title_tag <- "Addict individuals"

## Only no addicts
# microbiota_by_taxon <- subset(microbiota_by_taxon, Grouping=="Non-Addict"); suffix <- "_no_addict"; title_tag <- "No addict individuals"

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
## show rows with zeros
# dd <-miRNAs_data_to_transp
# filter_all(dd, any_vars(. == 0))
# miRNAs_data_to_transp[5,]
# Hotelling::clr(miRNAs_data_to_transp[5,-1])

miRNAs_data <- transpose_df(miRNAs_data_to_transp)
# miRNAs_data_pseudoCts <- miRNAs_data
# miRNAs_data_transf <- Hotelling::clr(as.matrix(miRNAs_data_pseudoCts[1,]))
# 
# exp(mean(log(as.matrix(miRNAs_data_pseudoCts[1,]))))
# 
# v<- as.vector(log(miRNAs_data_pseudoCts[1,]))
# mean(v)
# class(v)
# exp(mean(log(miRNAs_data_pseudoCts[1,])))
# # miRNAs_data_pseudoCts[miRNAs_data_pseudoCts == 0] <- 0.001
# 
# ## Esta bien by rows porque cada row es un sample, un raton
# miRNAs_data_transf <- Hotelling::clr(miRNAs_data_pseudoCts)
# # class(miRNAs_data_pseudoCts[3,5])

# miRNAs_data_transf$mouse_id <- row.names(miRNAs_data)
miRNAs_data$mouse_id <- row.names(miRNAs_data)

#################
# Selected miRNAs
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
# Hotelling::clr(miRNAs_data_selected)

miRNAs_data <- miRNAs_data_selected

first_miRNA <- 'mmu-miR-876-5p'; last_miRNA <- 'mmu-miR-192-5p';
axis_text_size_x <- 16; size_p_values <- 5; angle_reg <- 0; microbio_set <- "selected"
width_p <- 20; height_p <- 12

##################
# #################
# # All miRNAs
# first_miRNA <- 'mmu-miR-34c-3p'; last_miRNA <- 'mmu-miR-7666-3p';
# axis_text_size_x <- 8; size_p_values <- 4; angle_reg<-270; microbio_set <- "all"
# width_p <- 45; height_p <- 14

###############
microbiota_relAbund <- subset(microbiota_by_taxon_filt_transf,
                              select=-c(Grouping))

microbio_behavioral_merged <- merge (microbiota_relAbund,
                                     miRNAs_data,
                                     by= "mouse_id")
first_taxon<-"Alistipes"
last_taxon <- "Tyzzerella"
first_miRNA <- "bta-miR-2478"; 
last_miRNA <- "xtr-miR-9b-5p"

# last
first_taxon<-"Anaerotruncus"
last_taxon <- "Ruminococcaceae_UCG014"
first_miRNA <- "mmu-miR-876-5p"; 
last_miRNA <- "mmu-miR-192-5p"

## Test con selected
# first_miRNA <- "bta-miR-2478"; 
# last_miRNA <- "mmu-let-7e-3p"

## Variables
data <- gather(microbio_behavioral_merged, taxon, microbio_rel_ab, first_taxon:last_taxon)%>%
        gather(miRNA, value, first_miRNA:last_miRNA)

data_nest <- group_by(data, taxon, miRNA) %>% nest()


# data_nest <- group_by(data, city, telecon) %>% nest()
# cor_method <- "pearson"
cor_method <- "spearman"
cor_fun <- function(df) cor.test(df$microbio_rel_ab, df$value, method = cor_method) %>% tidy()


# library(fs)
# library(lubridate)

data_nest <- mutate(data_nest, model = map(data, cor_fun))
data_nest

# str(slice(data_nest, 1))
dim(data_nest[[3]][[1]]$value)

corr_pr <- select(data_nest, -data) %>% unnest()
corr_pr <- mutate(corr_pr, sig = ifelse(p.value <0.05, "Sig.", "Non Sig."))

###########
## FDR
fdr_cutoff <- 0.2
test_p <- corr_pr$p.value 
min(test_p)
# class (corr_pr$p.value)
head(test_p)
test_p<- test_p[order(test_p)]
test_p <- c(0.00001459265, test_p )
qvalue(test_p)
class (corr_pr$p.value)

min(p.adjust(test_p, method = 'BH', n = length(test_p)))
plot(qvalue(corr_pr$p.value)$qvalues)

corr_pr$fdr <- p.adjust(corr_pr$p.value, method = 'BH', n = length(corr_pr$p.value))
corr_pr <- mutate(corr_pr, sig = ifelse(fdr < fdr_cutoff, "Sig.", "Non Sig."))
corr_pr$fdr

## The result of p.adjust and qvalue is the same
# BiocManager::install("qvalue")
# library(qvalue)
## Plot
q_value <- qvalue(corr_pr$p.value)
plot(q_value)
corr_pr$qvalue <- q_value$qvalues
corr_pr <- mutate(corr_pr, sig = ifelse(qvalue <0.2, "Sig.", "Non Sig."))
corr_pr$qvalue

#######
sign <- subset(corr_pr, fdr<0.2)
sign
title_tag <- "log transf"
title_p <- paste("Correlations between miRNA expression and", 
                 first_up(taxon), "relative abundances\n", title_tag)
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
        axis.text.x = element_text(size=axis_text_size_x, angle=90),
        axis.text.y = element_text(size=axis_text_size_y),
        plot.title = element_text(size=24, hjust = 0.5),
        legend.text = element_text( size=14))
hm

out_dir <- paste0(home_dir, "/git/food_addiction_analysis/figures/cross_talking_microbio_miRNAs/")
dpi_q <- 200
extension_img <- ".png"
suffix <- "testttttt"
microbio_set <- "sssss"
ggsave (hm, file=paste0(out_dir, "heatmap_", microbio_set ,"_microbio_logTransform", taxon, suffix, extension_img), 
        width = width_p, height = height_p, dpi=dpi_q)


## Is done by row
data(bottle.df)

clr(bottle.df, "Number")

clr(bottle.df[1,],"Number")

# t(sapply(split(microbiota_by_taxon_filt[,c(-1,-2,-3)], c("Addict","Non-Addict")), function (x) colSums(x)>0.1)) 
# t(sapply(split(microbiota_by_taxon_filt[,c(-1,-2,-3)], c("Addict","Non-Addict")), function (x) print(x)))
# function(x) length(x[x<0])