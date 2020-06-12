#!/usr/bin/env Rscript

#############################################################################
### Jose Espinosa-Carrasco. CB-CRG Group. June 2020                       ###
#############################################################################
### PCA reinstatement experiment from Rafael's lab                        ###
### inactive presses during time out                                      ### 
### This version was done when all the variables were not annotated and   ###
### it is the first version that includes impulsivity/compulsivity        ###
#############################################################################

## https://github.com/cbcrg/phecomp/blob/20dcf868673836f62bd28c1062ba19986114cbf3/lib/R/pca_reinstament_impulsivity.R

#######################
## Example commands
# ./pca_addiction_humans.R --path_tbl_pca="/home/kadomu/projects/20200421_pca_behavior_elena/data/human_questionary/pca_humans.csv" --path_tbl_annotations="/home/kadomu/git/food_addiction_analysis/data/human_questionary/annot_descriptors_pca_human.csv" --plots

# ./pca_addiction_humans.R --path_tbl_pca="/home/kadomu/projects/20200421_pca_behavior_elena/data/human_questionary/PCA_CD1.csv" --path_tbl_annotations="/home/kadomu/git/food_addiction_analysis/data/human_questionary/annot_descriptors_pca_cd1.csv" --plots

# ./pca_addiction_humans.R --path_tbl_pca="/home/kadomu/projects/20200421_pca_behavior_elena/data/human_questionary/PCA_C57.csv" --path_tbl_annotations="/home/kadomu/git/food_addiction_analysis/data/human_questionary/annot_descriptors_pca_c57.csv" --plots

#####################
## VARIABLES
## Reading arguments
args <- commandArgs (TRUE)
# write (paste0("...........................", getwd()), stdout())

## Default setting when no arguments passed
# if ( length(args) < 1) {
#   args <- c("--help")
# }

## Help section
if("--help" %in% args) {
  cat("
        pca_addiction
        Arguments:
        --path_tbl_pca=path          - character, path to read tbl files
        --tbl_annotations=path       - character, tbl containing annotations
        --plots                      - boolean
        --image_format=image_format  - character
        --help                       - print this text
        Example:
        ./pca_addiction.R --path_tbl_pca=\"tbl.txt\" --tbl_annotations=\"tbl_annotations.txt\" --plots\n")
  
  q (save="no")
}

## Use to parse arguments beginning by --
parseArgs <- function(x)
{
  strsplit (sub ("^--", "", x), "=")
}

## Parsing arguments
argsDF <- as.data.frame (do.call("rbind", parseArgs(args)))
argsL <- as.list (as.character(argsDF$V2))
names (argsL) <- argsDF$V1

# path to variables bedgraph files
{
  if (is.null (argsL$path_tbl_pca))
  {
    write ("[Warning]: Path to table containing PCA table set to default", stderr())
    # path_tbl_pca <- '20200517_pca_behavioral_test.csv' 
    stop ("[FATAL]: Path to table containing data for the pca is mandatory")
  }
  else
  {
    path_tbl_pca <- argsL$path_tbl_pca
  }
}

# path to table containing annotations
{
  if (is.null (argsL$path_tbl_annotations))
  {
    write ("[Warning]: Path to table containing annotations for plots set to default", stderr())
    path_tbl_annotations <- "annot_descriptors_24_04_20.csv"
    stop ("[FATAL]: Path to annotations table is mandatory")
  }
  else
  {
    path_tbl_annotations <- argsL$path_tbl_annotations
  }
}
{
  if (is.null (argsL$plots))
  {
    save_plot = FALSE
  }
  else
  {
    save_plot = TRUE
  }
}

#####################
#####################
## PCA of reinstatement matrix
install.packages("FactoMineR")
# Calling libraries
library(Hmisc) # arrow function
# library(calibrate)
# library(multcomp)
library(ggplot2)
library(FactoMineR)

##Getting HOME directory 
home <- Sys.getenv("HOME") 

path_assets <- paste0(home, "/git/food_addiction_analysis")
out_folder <- paste0(path_assets, "/figures/var_behavioral_test/")

# Loading graph params:
source (paste0(path_assets,"/scripts/r/graph_parameters.R"))
# Loading functions:
source (paste0(path_assets,"/scripts/r/functions.R"))

# Parameter to set plot qualities
# save_plot = FALSE
# save_plot = TRUE

v_colours <- c("red", "gray", "blue", "lightblue", "magenta", "orange", "darkgreen")
cb_palette_adapt <- rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000"), 3)
cb_palette_mice <- c("#E69F00","#D55E00", "#56B4E9", "#009E73")

# Parameter to set plot qualities
dpi_q <- 100
size_text_circle <- 5.5
title_var_loadings <- "\nVariables factor map\n" 

path_assets <- paste0(home, "/git/food_addiction_analysis")

# out_folder <- paste0(pwd, "./figures/") 
out_folder <- paste0(path_assets, "/figures/human_questionary/")

name_out <- "pca_c57"
extension_img <- ".png"
# file_path <- path_tbl_pca <- "/home/kadomu/projects/20200421_pca_behavior_elena/data/human_questionary/pca_humans.csv"
file_path <- path_tbl_pca
data <- read.csv (file_path,
                  dec=".", 
                  sep=",",
                  na.strings = "NA",
                  stringsAsFactors = F)
# head (data)
# Remove categorical
data_filtered <-subset (data, select=-c(Grupo))
head (data_filtered)
res = PCA (data_filtered, scale.unit=TRUE)

head (data)

# I set as vector with colors for all the plots I just have to set the number of colours that I need for the plots
v_colours <- c("red", "gray", "blue", "lightblue", "magenta", "orange", "darkgreen")

data$Grupo[is.na(data$Grupo)]  <- "N"
group_mouse <- data$Grupo
title_2_plot <- "PCA humans"
sample_id <- c(1:length(data$Grupo))

list_var <- list()
list_var <- get_var_pcs(res)
pca_r <- res
list_variance <- list_var
# min_x=NA
# max_x=NA
# min_y=NA
# max_y=NA

pca_addiction_PC1_PC2 <- pca_plot (res, pc_x="Dim.1", pc_y="Dim.2", group_v=group_mouse, 
                                   title_pca=title_2_plot, id_v = sample_id,
                                   list_variance = list_var)

pca_addiction_PC1_PC3 <- pca_plot (res, pc_x="Dim.1", pc_y="Dim.3", group_v=group_mouse, 
                                   title_pca=title_2_plot, id_v = sample_id,
                                   list_variance = list_var)

pca_addiction_PC2_PC3 <- pca_plot (res, pc_x="Dim.2", pc_y="Dim.3", group_v=group_mouse, 
                                   title_pca=title_2_plot, id_v = sample_id,
                                   list_variance = list_var)
if (save_plot) {
  ggsave (pca_addiction_PC1_PC2, 
          file=paste(out_folder, name_out, "_PCA_PC1_vs_PC2", extension_img, sep=""), 
          width = 10, height = 10, dpi=dpi_q)
} else {
  pca_addiction_PC1_PC2
}

if (save_plot) {
  ggsave (pca_addiction_PC1_PC3, 
          file=paste(out_folder, name_out, "_PCA_PC1_vs_PC3", extension_img, sep=""), 
          width = 10, height = 10, dpi=dpi_q)
} else {
  pca_addiction_PC1_PC3
}

if (save_plot) {
  ggsave (pca_addiction_PC2_PC3, 
          file=paste(out_folder, name_out, "_PCA_PC2_vs_PC3", extension_img, sep=""), 
          width = 10, height = 10, dpi=dpi_q)
} else {
  pca_addiction_PC2_PC3
}

bars_plot_pc1 <- pca_barPlot(res$var$coord, sel_pc="Dim.1")
bars_plot_pc2 <- pca_barPlot(res$var$coord, sel_pc="Dim.2")
bars_plot_pc3 <- pca_barPlot(res$var$coord, sel_pc="Dim.3")

if (save_plot) {
  ggsave (bars_plot_pc1,
          file=paste(out_folder, name_out, "_bars_PC1", extension_img, sep=""), 
          width = 15, height = 12, dpi=dpi_q)
} else {
  bars_plot_pc1
}

if (save_plot) {
  ggsave (bars_plot_pc2,
          file=paste(out_folder, name_out, "_bars_PC2", extension_img, sep=""), 
          width = 15, height = 12, dpi=dpi_q)
} else {
  bars_plot_pc2
}

if (save_plot) {
  ggsave (bars_plot_pc3,
          file=paste(out_folder, name_out, "_bars_PC3", extension_img, sep=""), 
          width = 15, height = 12, dpi=dpi_q)
} else {
  bars_plot_pc3
}  
  
### Circle Plot

# reinst_annotation
# path_tbl_annotations <- paste0(path_assets, "/data/human_questionary/", "annot_descriptors_24_04_20.csv")

reinst_annotation <- read.csv (path_tbl_annotations,
                               dec=",",
                               sep=",",
                               stringsAsFactors=FALSE)

circle_plot <- as.data.frame (res$var$coord)
circle_plot$var <- rownames (circle_plot)
# reinst_annotation$Label_variable

# merging with annotation tbl
circle_plot_annotation_merged <- merge (circle_plot, reinst_annotation, by.x= "var", by.y = "Label_variable")

circle_plot_pc1_pc2 <- circle_plot_by_gr(circle_plot_annotation_merged, list_variance = list_var)
circle_plot_pc1_pc3 <- circle_plot_by_gr(circle_plot_annotation_merged, pc_x="Dim.1", pc_y="Dim.3", list_variance = list_var)
circle_plot_pc2_pc3 <- circle_plot_by_gr(circle_plot_annotation_merged, pc_x="Dim.2", pc_y="Dim.3", list_variance = list_var)

if (save_plot) {
  ggsave (circle_plot_pc1_pc2,
          file=paste(out_folder, name_out, "_circle_PC1_vs_PC2", extension_img, sep=""), 
          width = 15, height = 15, dpi=dpi_q)
} else {
  circle_plot_pc1_pc2
}

if (save_plot) {
  ggsave (circle_plot_pc2_pc3,
          file=paste(out_folder, name_out, "_circle_PC2_vs_PC3", extension_img, sep=""), 
          width = 15, height = 15, dpi=dpi_q)
} else {
  circle_plot_pc2_pc3
}

if (save_plot) {
  ggsave (circle_plot_pc1_pc3,
          file=paste(out_folder, name_out, "_circle_PC1_vs_PC3", extension_img, sep=""), 
          width = 15, height = 15, dpi=dpi_q)
} else {
  circle_plot_pc1_pc3
}


