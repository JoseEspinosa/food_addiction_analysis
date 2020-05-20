#############################################################################
### Jose A Espinosa. NPMMD/CB-CRG Group. Jan 2016                         ###
#############################################################################
### PCA reinstatement experiment from Rafael's lab                        ###
### inactive presses during time out                                      ### 
### This version was done when all the variables were not annotated and   ###
### it is the first version that includes impulsivity/compulsivity        ###
#############################################################################

#####################
#####################
## PCA of behavioral measures

# Calling libraries
library(Hmisc) # arrow function
library(ggplot2) 
library(FactoMineR)

##Getting HOME directory 
home <- Sys.getenv("HOME")

# Loading functions:
source ("./scripts/r/graph_parameters.R")

# Parameter to set plot qualities
# save_plot = FALSE
save_plot = TRUE

dpi_q <- 100
size_text_circle <- 5.5
title_var_loadings <- "\nVariables factor map\n" 
out_folder <- "./figures/" 
out_folder <- "./figures/var_behavioral_test/" 

name_out <- "var_behavioral_test"
extension_img <- ".png"
# file_path <- "./data/tbl_PCA_all_variables.csv"
# file_path <- "./data/tbl_PCA_main_variables.csv"
# file_path <- './data/20200429_tbl_1st_analysis_PCA_selection_variables.csv'
## Only variables used for the behavioral test 
file_path <- './data/20200517_pca_behavioral_test.csv' 

data <- read.csv (file_path,
                  dec=".", 
                  sep=",",
                  na.strings = "NA",
                  stringsAsFactors = F)

# reinst_annotation
reinst_annotation <- read.csv ("./data/annot_descriptors_24_04_20.csv",
                               dec=",",
                               sep=",",
                               stringsAsFactors=FALSE)

## Remove variables that all tables have
data_filtered <- subset (data, select=-c(Mice, 
                                         Genotype, 
                                         Food.pellets, 
                                         Addiction_categorization_LP,
                                         Criteria_EP,
                                         Criteria_MP,
                                         Criteria_LP)) 

data_filtered <- subset (data, select=-c(Mice, 
                                         Genotype, 
                                         Food.pellets, 
                                         Addiction_categorization_LP,
                                         Criteria_EP,
                                         Criteria_MP,
                                         Criteria_LP)) 

# # subset(dt, select = grep("bar|baz", names(dt)))
# # phase <- "LP"
# # phase <- "EP"
# phase <- "MP"
# data_filtered <- subset(data_filtered, select = grep(phase, names(data_filtered)))
# names(data_filtered)
# out_folder <- paste0 ("./figures/var_behavioral_test_", phase, "/") 
# name_out <- paste0 ("var_behavioral_test_", phase)

# Remove categorical
# data_filtered <-subset (data, select=-c(Mice, 
#                                         Genotype,
#                                         Food.pellets, 
#                                         Addiction_categorization_LP,
#                                         Appetitive_extinction_crit_EP, 
#                                         Appetitive_extinction_crit_MP, 
#                                         Appetitive_extinction_crit_LP,
#                                         Appetitive_extinction_crit_EP,
#                                         Appetitive_reinstatement_crit_EP,
#                                         Appetitive_reinstatement_crit_MP,
#                                         Appetitive_reinstatement_crit_LP,
#                                         Criteria_EP,
#                                         Criteria_MP,
#                                         Criteria_LP,
#                                         Addiction_prediction_EP,
#                                         Addiction_prediction_MP,
#                                         Addiction_categorization_LP.1))

# data_filtered <- subset (data, select=-c(Mice, 
#                                          Genotype, 
#                                          Food.pellets, 
#                                          Addiction_categorization_LP, 
#                                          Appetitive_extinction_crit_EP, 
#                                          Appetitive_extinction_crit_MP,
#                                          Appetitive_extinction_crit_LP,
#                                          Appetitive_extinction_crit_EP,
#                                          Appetitive_reinstatement_crit_EP,
#                                          Appetitive_reinstatement_crit_MP,
#                                          Appetitive_reinstatement_crit_LP,
#                                          Addiction_prediction_EP,
#                                          Addiction_prediction_MP,
#                                          Addiction_categorization_LP.1))

## Functions to assign labels and variance
label_pc <- function (name_dim="Dim.1") {
  return (switch(name_dim, 
                 "Dim.1" = "PC1", 
                 "Dim.2" = "PC2", 
                 "Dim.3" = "PC3",
                 "Dim.4" = "PC4",
                 "Dim.5" = "PC5",
                 "Dim.6" = "PC6",
                 "Dim.7" = "PC7"))
}

var_pc <- function (name_dim="Dim.1") {
  return (switch(name_dim, 
                 "Dim.1" = var_PC1, 
                 "Dim.2" = var_PC2,
                 "Dim.3" = var_PC3,
                 "Dim.4" = var_PC4,
                 "Dim.5" = var_PC5,
                 "Dim.6" = var_PC6,
                 "Dim.7" = var_PC7))
}

pca_plot <- function (df_to_pca, pc_x="Dim.1", pc_y="Dim.2", id_v=NA, group_v=NA, title_pca="", 
                      min_x=NA, max_x=NA, min_y=NA, max_y=NA) {
  
  res = PCA (data_filtered, scale.unit=TRUE)
  cb_palette <- c("#E69F00","#D55E00", "#56B4E9", "#009E73", "#999999", "#CC79A7")
  
  list_var <- list()
  dim_name <- "Dim."
  
  # Variance of PCs
  for (i in 1 : length(res$eig [,2])) {
    id = paste0(dim_name, i)
    list_var [[id]] <- round (res$eig [i,2], 1)
    
  }
  
  # PC coordinates are store here
  # Convert PCA results to data frame
  pca2plot <- as.data.frame (res$ind$coord)
  
  # Set axis minima and maxima
  if (is.na(min_x)) {
    min_x <- floor(min(pca2plot[pc_x]))
  }
  if (is.na(max_x)) {
    max_x <- ceiling(max(pca2plot[pc_x]))
  }
  if (is.na(min_y)) {
    min_y <- floor(min(pca2plot[pc_y]))
  }
  if (is.na(max_y)) {
    max_y <- ceiling(max(pca2plot[pc_y]))
  }

  n_samples <- length (pca2plot[,1])
  
  if (!is.na(id_v[1])) {
    pca2plot$id <- group_v
  }
  else {
    pca2plot$id <- 1:n_samples
  }
  
  if (!is.na(group_v[1])) {
    pca2plot$Group <- group_v
  }
  else {
    pca2plot$Group <- "NA"
  }
  
  title_p <- paste(title_pca, label_pc(pc_x), "vs.", label_pc(pc_y), "\n", sep=" ")
  print (title_p)
  title_x <- paste("\n", label_pc(pc_x), "(", list_var[[pc_x]], "% of variance)", sep="")
  title_y <- paste(label_pc(pc_x), " (", list_var[[pc_y]], "% of variance)\n", sep = "")
  
  pca_p <- ggplot (pca2plot, aes(x=get(pc_x), y=get(pc_y), colour=Group)) +
    geom_point (size = 3.5, show.legend = T) +
    scale_color_manual (values=cb_palette) +
    geom_text (aes(label=id), vjust=-0.5, hjust=1, size=4, show.legend = F)+
    theme(legend.key=element_rect(fill=NA)) +
    scale_x_continuous (limits=c(min_x, max_x), breaks=min_x:max_x) +
    scale_y_continuous (limits=c(min_y, max_y), breaks=min_y:max_y) +
    labs(title = title_p, x =  title_x,
         y=title_y) +
    guides(colour = guide_legend(override.aes = list(size = 3)))+
    theme(legend.key=element_rect(fill=NA)) +
    coord_fixed()

  return (pca_p)
}




group_mouse <- paste0(data$Addiction_categorization_LP, "_", data$Criteria_LP)
group_mouse <- factor(group_mouse, levels = c("N_0", "N_1", "A_2", "A_3"))
title_2_plot <- "PCA addiction"
sample_id <- data$Mice

pca_addiction_PC1_PC2 <- pca_plot (data_filtered, group_v=group_mouse, title_pca=title_2_plot, id_v = sample_id)
pca_addiction_PC1_PC3 <- pca_plot (data_filtered, group_v=group_mouse, title_pca=title_2_plot, id_v = sample_id)


pca_plot (data_filtered)



res = PCA (data_filtered, scale.unit=TRUE)

# I set a vector of colors for all the plots I just have to set the number of colours that I need for the plots
v_colours <- c("red", "gray", "blue", "lightblue", "magenta", "orange", "darkgreen")
cb_palette_adapt <- rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000"), 3)

cb_palette_mice <- c("#E69F00","#D55E00", "#56B4E9", "#009E73")

# Variance of PCs
var_PC1 <- round (res$eig [1,2], 1)
var_PC2 <- round (res$eig [2,2], 1)
var_PC3 <- round (res$eig [3,2], 1)
var_PC4 <- round (res$eig [4,2], 1)
var_PC5 <- round (res$eig [5,2], 1)
var_PC6 <- round (res$eig [6,2], 1)
var_PC7 <- round (res$eig [7,2], 1)

# PC coordinates are store here
# Convert PCA results to data frame
pca2plot <- as.data.frame (res$ind$coord)
max(res$ind$coord)
ceiling(min(pca2plot["Dim.1"],pca2plot["Dim.2"]))
floor(min(pca2plot["Dim.1"],pca2plot["Dim.2"]))
pca2plot$id <- data$Mice
# head(data)

# Changes labels of the groups, some mice are categorized as addicted and the rest as NA (non-addicted)
data$Addiction_categorization_LP[is.na(data$Addiction_categorization_LP)]  <- "N"
# pca2plot$group <- data$Addiction_categorization_LP
pca2plot$group <- paste0 (data$Addiction_categorization_LP, "_", data$Criteria_LP)
class(pca2plot$group)
# pca2plot$group[order(pca2plot$group, decreasing=FALSE)]

pca2plot$group <- factor(pca2plot$group, levels = c("N_0", "N_1", "A_2", "A_3"))

## PC1 vs. PC3
title_p <- paste ("PCA addiction, PC1 vs. PC2\n", sep="")

pca_addiction_PC1_PC2 <- ggplot (pca2plot, aes(x=Dim.1, y=Dim.2, colour=group)) + 
                 geom_point (size = 3.5, show.legend = T) + 
                 # scale_color_manual(values=c("orange", "blue")) +
                 # scale_color_manual (values=cb_palette_adapt) +
                 scale_color_manual (values=cb_palette_mice) +
                 geom_text (aes(label=id), vjust=-0.5, hjust=1, size=4, show.legend = F)+
                 theme(legend.key=element_rect(fill=NA)) +
                 scale_x_continuous (limits=c(-10, 12), breaks=-10:12) + 
                 scale_y_continuous (limits=c(-10, 12), breaks=-10:12) +
                 labs(title = title_p, x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
                      y=paste("PC2 (", var_PC2, "% of variance)\n", sep = "")) +
                 guides(colour = guide_legend(override.aes = list(size = 3)))+
                 theme(legend.key=element_rect(fill=NA))

# keeping aspect ratio
pca_addiction_PC1_PC2_aspect_ratio <- pca_addiction_PC1_PC2 + coord_fixed()

if (save_plot) {
  ggsave (pca_addiction_PC1_PC2_aspect_ratio, 
          file=paste(out_folder, name_out, "_PCA_PC1_vs_PC2", extension_img, sep=""), 
          width = 10, height = 10, dpi=dpi_q)
} else {
  pca_addiction_PC1_PC2_aspect_ratio
}

## PC1 vs. PC3
title_p <- paste ("PCA addiction, PC1 vs. PC3\n", sep="")

pca_addiction_PC1_PC3 <- ggplot (pca2plot, aes(x=Dim.1, y=Dim.3, colour=group)) + 
  geom_point (size = 3.5, show.legend = T) + 
  # scale_color_manual(values=c("orange", "blue")) +
  # scale_color_manual (values=cb_palette_adapt) +
  scale_color_manual (values=cb_palette_mice) +
  geom_text (aes(label=id), vjust=-0.5, hjust=1, size=4, show.legend = F)+
  theme(legend.key=element_rect(fill=NA)) +
  scale_x_continuous (limits=c(-10, 12), breaks=-10:12) + 
  scale_y_continuous (limits=c(-10, 12), breaks=-10:12) +
  labs(title = title_p, x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
       y=paste("PC3 (", var_PC3, "% of variance)\n", sep = "")) +
  guides(colour = guide_legend(override.aes = list(size = 3)))+
  theme(legend.key=element_rect(fill=NA))

# keeping aspect ratio
pca_addiction_PC1_PC3_aspect_ratio <- pca_addiction_PC1_PC3 + coord_fixed()

if (save_plot) {
  ggsave (pca_addiction_PC1_PC3_aspect_ratio, 
          file=paste(out_folder, name_out, "_PCA_PC1_vs_PC3", extension_img, sep=""), 
          width = 10, height = 10, dpi=dpi_q)
} else {
  pca_addiction_PC1_PC3_aspect_ratio
}

## PC2 vs. PC3
title_p <- paste ("PCA addiction, PC2 vs. PC3\n", sep="")

pca_addiction_PC2_PC3 <- ggplot (pca2plot, aes(x=Dim.2, y=Dim.3, colour=group)) + 
  geom_point (size = 3.5, show.legend = T) + 
  # scale_color_manual(values=c("orange", "blue")) +
  # scale_color_manual (values=cb_palette_adapt) +
  scale_color_manual (values=cb_palette_mice) +
  geom_text (aes(label=id), vjust=-0.5, hjust=1, size=4, show.legend = F)+
  theme(legend.key=element_rect(fill=NA)) +
  scale_x_continuous (limits=c(-10, 12), breaks=-10:12) + 
  scale_y_continuous (limits=c(-10, 12), breaks=-10:12) +
  labs(title = title_p, x = paste("\nPC2 (", var_PC2, "% of variance)", sep=""), 
       y=paste("PC3 (", var_PC3, "% of variance)\n", sep = "")) +
  guides(colour = guide_legend(override.aes = list(size = 3)))+
  theme(legend.key=element_rect(fill=NA))

# keeping aspect ratio
pca_addiction_PC2_PC3_aspect_ratio <- pca_addiction_PC2_PC3 + coord_fixed()

if (save_plot) {
  ggsave (pca_addiction_PC2_PC3_aspect_ratio, 
          file=paste(out_folder, name_out, "_PCA_PC2_vs_PC3", extension_img, sep=""), 
          width = 10, height = 10, dpi=dpi_q)
} else {
  pca_addiction_PC2_PC3_aspect_ratio
}

### Circle Plot
circle_plot <- as.data.frame (res$var$coord)
circle_plot$var <- rownames (circle_plot)
# reinst_annotation$Label_variable

# merging with annotation tbl
circle_plot_annotation_merged <- merge (circle_plot, reinst_annotation, by.x= "var", by.y = "Label_variable")

# data_annotation  <- circle_plot_annotation_merged
# neg_labels <- labels_v [which (data_annotation [[ "Dim.1"]] < 0)]
`$`(data_annotation , "Dim.1")
# data_annotation [[ pc_x]]

### Function
circle_plot_by_gr <- function (data_annotation, pc_x="Dim.1", pc_y="Dim.2", grouping_var="Variable", gr_plot_type = "arrows") {
  
  label_axis_x <- label_pc (pc_x)
  label_axis_y <- label_pc (pc_y)
  var_x <- var_pc (pc_x)
  var_y <- var_pc (pc_y)
  
  labels_v <- data_annotation$Variable
  neg_labels <- labels_v [which (data_annotation [[ pc_x ]] < 0)]
  neg_positions <- data_annotation [which (data_annotation [[ pc_x ]] < 0), c(pc_x, pc_y)]
  neg_positions$Variable <- neg_labels
  
  pos_labels <- labels_v [which (data_annotation [[ pc_x ]] >= 0)]
  pos_positions <- data_annotation [which (data_annotation [[ pc_x ]] >= 0), c(pc_x, pc_y)]
  pos_positions$Variable <- pos_labels
  
  angle <- seq(-pi, pi, length = 50)
  df.circle <- data.frame(x = sin(angle), y = cos(angle))
  
  pos_positions_plot <- pos_positions

  pos_positions_plot [[ pc_x]] <- pos_positions [[ pc_x]] - 0.025
  pos_positions_plot [[ pc_y]] <- pos_positions [[ pc_y]] - 0.02
  
  neg_positions_plot <- neg_positions
  neg_positions_plot [[ pc_x]] <- neg_positions [[ pc_x]] - 0.01
  neg_positions_plot [[ pc_y]] <- neg_positions [[ pc_y]] - 0.05
  list_labels_pc <- list("color" = "red", "shape" = "square", "length" = 5)
  
  if (grouping_var=="Variable") { 
    # data_annotation$Variable <- as.factor(data_annotation$Variable)
    p_circle_plot <- ggplot(data_annotation) + 
      # geom_segment (data=data_annotation, aes_string(x=0, y=0, xend=pc_x, yend=pc_y), 
                    # arrow=arrow(length=unit(0.2,"cm")), alpha=1, size=1, colour="red") +
      geom_segment (data=data_annotation, 
                    aes_string(colour=grouping_var, x=0, y=0, xend=pc_x, yend=pc_y),
                    arrow=arrow(length=unit(0.2,"cm")), 
                    alpha=1, 
                    size=1, 
                    show.legend = FALSE) +
      xlim (c(-1.2, 1.2)) + ylim (c(-1.2, 1.2)) +
      # geom_text (data=neg_positions_plot,
      #            aes_string (colour=grouping_var, x=pc_x, y=pc_y, label=grouping_var, hjust=0.1),
      #            show.legend = FALSE,
      #            size=size_text_circle) +
      # geom_text (data=pos_positions_plot,
      #            aes_string (colour=grouping_var, x=pc_x, y=pc_y, label=grouping_var, hjust=0.1),
      #            show.legend = FALSE,
      #            size=size_text_circle) +
      geom_text (data=data_annotation,
                 aes_string (colour=grouping_var, x=pc_x, y=pc_y, label=grouping_var, hjust=0.1),
                 show.legend = FALSE,
                 size=size_text_circle) +
      geom_vline (xintercept = 0, linetype="dotted") +
      geom_hline (yintercept = 0, linetype="dotted") +
      labs (title = title_var_loadings, 
            x = paste("\n", label_axis_x, "(", var_x, "% of variance)", sep=""), 
            y=paste(label_axis_y, " (", var_y, "% of variance)\n", sep = "")) +
      geom_polygon (data = df.circle, aes(x, y), alpha=1, colour="black", fill=NA, size=1) +
      scale_color_manual (values =  cb_palette_adapt) +
      coord_fixed()
    
    return(p_circle_plot) 
    
  } else {
    if (gr_plot_type == "arrows") {
      p_circle_arrows_by_gr <- ggplot(data_annotation) +
                                geom_segment (data=circle_plot_annotation_merged,
                                              aes_string (colour=grouping_var, x=0, y=0, xend=pc_x, yend=pc_y),
                                              arrow=arrow(length=unit(0.35,"cm")), alpha=1, size=2) +
                                scale_x_continuous(limits=c(-1.3, 1.3), breaks=(c(-1,0,1))) +
                                scale_y_continuous(limits=c(-1.3, 1.3), breaks=(c(-1,0,1))) +
                                scale_color_manual(values = cb_palette_adapt) +
                                geom_text (data=neg_positions_plot, aes (x=get(pc_x), y=get(pc_y), hjust=0, label=neg_labels), show.legend = FALSE, size=size_text_circle) +
                                geom_text (data=pos_positions_plot, aes (x=get(pc_x), y=get(pc_y), hjust=0, label=pos_labels), show.legend = FALSE, size=size_text_circle) +
                                geom_vline (xintercept = 0, linetype="dotted") +
                                geom_hline (yintercept=0, linetype="dotted") +
                                labs (title = title_var_loadings, x = paste("\n", label_axis_x, "(", var_x, "% of variance)", sep=""),
                                      y=paste(label_axis_y, " (", var_y, "% of variance)\n", sep = "")) +
                                geom_polygon (data = df.circle, aes(x, y), alpha=1, colour="black", fill=NA, size=1) +
                                guides(color=guide_legend(guide_legend(title = "Annotation"))) +
                                theme (legend.key = element_blank()) +
                                coord_fixed()

            return (p_circle_arrows_by_gr)
  
    } else if (gr_plot_type == "labels") {

        p_circle_arrows_labels <- ggplot(circle_plot_annotation_merged,) +
          geom_text (aes(colour=Annotation, x=get(pc_x), y=get(pc_x), label=labels_v), show.legend = FALSE, size=size_text_circle, fontface="bold", vjust=-0.4) +
          # geom_label (aes(fill=Annotation, x=get(pc_x), y=get(pc_x), label=labels_v), colour="white",show.legend = FALSE, size=size_text_circle, fontface="bold", vjust=-0.4) +
          scale_fill_manual(values = cb_palette_adapt) +
          geom_point(aes(colour=Annotation, x=get(pc_x), y=get(pc_y)), size=0) +
          scale_color_manual(values = cb_palette_adapt) +
          xlim (c(-1.2, 1.2)) + ylim (c(-1.2, 1.2)) +
          labs (title = title_var_loadings,
                x = paste("\n", label_axis_x, "(", var_x, "% of variance)", sep=""),
                y=paste(label_axis_y, " (", var_y, "% of variance)\n", sep = "")) +
          geom_vline(xintercept = 0, linetype = "longdash") +
          geom_hline(yintercept = 0, linetype = "longdash") +
          theme (legend.key = element_blank(), legend.key.height = unit (0.8, "line"), legend.title=element_blank()) +
          guides (colour = guide_legend (override.aes = list(size = 3))) +
          coord_fixed()
        
        return (p_circle_arrows_labels)
    }
  }
}

circle_plot_pc1_pc2 <- circle_plot_by_gr(circle_plot_annotation_merged)
circle_plot_pc1_pc3 <- circle_plot_by_gr(circle_plot_annotation_merged, pc_x="Dim.1", pc_y="Dim.3")
circle_plot_pc2_pc3 <- circle_plot_by_gr(circle_plot_annotation_merged, pc_x="Dim.2", pc_y="Dim.3")

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

circle_plot_by_behavior_pc1_pc2 <- circle_plot_by_gr(circle_plot_annotation_merged, grouping_var="Annotation")
circle_plot_by_behavior_pc1_pc3 <- circle_plot_by_gr(circle_plot_annotation_merged, pc_x="Dim.1", pc_y="Dim.3", grouping_var="Annotation")
circle_plot_by_behavior_pc2_pc3 <- circle_plot_by_gr(circle_plot_annotation_merged, pc_x="Dim.2", pc_y="Dim.3", grouping_var="Annotation")

if (save_plot) {
  ggsave (circle_plot_by_behavior_pc1_pc2, 
          file=paste(out_folder, name_out, "_circle_behavior_PC1_vs_PC2", extension_img, sep=""), 
          width = 15, height = 15, dpi=dpi_q)
} else {
  circle_plot_by_behavior_pc1_pc2
}

if (save_plot) {
  ggsave (circle_plot_by_behavior_pc1_pc3, 
          file=paste(out_folder, name_out, "_circle_behavior_PC1_vs_PC3", extension_img, sep=""), 
          width = 15, height = 15, dpi=dpi_q)
} else {
  circle_plot_by_behavior_pc1_pc3
}

if (save_plot) {
  ggsave (circle_plot_by_behavior_pc2_pc3, 
          file=paste(out_folder, name_out, "_circle_behavior_PC2_vs_PC3", extension_img, sep=""), 
          width = 15, height = 15, dpi=dpi_q)
} else {
  circle_plot_by_behavior_pc2_pc3
}


# circle_plot_by_behavior_pc1_pc2_labels <- circle_plot_by_gr(circle_plot_annotation_merged, grouping_var="Annotation", gr_plot_type="labels")
# circle_plot_by_behavior_pc1_pc3_labels <- circle_plot_by_gr(circle_plot_annotation_merged, pc_x="Dim.1", pc_y="Dim.3", grouping_var="Annotation", gr_plot_type="labels")
# circle_plot_by_behavior_pc2_pc3_labels <- circle_plot_by_gr(circle_plot_annotation_merged, pc_x="Dim.2", pc_y="Dim.3", grouping_var="Annotation", gr_plot_type="labels")

circle_plot_by_period_pc1_pc2 <- circle_plot_by_gr(circle_plot_annotation_merged, grouping_var="Period")
circle_plot_by_period_pc1_pc3 <- circle_plot_by_gr(circle_plot_annotation_merged, pc_x="Dim.1", pc_y="Dim.3", grouping_var="Period")
circle_plot_by_period_pc2_pc3 <- circle_plot_by_gr(circle_plot_annotation_merged, pc_x="Dim.2", pc_y="Dim.3", grouping_var="Period")

# plot_name <- "PCA_factors_addiction"
# plot_name <- "PCA_factors_addiction_1st_analysis_PCA_selection_variables"

######## AQUI
if (save_plot) {
  ggsave (circle_plot_by_period_pc1_pc2, 
          file=paste(out_folder, name_out, "_circle_period_PC1_vs_PC2", extension_img, sep=""), 
          width = 15, height = 15, dpi=dpi_q)
} else {
  circle_plot_by_period_pc1_pc2
}

if (save_plot) {
  ggsave (circle_plot_by_period_pc1_pc3, 
          file=paste(out_folder, name_out, "_circle_period_PC1_vs_PC3", extension_img, sep=""), 
          width = 15, height = 15, dpi=dpi_q)
} else {
  circle_plot_by_period_pc1_pc3
}

if (save_plot) {
  ggsave (circle_plot_by_period_pc2_pc3, 
          file=paste(out_folder, name_out, "_circle_period_PC2_vs_PC3", extension_img, sep=""), 
          width = 15, height = 15, dpi=dpi_q)
} else {
  circle_plot_by_period_pc2_pc3
}

############
## BARPLOT
pca_barPlot <- function (pca_coord="", sel_pc="Dim.1") {
  
  df.bars <- cbind (as.numeric(sort(pca_coord[ ,sel_pc ]^2/sum(pca_coord[,sel_pc ]^2)*100,decreasing=TRUE)), 
                    names(pca_coord[ ,sel_pc ])[order(pca_coord[ ,sel_pc ]^2,decreasing=TRUE)])
  df.bars_to_plot <- as.data.frame (df.bars)
  df.bars_to_plot$index <- as.factor (df.bars_to_plot$V2)
  df.bars_to_plot$value <- as.numeric(sort(res$var$coord[ ,sel_pc ]^2/sum(res$var$coord[ ,sel_pc ]^2)*100,decreasing=TRUE))
  df.bars_to_plot$index <- factor(df.bars_to_plot$index, levels = df.bars_to_plot$index[order(df.bars_to_plot$value, decreasing=TRUE)])
  
  # Filtering only the top contributors more than 2 %
  threshold <- 2
  df.bars_to_plot <- df.bars_to_plot [df.bars_to_plot$value > threshold, ]
  
  label_pc <- label_pc (sel_pc)
  var_pc <- var_pc (sel_pc)

  title_bars <- paste ("Variable contribution to ", label_pc, "\n", sep="")
  max_y <- ceiling(max(df.bars_to_plot$value) * 1.20)
  # print (df.bars_to_plot)
  bars_plot <- ggplot (data=df.bars_to_plot, aes(x=index, y=value)) + 
    scale_y_continuous (name="", breaks=seq(0, max_y, 5), limits=c(0, max_y)) +
    geom_bar (stat="identity", fill="gray", width=0.8) + 
    labs (title = title_bars, x = "", y="Contribution in %\n") +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  
  return (bars_plot)
}

bars_plot_pc1 <- pca_barPlot(res$var$coord, sel_pc="Dim.1")
bars_plot_pc2 <- pca_barPlot(res$var$coord, sel_pc="Dim.2")
bars_plot_pc3 <- pca_barPlot(res$var$coord, sel_pc="Dim.3")

bars_plot_pc1 
bars_plot_pc2
bars_plot_pc3 

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















dailyInt_theme <- theme_update (axis.title.x = element_text (size=base_size * 2, face="bold"),
                                axis.title.y = element_text (size=base_size * 2, angle = 90, face="bold"),
                                plot.title = element_text (size=base_size * 2, face="bold"))





###############
### Circle Plot



####################
####################
####################
####################




