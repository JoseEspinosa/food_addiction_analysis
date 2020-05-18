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
save_plot = FALSE

dpi_q <- 100
size_text_circle <- 5.5
title_var_loadings <- "\nVariables factor map\n" 
out_folder <- "./figures/"  
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
data_filtered <- subset (data, select=-c(Mice, Genotype, Food.pellets, Addiction_categorization_LP)) 

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

res = PCA (data_filtered, scale.unit=TRUE)

# I set a vector of colors for all the plots I just have to set the number of colours that I need for the plots
v_colours <- c("red", "gray", "blue", "lightblue", "magenta", "orange", "darkgreen")
cb_palette_adapt <- rep(c("#999999", "#009E73", "#0072B2","#E69F00", "#0072B2", "#D55E00", "#CC79A7"), 3)

# Variance of PC1 and PC2
var_PC1 <- round (res$eig [1,2], 1)
var_PC2 <- round (res$eig [2,2], 1)
var_PC3 <- round (res$eig [3,2], 1)

# PC coordinates are store here
# Convert PCA results to data frame
pca2plot <- as.data.frame (res$ind$coord)
pca2plot$id <- data$Mice
head(data)

# Changes labels of the groups, some mice are categorized as addicted and the rest as NA (non-addicted)
data$Addiction_categorization_LP[is.na(data$Addiction_categorization_LP)]  <- "N"
pca2plot$group <- data$Addiction_categorization_LP

## PC1 vs. PC3
title_p <- paste ("PCA addiction, PC1 vs. PC2\n", sep="")

pca_addiction_PC1_PC2 <- ggplot (pca2plot, aes(x=Dim.1, y=Dim.2, colour=group)) + 
                 geom_point (size = 3.5, show.legend = T) + 
                 scale_color_manual(values=c("orange", "blue")) +
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
  scale_color_manual(values=c("orange", "blue")) +
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
  scale_color_manual(values=c("orange", "blue")) +
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

###############
### Circle Plot
circle_plot <- as.data.frame (res$var$coord)
circle_plot$var <- rownames (circle_plot)

# merging with annotation tbl
circle_plot_annotation_merged <- merge (circle_plot, reinst_annotation, by.x= "var", by.y = "Label_variable")
labels_v <- circle_plot_annotation_merged$Variable
neg_labels <- labels_v [which (circle_plot_annotation_merged$Dim.1 < 0)]
neg_positions <- circle_plot_annotation_merged [which (circle_plot_annotation_merged$Dim.1 < 0), c("Dim.1", "Dim.2")]

pos_labels <- labels_v [which (circle_plot_annotation_merged$Dim.1 >= 0)]
pos_positions <- circle_plot_annotation_merged [which (circle_plot_annotation_merged$Dim.1 >= 0), c("Dim.1", "Dim.2")]

angle <- seq(-pi, pi, length = 50)
df.circle <- data.frame(x = sin(angle), y = cos(angle))

pos_positions_plot <- pos_positions
pos_positions_plot$Dim.1 <- pos_positions$Dim.1 - 0.025
pos_positions_plot$Dim.2 <- pos_positions$Dim.2 + 0.02

neg_positions_plot <- neg_positions
neg_positions_plot$Dim.1 <- neg_positions$Dim.1 #- 0.01
neg_positions_plot$Dim.2 <- neg_positions$Dim.2 + 0.05

p_circle_plot <- ggplot(circle_plot_annotation_merged) + 
  geom_segment (data=circle_plot, aes(x=0, y=0, xend=Dim.1, yend=Dim.2), 
                arrow=arrow(length=unit(0.2,"cm")), alpha=1, size=1, colour="red") +
  xlim (c(-1.2, 1.2)) + ylim (c(-1.2, 1.2)) +
  # geom_text (data=neg_positions_plot, aes (x=Dim.1, y=Dim.2, label=neg_labels, hjust=1.2), show.legend = FALSE, size=size_text_circle) +
  geom_text (data=neg_positions_plot, aes (x=Dim.1, y=Dim.2, label=neg_labels, hjust=0.1), show.legend = FALSE, size=size_text_circle) +
  # geom_text (data=pos_positions_plot, aes (x=Dim.1, y=Dim.2, label=pos_labels, hjust=-0.3), show.legend = FALSE, size=size_text_circle) +
  geom_text (data=pos_positions_plot, aes (x=Dim.1, y=Dim.2, label=pos_labels, hjust=0.1), show.legend = FALSE, size=size_text_circle) +
  geom_vline (xintercept = 0, linetype="dotted") +
  geom_hline (yintercept=0, linetype="dotted") +
  labs (title = title_var_loadings, x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
        y=paste("PC2 (", var_PC2, "% of variance)\n", sep = "")) +
  geom_polygon (data = df.circle, aes(x, y), alpha=1, colour="black", fill=NA, size=1) #+
#   theme(axis.title.x = element_text(size = size_axis)) +
#   theme(axis.title.y = element_text(size = size_axis))

p_circle_plot_coord_fixed <- p_circle_plot + coord_fixed() #+ 
#   theme(plot.title = element_text(size=22)) + 
#   theme(axis.title.x = element_text(size =22)) +
#   theme(axis.title.y = element_text(size =22))
p_circle_plot_coord_fixed

if (save_plot) {
  ggsave (p_circle_plot_coord_fixed, file=paste(home, dir_plots, "circle_annotated_behavior", img_format, sep=""),
          width = 15, height = 15, dpi=dpi_q)
} else {
  p_circle_plot_coord_fixed 
}

## Plotting by type of behavioral annotation
circle_plot_annotation_merged$Annotation
p_circle_plot_by_gr <- ggplot(circle_plot_annotation_merged) + 
  geom_segment (data=circle_plot_annotation_merged, aes(colour=Annotation, x=0, y=0, xend=Dim.1, yend=Dim.2), 
                arrow=arrow(length=unit(0.35,"cm")), alpha=1, size=2) +
  scale_x_continuous(limits=c(-1.3, 1.3), breaks=(c(-1,0,1))) +
  scale_y_continuous(limits=c(-1.3, 1.3), breaks=(c(-1,0,1))) +
  #                        xlim (c(-1.2, 1.2)) + ylim (c(-1.2, 1.2)) +
  scale_color_manual(values = cb_palette_adapt) +
  geom_text (data=neg_positions_plot, aes (x=Dim.1, y=Dim.2, label=neg_labels, hjust=1.2), show.legend = FALSE, size=size_text_circle) + 
  geom_text (data=pos_positions_plot, aes (x=Dim.1, y=Dim.2, label=pos_labels, hjust=-0.3), show.legend = FALSE, size=size_text_circle) +
  geom_vline (xintercept = 0, linetype="dotted") +
  geom_hline (yintercept=0, linetype="dotted") +
  labs (title = title_var_loadings, x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
        y=paste("PC2 (", var_PC2, "% of variance)\n", sep = "")) +
  geom_polygon (data = df.circle, aes(x, y), alpha=1, colour="black", fill=NA, size=1) +
  guides(color=guide_legend(guide_legend(title = "Annotation"))) +
  theme (legend.key = element_blank())

####################################
## Same thing but without arrows
# aes(colour=annot_gr,
p_circle_points <- ggplot(circle_plot_annotation_merged,) + 
  geom_text (aes(colour=Annotation, x=Dim.1, y=Dim.2,label=labels_v), show.legend = FALSE, size=size_text_circle, fontface="bold", vjust=-0.4) +
  # geom_label (aes(fill=Annotation, x=Dim.1, y=Dim.2,label=labels_v), colour="white",show.legend = FALSE, size=size_text_circle, fontface="bold", vjust=-0.4) +
  scale_fill_manual(values = cb_palette_adapt) +
  geom_point(aes(colour=Annotation, x=Dim.1, y=Dim.2), size=0) +
  scale_color_manual(values = cb_palette_adapt) +
  xlim (c(-1.2, 1.2)) + ylim (c(-1.2, 1.2)) + 
  labs (title = title_var_loadings) +
  labs (x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
        y=paste("PC2 (", var_PC2, "% of variance)\n", sep = "")) +
  geom_vline(xintercept = 0, linetype = "longdash") +
  geom_hline(yintercept = 0, linetype = "longdash") +
  theme (legend.key = element_blank(), legend.key.height = unit (0.8, "line"), legend.title=element_blank()) +
  guides (colour = guide_legend (override.aes = list(size = 3)))

# p_circle_points_leg <- p_circle_points + theme(legend.text = element_text(size = 20))

p_circle_points_coord_fixed <-p_circle_points + coord_fixed()
# p_circle_points_coord_fixed

# plot_name <- "PCA_factors_addiction"
plot_name <- "PCA_factors_addiction_1st_analysis_PCA_selection_variables"
if (save_plot) {
  ggsave (p_circle_points_coord_fixed, file=paste('/home/kadomu/projects/20200421_pca_behavior_elena/results/',
                                                  plot_name, extension_img, sep=""), width = 15, height = 15, dpi=dpi_q)
} else {
  p_circle_points_coord_fixed
}

############
## BARPLOT
df.bars <- cbind (as.numeric(sort(res$var$coord[,1]^2/sum(res$var$coord[,1]^2)*100,decreasing=TRUE)), names(res$var$coord[,1])[order(res$var$coord[,1]^2,decreasing=TRUE)])
df.bars_to_plot <- as.data.frame(df.bars)
df.bars_to_plot$index <- as.factor (df.bars_to_plot$V2)
# class (df.bars_to_plot$V1)
df.bars_to_plot$value <- as.numeric(sort(res$var$coord[,1]^2/sum(res$var$coord[,1]^2)*100,decreasing=TRUE))
df.bars_to_plot$index <- factor(df.bars_to_plot$index, levels = df.bars_to_plot$index[order(df.bars_to_plot$value, decreasing=TRUE)])

# PC1
# Filtering only the top contributors more than 2 %
threshold <- 2
df.bars_to_plot <- df.bars_to_plot [df.bars_to_plot$value > threshold, ]

title_b_pc1 <- paste ("Variable contribution to PC1", "\n", sep="")
max_y_pc1 <- ceiling(max(df.bars_to_plot$value) * 1.20)

bars_plot <- ggplot (data=df.bars_to_plot, aes(x=index, y=value)) + 
  scale_y_continuous (name="", breaks=seq(0, max_y_pc1, 5), limits=c(0, max_y_pc1)) +
  geom_bar (stat="identity", fill="gray", width=0.8) + 
  labs (title = title_b_pc1, x = "", y="Contribution in %\n") +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

if (save_plot) {
  ggsave (bars_plot, file=paste(home, "/old_data/figures/", 
                                "bars_PC1_",  "impulsivity.tiff", sep=""), width = 15, height = 12, dpi=dpi_q)
} else {
  bars_plot
}
  
# PC2
title_b <- paste ("Variable contribution to PC2", "\n", sep="")

df.bars_PC2 <- cbind (as.numeric(sort(res$var$coord[,2]^2/sum(res$var$coord[,2]^2)*100,decreasing=TRUE)), names(res$var$coord[,2])[order(res$var$coord[,2]^2,decreasing=TRUE)])
df.bars_to_plot_PC2 <- as.data.frame(df.bars_PC2)
df.bars_to_plot_PC2$index <- as.factor (df.bars_to_plot_PC2$V2)

# df.bars_to_plot_PC2$value <- as.numeric(sort(res$var$coord[,2]^2/sum(res$var$coord[,2]^2)*100,decreasing=TRUE))
df.bars_to_plot_PC2$value <- as.numeric(sort(res$var$coord[,2]^2/sum(res$var$coord[,2]^2)*100,decreasing=TRUE))

# Filtering only the top contributors more than 2 %
threshold_pc2 <- 2
df.bars_to_plot_PC2 <- df.bars_to_plot_PC2 [df.bars_to_plot_PC2$value > threshold_pc2, ]
df.bars_to_plot_PC2$index
df.bars_to_plot_PC2$index <- factor(df.bars_to_plot_PC2$index, levels = df.bars_to_plot_PC2$index[order(df.bars_to_plot_PC2$value, decreasing=TRUE)])

title_b_pc2 <- paste ("Variable contribution to PC2", "\n", sep="")
max_y_pc2 <- ceiling(max(df.bars_to_plot_PC2$value) * 1.20)

bars_plot_PC2 <- ggplot (data=df.bars_to_plot_PC2, aes(x=index, y=value)) +
  scale_y_continuous (name="", breaks=seq(0, max_y_pc2, 5), limits=c(0, max_y_pc2)) +
  geom_bar (stat="identity", fill="gray", width=0.8) + 
  labs (title = title_b_pc2, x = "", y="Contribution in %\n") +
  theme (axis.text.x=element_text(angle=45, vjust=1, hjust=1))

if (save_plot) {
  ggsave (bars_plot_PC2, file=paste(home, "/old_data/figures/", 
                                    "bars_PC2_",  phase, "Phase.tiff", sep=""), width = 15, height = 12, dpi=dpi_q)
} else {
  bars_plot_PC2
}

# PC3
title_b <- paste ("Variable contribution to PC3", "\n", sep="")
df.bars_PC3 <- cbind (as.numeric(sort(res$var$coord[,3]^2/sum(res$var$coord[,3]^2)*100,decreasing=TRUE)), names(res$var$coord[,3])[order(res$var$coord[,3]^2,decreasing=TRUE)])
df.bars_to_plot_PC3 <- as.data.frame(df.bars_PC3)
df.bars_to_plot_PC3$index <- as.factor (df.bars_to_plot_PC3$V2)
df.bars_to_plot_PC3$value <- as.numeric(sort(res$var$coord[,3]^2/sum(res$var$coord[,3]^2)*100,decreasing=TRUE))

# Filtering only the top contributors more than 2 %
threshold_pc3 <- 2
df.bars_to_plot_PC3 <- df.bars_to_plot_PC3 [df.bars_to_plot_PC3$value > threshold_pc3, ]
df.bars_to_plot_PC3$index
df.bars_to_plot_PC3$index <- factor(df.bars_to_plot_PC3$index, levels = df.bars_to_plot_PC3$index[order(df.bars_to_plot_PC3$value, decreasing=TRUE)])

title_b_pc3 <- paste ("Variable contribution to PC3", "\n", sep="")
max_y_pc3 <- ceiling(max(df.bars_to_plot_PC3$value) * 1.20)

# Variability explained by PC3
bars_plot_PC3 <- ggplot (data=df.bars_to_plot_PC3, aes(x=index, y=value)) +
  scale_y_continuous (name="", breaks=seq(0, max_y_pc3, 5), limits=c(0, max_y_pc3)) +
  geom_bar (stat="identity", fill="gray", width=0.8) + 
  labs (title = title_b_pc3, x = "", y="Contribution in %\n") +
  theme (axis.text.x=element_text(angle=45, vjust=1, hjust=1))

if (save_plot) {
  ggsave (bars_plot_PC3, file=paste(home, "/old_data/figures/", 
                                    "bars_PC3_",  phase, "Phase.tiff", sep=""), width = 15, height = 12, dpi=dpi_q)
} else {
  bars_plot_PC3
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



# Plotting the variables by experimental phase
circle_plot$var <- rownames (circle_plot)

circle_plot$var <- gsub ("day", "", circle_plot$var)
circle_plot$var <- gsub ("dep_IMUL", "dep", circle_plot$var)
circle_plot$var <- gsub ("adlib_IMPUL", "adlib" , circle_plot$var)

circle_plot$varGroup <- circle_plot$var
circle_plot$varGroup [grep("^dep", circle_plot$var)] <- "dep_impul"
circle_plot$varGroup [grep("^adlib", circle_plot$var)] <- "ad_impul"

# circle_plot$varGroup [grep("^adlib_act", circle_plot$var)] <- "adlib_act"
# circle_plot$varGroup [grep("^adlib_inact", circle_plot$var)] <- "adlib_in"
# circle_plot$varGroup [grep("^ex_act", circle_plot$var)] <- "ex_act"
# circle_plot$varGroup [grep("^ex_inact", circle_plot$var)] <- "ex_inact"
# circle_plot$varGroup [c(81:length(circle_plot$varGroup))] <- "others"
as.factor(circle_plot$varGroup)
colnames (circle_plot) <- c("Dim.1", "Dim.2", "Dim.3", "Dim.4", "Dim.5", "var", "varGroup")

# I only need the nummber of session for each of them
circle_plot$session <- gsub("^dep_", "", circle_plot$var) 
circle_plot$session <- gsub("^adlib_", "", circle_plot$session)

# circle_plot$session <- gsub("^adlib_inact_", "", circle_plot$session)
# circle_plot$session <- gsub("^ex_act_", "", circle_plot$session)
# circle_plot$session <- gsub("^ex_inact_", "", circle_plot$session)

############
# Doing a circle plot with arrows coloured by experimental phase

# Adding session to the circle_plot df to plot them
neg_labels <- labels_v [which (circle_plot$Dim.1 < 0)]
neg_positions <- circle_plot [which (circle_plot$Dim.1 < 0), c(1,2,8)]

pos_labels <- labels_v [which (circle_plot$Dim.1 >= 0)]
pos_positions <- circle_plot [which (circle_plot$Dim.1 >= 0), c(1,2,8)]

title_c <- paste ("PCA of the variables ","- presses during time-out\n", sep="")
p_circle_plot_colors <- ggplot(circle_plot) + 
  geom_segment (data=circle_plot, aes(colour=varGroup, x=0, y=0, xend=Dim.1, yend=Dim.2), 
                arrow=arrow(length=unit(0.2,"cm")), alpha=1, size=1) +
  scale_color_manual (values = v_colours) +
  xlim (c(-1.2, 1.2)) + ylim (c(-1.2, 1.2)) +
  geom_text (data=neg_positions, aes (x=Dim.1, y=Dim.2, label=session, hjust=0.9, vjust=-0.4), 
             show.legend = FALSE, size=5) + 
  geom_text (data=pos_positions, aes (x=Dim.1, y=Dim.2, label=session, hjust=-0.2), 
             show.legend = FALSE, size=5) +
  geom_vline (xintercept = 0, linetype="dotted") +
  geom_hline (yintercept=0, linetype="dotted") +
  labs (title = title_c, x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
        y=paste("PC2 (", var_PC2, "% of variance)\n", sep = "")) +
  geom_polygon (data = df.circle, aes(x, y), alpha=1, colour="black", fill=NA, size=1) +
  theme (legend.key = element_blank(), legend.key.height = unit (1.5, "line"), 
         legend.title=element_blank()) 
base_size <- 10

p_circle_plot_colors
dailyInt_theme <- theme_update (axis.title.x = element_text (size=base_size * 2, face="bold"),
                                axis.title.y = element_text (size=base_size * 2, angle = 90, face="bold"),
                                plot.title = element_text (size=base_size * 2, face="bold"))
p_circle_plot_colors

if (save_plot) {
  ggsave (p_circle_plot_colors, 
          file=paste(home, "/old_data/figures/", 
                     "circle_color_act_",  "impulsivity.tiff", sep=""), width = 15, height = 12, dpi=dpi_q)
}

# # Colour circle plot by session 1-5 5-10 10-15 15-20
# circle_plot 
# circle_plot$session_bin <- "" 
# 
# circle_plot [which (as.numeric (circle_plot$session) < 6), "session_bin"] <- "1_5"
# circle_plot [which (as.numeric (circle_plot$session) > 5 & as.numeric (circle_plot$session)< 11), "session_bin"] <- "6_10"
# circle_plot [which (as.numeric (circle_plot$session) > 10 & as.numeric (circle_plot$session)< 16), "session_bin"] <- "11_15"
# circle_plot [which (as.numeric (circle_plot$session) > 15), "session_bin"] <- "16_20"
# 
# # Plot with arrows coloured by session bin
# p_circle_plot_colors_bin <- ggplot(circle_plot) + 
#   geom_segment (data=circle_plot, aes(colour=session_bin, x=0, y=0, xend=Dim.1, yend=Dim.2), 
#                 arrow=arrow(length=unit(0.2,"cm")), alpha=1, size=1) +
#   scale_color_manual (values = v_colours) +
#   xlim (c(-1.2, 1.2)) + ylim (c(-1.2, 1.2)) +
#   geom_text (data=neg_positions, aes (x=Dim.1, y=Dim.2, label=session, hjust=0.9, vjust=-0.4), 
#              show.legend = FALSE, size=5) + 
#   geom_text (data=pos_positions, aes (x=Dim.1, y=Dim.2, label=session, hjust=-0.2), 
#              show.legend = FALSE, size=5) +
#   geom_vline (xintercept = 0, linetype="dotted") +
#   geom_hline (yintercept=0, linetype="dotted") +
#   labs (title = title_c, x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
#         y=paste("PC2 (", var_PC2, "% of variance)\n", sep = "")) +
#   geom_polygon (data = df.circle, aes(x, y), alpha=1, colour="black", fill=NA, size=1) +
#   theme (legend.key = element_blank(), legend.key.height = unit (1.5, "line"), 
#          legend.title=element_blank()) 
# 
# base_size <- 10
# 
# dailyInt_theme <- theme_update (axis.title.x = element_text (size=base_size * 2, face="bold"),
#                                 axis.title.y = element_text (size=base_size * 2, angle = 90, face="bold"),
#                                 plot.title = element_text (size=base_size * 2, face="bold"))
# p_circle_plot_colors_bin

# ggsave (p_circle_plot_colors_bin, file=paste(home, "/old_data/figures/", 
#                                              "circle_color_bin_",  phase, "Phase.tiff", sep=""), width = 15, height = 12, dpi=dpi_q)

# Doing the same plot as above by colours but in this case facet
p_var_by_group_scale_free <- ggplot(circle_plot) + 
  labs (title = title_c, x = paste("PC1 (", var_PC1, "% of variance)", sep=""), 
        y=paste("PC2 (", var_PC2, "% of variance)", sep = "")) +
  geom_text (aes(colour=varGroup, x=Dim.1, y=Dim.2, label=session), show.legend = FALSE, size=5) +
  scale_color_manual (values = c("red", "gray", "blue", "lightblue", "magenta", "orange", "darkgreen")) +
  facet_wrap(~varGroup, scales="free")

p_var_by_group_scale_free
# ggsave (p_var_by_group_scale_free, file=paste(home, "/old_data/figures/", 
#                                    "points_act_facet_",  phase, "Phase.tiff", sep=""), width = 15, height = 12, dpi=dpi_q)

p_var_by_group <- ggplot(circle_plot) + 
  xlim (c(-1, 1)) + ylim (c(-1, 1)) +
  geom_text (aes(colour=varGroup, x=Dim.1, y=Dim.2, label=session), show.legend = FALSE, size=5, vjust=-0.4) +
  geom_point(aes(colour=varGroup, x=Dim.1, y=Dim.2), size=3)+
  scale_color_manual (values = v_colours) +
  labs (title = title_c)
labs (x = paste("\nPC1 (", var_PC1, "% of variance)", sep=""), 
      y=paste("PC2 (", var_PC2, "% of variance)\n", sep = "")) +
  theme (legend.key = element_blank(), legend.key.height = unit (1.5, "line"), 
         legend.title=element_blank()) 

p_var_by_group

if (save_plot) {
  # ggsave (p_var_by_group, file=paste(home, "/old_data/figures/", 
  #                                          "points_act_",  phase, "Phase.tiff", sep=""), width = 15, height = 12, dpi=dpi_q)
}

