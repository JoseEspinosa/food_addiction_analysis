# Loading functions:
source ("./scripts/r/graph_parameters.R")

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

# Return a list with the variance for each of the PCs
get_var_pcs <- function (pca_r) {
  
  list_v <- list()
  dim_name <- "Dim."
  # Variance of PCs
  for (i in 1 : length(pca_r$eig [,2])) {
    id = paste0(dim_name, i)
    list_v [[id]] <- round (pca_r$eig [i,2], 1)
  }
  
  return (list_v)
}

pca_plot <- function (pca_r, pc_x="Dim.1", pc_y="Dim.2", id_v=NA, group_v=NA, title_pca="", 
                      min_x=NA, max_x=NA, min_y=NA, max_y=NA, list_variance=NA) {
  
  cb_palette <- c("#E69F00","#D55E00", "#56B4E9", "#009E73", "#999999", "#CC79A7")
  
  dim_name <- "Dim."
  
  # PC coordinates are store here
  # Convert PCA results to data frame
  pca2plot <- as.data.frame (pca_r$ind$coord)
  
  list_var <- list()
  
  if (is.na(list_variance[[1]])) {
    list_var <- get_var_pcs(pca_r)
  } else {
    list_var <- list_variance
  }
  
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
    pca2plot$id <- id_v
  } else {
    pca2plot$id <- 1:n_samples
  }
  
  if (!is.na(group_v[1])) {
    pca2plot$Group <- group_v
  } else {
    pca2plot$Group <- "NA"
  }
  
  title_p <- paste(title_pca, label_pc(pc_x), "vs.", label_pc(pc_y), "\n", sep=" ")
  title_x <- paste("\n", label_pc(pc_x), "(", list_var[[pc_x]], "% of variance)", sep="")
  title_y <- paste(label_pc(pc_y), " (", list_var[[pc_y]], "% of variance)\n", sep = "")
  
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

library(FactoMineR)

##############
## Behavior
bf <- paste0(home_dir, "/git/food_addiction_analysis/forJose/behavior.txt")
mb=read.table(bf,header=TRUE,sep="\t",row.names=1)
mb
res = PCA (mb, scale.unit=TRUE, ncp=7)

home_dir <- Sys.getenv("HOME")
taxon <- "genus"
ra <- paste0(home_dir, "/git/food_addiction_analysis/forJose/relative_abundances_by_", taxon, ".txt")
ma <- read.table(ra, sep="\t", skip=1,
                 row.names=1, colClasses = "character")
m_gr <- ma[c(1,2),]
colnames(m_gr)=ma[1,]

group <- as.data.frame(t(m_gr[c(1,2), rownames(res$ind$coord)]))
group_v <- group$Grouping
sample_id = as.vector(row.names(as.data.frame (res$ind$coord)))

title_2_plot = "PCA"
list_var <- list()
list_var <- get_var_pcs(res)

pca_addiction_PC1_PC2 <- pca_plot (res, pc_x="Dim.1", pc_y="Dim.2", group_v=group_v, 
                                   title_pca=title_2_plot, id_v = sample_id,
                                   list_variance = list_var)

pca_dim <- as.data.frame (res$ind$coord)
pca_dim$mouse_id <- row.names(pca_dim)



               