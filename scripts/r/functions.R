#!/usr/bin/env Rscript

#############################################################################
### Jose Espinosa-Carrasco CB-CRG Group. June 2020                        ###
#############################################################################
### Functions to plot PCA                                                 ###
#############################################################################

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

### Function
circle_plot_by_gr <- function (data_annotation, pc_x="Dim.1", pc_y="Dim.2", grouping_var="var",
                               list_variance=NA, gr_plot_type = "arrows") {

  label_axis_x <- label_pc (pc_x)
  label_axis_y <- label_pc (pc_y)

  # if (is.na(list_variance[[1]])) {
  #   list_var <- get_var_pcs(pca_r)
  # } else {
  #   list_var <- list_variance
  # }

  var_x <- list_variance [pc_x]
  var_y <- list_variance [pc_y]

  labels_v <- data_annotation$var
  neg_labels <- labels_v [which (data_annotation [[ pc_x ]] < 0)]
  neg_positions <- data_annotation [which (data_annotation [[ pc_x ]] < 0), c(pc_x, pc_y)]
  neg_positions$var <- neg_labels

  pos_labels <- labels_v [which (data_annotation [[ pc_x ]] >= 0)]
  pos_positions <- data_annotation [which (data_annotation [[ pc_x ]] >= 0), c(pc_x, pc_y)]
  pos_positions$var <- pos_labels

  angle <- seq(-pi, pi, length = 50)
  df.circle <- data.frame(x = sin(angle), y = cos(angle))

  pos_positions_plot <- pos_positions

  pos_positions_plot [[ pc_x]] <- pos_positions [[ pc_x]] - 0.025
  pos_positions_plot [[ pc_y]] <- pos_positions [[ pc_y]] - 0.02

  neg_positions_plot <- neg_positions
  neg_positions_plot [[ pc_x]] <- neg_positions [[ pc_x]] - 0.01
  neg_positions_plot [[ pc_y]] <- neg_positions [[ pc_y]] - 0.05
  list_labels_pc <- list("color" = "red", "shape" = "square", "length" = 5)

  if (grouping_var=="var") {
    # data_annotation$var <- as.factor(data_annotation$var)
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


                                # geom_text (data=neg_positions_plot, aes (x=get(pc_x), y=get(pc_y), hjust=0, label=neg_labels), show.legend = FALSE, size=size_text_circle) +
                                geom_text (data=data_annotation, aes_string (x=pc_x, y=pc_y, hjust=0, label="var"), show.legend = FALSE, size=size_text_circle) +
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
  
  title_bars <- paste ("Variable contribution to ", label_pc, "\n", sep="")
  max_y <- ceiling(max(df.bars_to_plot$value) * 1.20)
  # print (df.bars_to_plot)
  bars_plot <- ggplot (data=df.bars_to_plot, aes(x=index, y=value)) + 
    scale_y_continuous (breaks=seq(0, max_y, 5), limits=c(0, max_y)) +
    geom_bar (stat="identity", fill="gray", width=0.8) + 
    labs (title = title_bars, x = "", y="Contribution in %\n") +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  
  return (bars_plot)
}