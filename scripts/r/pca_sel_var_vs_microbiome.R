#############################################################################
### Jose Espinosa-Carrasco CB-CRG Group. Nov 2020                         ###
#############################################################################
### Cross-talking between behavior and microbiome                         ###
#############################################################################

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

first_up <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
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
home_dir <- Sys.getenv("HOME")

##############
## Behavior
bf <- paste0(home_dir, "/git/food_addiction_analysis/forJose/behavior.txt")
mb=read.table(bf,header=TRUE,sep="\t",row.names=1)
mb
res = PCA (mb, scale.unit=TRUE, ncp=7)

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

get_var_pcs(res)
res$var
pca_addiction_PC1_PC2 <- pca_plot (res, pc_x="Dim.1", pc_y="Dim.2", group_v=group_v, 
                                   title_pca=title_2_plot, id_v = sample_id,
                                   list_variance = list_var)

pca_dim <- as.data.frame (res$ind$coord)
pca_dim$mouse_id <- row.names(pca_dim)

### Microbiome data
taxon <- "genus"
ra <- paste0(home_dir, "/git/food_addiction_analysis/forJose/relative_abundances_by_", taxon, ".txt")
ma <- read.table(ra, sep="\t", skip=1,
                 row.names=1, colClasses = "character")

m_microbiome=read.table(ra,sep="\t",skip=3,row.names=1)
colnames(m_microbiome)=ma[1,]

###################################
## Filter zero values
M=t(m_microbiome)
rownames(M)=ma[1,]

# ## Any genera with zeros is removed i.e. any row with a zero
test=apply(M,2,min)
## test
M_f = M[,-which(test==0)]
m_microbiome_f <- t(M_f) 

first_taxon<-'Alistipes'; last_taxon<-'Tyzzerella';
m_microbiome_to_clr <- m_microbiome_f

## all genera
# first_taxon<-'Acetatifactor'; last_taxon<-'Tyzzerella';
# m_microbiome_to_clr <- m_microbiome

###################################
## Transform microbiome data
## microbiome samples (mouse) should be in the columns for microbiome::transform
microbiome_clr <- microbiome::transform(m_microbiome_to_clr, "clr")

#######
# merging tables
microbiome_clr_df <- as.data.frame(t(microbiome_clr))
microbiome_df <- microbiome_clr_df

microbiome_df$mouse_id <- rownames(microbiome_df)

pca_microbiome_merged <- merge (microbiome_df, pca_dim, by= "mouse_id")

## Merge tables
data <- gather(pca_microbiome_merged, taxon, microbio_rel_ab, first_taxon:last_taxon)%>%
        gather(pca_dim, value, Dim.1:Dim.7)

data_nest <- group_by(data, taxon, pca_dim) %>% nest()

# cor_method <- "pearson"
cor_method <- "spearman"
cor_fun <- function(df) cor.test(df$microbio_rel_ab, df$value, method = cor_method) %>% tidy()

library(broom) # Convert results of statistical functions (lm, t.test, cor.test, etc.) into tidy tables

data_nest <- mutate(data_nest, model = map(data, cor_fun))

corr_pr <- select(data_nest, -data) %>% unnest()
corr_pr <- mutate(corr_pr, sig = ifelse(p.value <0.05, "Sig.", "Non Sig."))

## Plot
title_p <- paste("Correlations between PCA vars and", 
                 first_up(taxon), "filtered zeros, transformed ","relative abundances\n")

axis_text_size_x<-18
axis_text_size_y <- 14
hm <- ggplot() + geom_tile(data = corr_pr,
                           # aes(taxon, miRNA, fill = estimate),
                           aes(pca_dim, taxon, fill = estimate),
                           size = 1,
                           colour = "white") +
  ## black lines around tiles
  geom_tile(data = filter(corr_pr, sig == "Sig."),
            aes(pca_dim, taxon),
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

out_dir <- paste0(home_dir, "/git/food_addiction_analysis/figures/pca_sel_behavior/")
dpi_q <- 200
extension_img <- ".png"
width_p <- 20; height_p <- 14
ggsave (hm, file=paste0(out_dir, "heatmap_", "filtered_microbio_transf_", taxon, "_pcaVars", extension_img), 
        width = width_p, height = height_p, dpi=dpi_q)

######################
## microRNA data:
miRNAs_1 <- paste0(home_dir, "/git/food_addiction_analysis/forJose/miRNAs_1.txt")
miRNAs_2 <- paste0(home_dir, "/git/food_addiction_analysis/forJose/miRNAs_2.txt")

l1 = read.table(miRNAs_1, header=FALSE, row.names=1, colClasses = "character")
l2 = read.table(miRNAs_2, header=FALSE, row.names=1, colClasses = "character")

mi1=read.table(miRNAs_1, skip=1, row.names=1)
mi2=read.table(miRNAs_2, header=TRUE, row.names=1)

colnames(mi1)=l1[1,]
colnames(mi2)=l2[1,]
smi1 = mi1[rownames(mi2),]
mi = cbind(smi1,mi2)

#we miss some control mice that are in the microbiome data:
M3 = M_f[intersect(colnames(mi),rownames(M_f)),]
MI = t(mi[,intersect(colnames(mi), rownames(M_f))])

test = apply(MI,2,min)
MI2 = MI[,which(test>0)]

miRNA_to_transf <- t(MI2)
miRNA_clr <- as.data.frame(t(microbiome::transform(miRNA_to_transf, "clr")))
miRNA_clr$mouse_id <- row.names(miRNA_clr)

## Select miRNAs from Elena's table
# miRNAs_data_selected <- subset(miRNA_clr, select=c('mouse_id',
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
# 
# miRNA_clr <- miRNAs_data_selected

## Merge tables
pca_miRNA_merged <- merge (miRNA_clr, pca_dim, by= "mouse_id")
colnames(pca_miRNA_merged)[2]

# All miRNAs
# first_miRNA <- 'mmu-miR-764-3p'; last_miRNA <- 'mmu-miR-26b-5p';
# axis_text_size_x <- 8; size_p_values <- 4; angle_reg<-270; microbio_set <- "zero_filt"
# width_p <- 45; height_p <- 14

# Selected miRNAs
first_miRNA <- 'mmu-miR-876-5p'; last_miRNA <- 'mmu-miR-192-5p';
axis_text_size_x <- 16; size_p_values <- 5; angle_reg <- 0; microbio_set <- "selected"
width_p <- 20; height_p <- 12

## Gather
data <- gather(pca_miRNA_merged, miRNA, value, first_miRNA:last_miRNA)%>%
        gather(pca_dim, pca_coord, Dim.1:Dim.7)
        
data_nest <- group_by(data, miRNA, pca_dim) %>% nest()

# cor_method <- "pearson"
cor_method <- "spearman"
cor_fun <- function(df) cor.test(df$value, df$pca_coord, method = cor_method) %>% tidy()

library(broom) # Convert results of statistical functions (lm, t.test, cor.test, etc.) into tidy tables

data_nest <- mutate(data_nest, model = map(data, cor_fun))

corr_pr <- select(data_nest, -data) %>% unnest()
corr_pr <- mutate(corr_pr, sig = ifelse(p.value <0.05, "Sig.", "Non Sig."))

## Plot
title_p <- paste("Correlations between PCA vars and miRNAS",
                 "filtered zeros, transformed\n")
# title_p <- paste("Correlations between PCA vars and selected miRNAS", 
#                  " transformed\n")

# axis_text_size_x<-18
# axis_text_size_y <- 14
hm <- ggplot() + geom_tile(data = corr_pr,
                           aes(pca_dim, miRNA, fill = estimate),
                           size = 1,
                           colour = "white") +
  ## black lines around tiles
  geom_tile(data = filter(corr_pr, sig == "Sig."),
            aes(pca_dim, miRNA),
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

out_dir <- paste0(home_dir, "/git/food_addiction_analysis/figures/pca_sel_behavior/")
dpi_q <- 200
extension_img <- ".png"


# All miRNAs
axis_text_size_x <- 8; size_p_values <- 4; angle_reg<-270; microbio_set <- "all"
width_p <- 45; height_p <- 30

# width_p <- 20; height_p <- 14
ggsave (hm, file=paste0(out_dir, "heatmap_", "filtered_miRNA_transf_", "_pcaVars", extension_img),
        width = width_p, height = height_p, dpi=dpi_q)

## selected miRNAs
# ggsave (hm, file=paste0(out_dir, "heatmap_", "filtered_miRNA_selected_", "pcaVars", extension_img), 
        # width = width_p, height = height_p, dpi=dpi_q)

