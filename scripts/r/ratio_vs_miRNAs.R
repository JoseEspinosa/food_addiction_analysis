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

m=read.table(ra,sep="\t",skip=3,row.names=1)
M=t(m)
rownames(M)=ma[1,]
gr=as.character(ma[2,])

library(propr)

#check differential proportionality:

test=apply(M,2,min)
M2=M[,-which(test==0)]

res2=propd(M2,group=gr)
res2=updateCutoffs(res2,seq(0.5,0.9,0.01))
res2

#Bacteroidales S24-7 group_unidentified"
#"Prevotellaceae_UCG001"   
M2_df <- data.frame(M2)
M2_df$mouse_id <- row.names(M2_df)
M2_df$ratio_murib_prevo <- M2_df[,5] / M2_df[,19]
M2_df_sel <- subset(M2_df, select=c("mouse_id","ratio_murib_prevo"))

## Copy code from the other script
# miRNA_clr <- miRNAs_data_selected

## Merge tables
ratio_miRNA_merged <- merge (miRNA_clr, M2_df_sel, by= "mouse_id")

# Selected miRNAs
# first_miRNA <- 'mmu-miR-876-5p'; last_miRNA <- 'mmu-miR-192-5p';
# axis_text_size_x <- 16; size_p_values <- 5; angle_reg <- 0; microbio_set <- "selected"
# width_p <- 20; height_p <- 12

# All miRNAs
first_miRNA <- 'mmu-miR-764-3p'; last_miRNA <- 'mmu-miR-26b-5p';
axis_text_size_x <- 8; size_p_values <- 4; angle_reg<-270; microbio_set <- "zero_filt"
width_p <- 45; height_p <- 14

## Gather
data <- gather(ratio_miRNA_merged, miRNA, value, first_miRNA:last_miRNA)%>%
        gather(ratio, ratio_value, ratio_murib_prevo)
data_nest <- group_by(data, miRNA, ratio) %>% nest()

# cor_method <- "pearson"
cor_method <- "spearman"
cor_fun <- function(df) cor.test(df$value, df$ratio_value, method = cor_method) %>% tidy()

library(broom) # Convert results of statistical functions (lm, t.test, cor.test, etc.) into tidy tables

data_nest <- mutate(data_nest, model = map(data, cor_fun))

corr_pr <- select(data_nest, -data) %>% unnest()
corr_pr <- mutate(corr_pr, sig = ifelse(p.value <0.05, "Sig.", "Non Sig."))

###########
## FDR
fdr_cutoff <- 0.2
corr_pr$fdr <- p.adjust(corr_pr$p.value, method = 'BH', n = length(corr_pr$p.value))
corr_pr <- mutate(corr_pr, sig = ifelse(fdr < fdr_cutoff, "Sig.", "Non Sig."))
corr_pr$fdr
subset(corr_pr, fdr < 0.2)
## Plot
# title_p <- paste("Correlations between ratio Muribaculaceae/Prevotellacea \nand miRNAS",
#                  "filtered zeros, transformed\n")

title_p <- paste("Correlations between ratio Muribaculaceae/Prevotellacea \nand all miRNAS",
                 "filtered zeros, transformed\n")

hm <- ggplot() + geom_tile(data = corr_pr,
                           aes(ratio, miRNA, fill = estimate),
                           size = 1,
                           colour = "white") +
  ## black lines around tiles
  geom_tile(data = filter(corr_pr, sig == "Sig."),
            aes(ratio, miRNA),
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
# axis_text_size_x <- 8; size_p_values <- 4; angle_reg<-270; microbio_set <- "all"
# width_p <- 45; height_p <- 30

# width_p <- 20; height_p <- 14
ggsave (hm, file=paste0(out_dir, "heatmap_", "all_filtered_miRNA_transf_", "ratio_murib_prevo", extension_img),
        width = width_p, height = height_p, dpi=dpi_q)

