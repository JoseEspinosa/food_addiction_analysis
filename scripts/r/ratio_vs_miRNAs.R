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
res2
which(res2@results$theta<0.8)
res2@results[158,1:3]
Partner Pair     theta
# 158      19    5 0.5093579
colnames(res2@counts)[c(19,5)]

res2@results[65,1:3]
colnames(res2@counts)[c(12,10)]
# "Lachnoclostridium" "Enterorhabdus"    

res2@results[88,1:3]
colnames(res2@counts)[c(14,10)]
# "Lachnospiraceae_NK4A136_group" "Enterorhabdus"   

res2@results[281,1:3]
colnames(res2@counts)[c(25,5)]
# "Ruminococcaceae_UCG014"                 "Bacteroidales S24-7 group_unidentified"

res2@results[302,1:3]
colnames(res2@counts)[c(26,2)]
# "Tyzzerella"     "Alloprevotella"

res2@results[310,1:3]
colnames(res2@counts)[c(26,10)]
# "Tyzzerella"    "Enterorhabdus"

#Bacteroidales S24-7 group_unidentified"
#"Prevotellaceae_UCG001"   
M2_df <- data.frame(M2)
M2_df$mouse_id <- row.names(M2_df)
M2_df$ratio_murib_prevo <- M2_df[,5] / M2_df[,19]
M2_df_sel <- subset(M2_df, select=c("mouse_id","ratio_murib_prevo"))

### miRNAs data
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

miRNA_ori <- as.data.frame(MI2)
miRNA_ori$mouse_id <- row.names(miRNA_ori)

## Transformed data
# miRNA_to_transf <- t(MI2)
# miRNA_clr <- as.data.frame(t(microbiome::transform(miRNA_to_transf, "clr")))
# miRNA_clr$mouse_id <- row.names(miRNA_clr)

miRNA_to_sel<-miRNA_ori
# miRNA_to_sel<-miRNA_clr

miRNAs_data_selected <- subset(miRNA_to_sel, select=c('mouse_id',
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


## Merge tables
ratio_miRNA_merged <- merge (miRNAs_data_selected, M2_df_sel, by= "mouse_id")

# Selected miRNAs
first_miRNA <- 'mmu-miR-876-5p'; last_miRNA <- 'mmu-miR-192-5p';
axis_text_size_x <- 16; size_p_values <- 5; angle_reg <- 0; microbio_set <- "selected"
width_p <- 20; height_p <- 12

# All miRNAs
# first_miRNA <- 'mmu-miR-764-3p'; last_miRNA <- 'mmu-miR-26b-5p';
# axis_text_size_x <- 8; size_p_values <- 4; angle_reg<-270; microbio_set <- "zero_filt"
# width_p <- 45; height_p <- 14

## Gather
data <- gather(ratio_miRNA_merged, miRNA, value, first_miRNA:last_miRNA)%>%
        gather(ratio, ratio_value, ratio_murib_prevo)
data_nest <- group_by(data, miRNA, ratio) %>% nest()

miRNA29c <- as.data.frame(data_nest[[3]][[7]])
gr_df=as.data.frame(t(ma[c(1,2),]))
gr_df
miRNA29c_gr <-  merge (miRNA29c, gr_df, by= "mouse_id")
addicts <- subset (miRNA29c_gr, Grouping=="Addict")
cor(addicts$value, addicts$ratio)
plot(data_nest[[3]][[7]]$value, data_nest[[3]][[7]]$ratio_value)
text(data_nest[[3]][[7]]$value+0.1, data_nest[[3]][[7]]$ratio_value,data_nest[[3]][[7]]$mouse_id,cex=0.7)
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
title_p <- paste("Correlations between ratio Muribaculaceae/Prevotellacea \nand selected miRNAS",
                 "filtered zeros, NOT transformed\n")

# title_p <- paste("Correlations between ratio Muribaculaceae/Prevotellacea \nand all miRNAS",
#                  "filtered zeros, transformed\n")

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
ggsave (hm, file=paste0(out_dir, "heatmap_", "all_filtered_miRNA_NOT_transf_", "ratio_murib_prevo", extension_img),
        width = width_p, height = height_p, dpi=dpi_q)

