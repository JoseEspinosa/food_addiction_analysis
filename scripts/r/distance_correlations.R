
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

m_microbiome=read.table(ra,sep="\t",skip=3,row.names=1)
colnames(m_microbiome)=ma[1,]

cor_m_microbiome <- cor(m_microbiome, method="pearson")
cor_m
dist_m_microbiome <- 1-cor_m

##############
## Behavior
bf <- paste0(home_dir, "/git/food_addiction_analysis/forJose/behavior.txt")
mb=read.table(bf,header=TRUE,sep="\t",row.names=1)
mb_t <- t(mb)

## reorder to match microbiome columns (mice)
col_names_microb <- colnames(m_microbiome)
mb_t_ord <- mb_t [, col_names_microb]
cor_m_beh <- cor(mb_t_ord, method="pearson")
dist_m_beh <- 1 - cor_m_beh

colnames (dist_m_beh)
colnames (dist_m_microbiome)

str(dist_m_beh)
str(dist_m_microbiome)

cor(c(dist_m_beh), c(dist_m_beh))
cor_distance_beh_microb = cor(c(dist_m_beh), c(dist_m_microbiome), method="pearson")

#########################################
# Filter microbiome abundances like Ionas
M=t(m_microbiome)
rownames(M)=ma[1,]

## Any genera with zeros is removed
test=apply(M,2,min)
test
M_f = M[,-which(test==0)]
mb_t_f_ord <- t(M_f) [, col_names_microb] 
mb_t_f_ord

cor_m_microbiome_filt <- cor(mb_t_f_ord, method="pearson")
cor_m_microbiome_filt
dist_m_microbiome_filt <- 1-cor_m_microbiome_filt

cor_distance_beh_microb_filt = cor(c(dist_m_beh), c(dist_m_microbiome_filt), method="pearson")





