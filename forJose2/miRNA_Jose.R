#############################################################################
### Jose Espinosa-Carrasco CB-CRG Group. Nov 2020                         ###
#############################################################################
### Ratio analysis Ionas (toJose2)                                        ###
#############################################################################

library("vegan")

ra <- paste0(home_dir, "/git/food_addiction_analysis/forJose/relative_abundances_by_", taxon, ".txt")

ma=read.table(ra,
              sep="\t",
              skip=1,
              row.names=1,
              colClasses = "character")

m=read.table(ra,
             sep="\t",
             skip=3,
             row.names=1)

M=t(m)
rownames(M)=ma[1,]

gr=as.character(ma[2,])

bf <- paste0(home_dir, "/git/food_addiction_analysis/forJose/behavior.txt")
mb=read.table(bf, header=TRUE,sep="\t",row.names=1)

MB=mb[rownames(M),]

miRNA1 <- paste0(home_dir, "/git/food_addiction_analysis/forJose/miRNAs_1.txt")
miRNA2 <- paste0(home_dir, "/git/food_addiction_analysis/forJose/miRNAs_2.txt")

l1=read.table(miRNA1, header=FALSE, row.names=1, colClasses = "character")
l2=read.table(miRNA2, header=FALSE, row.names=1, colClasses = "character")
mi1=read.table(miRNA1, skip=1, row.names=1)
mi2=read.table(miRNA2, header=TRUE, row.names=1)
colnames(mi1)=l1[1,]
colnames(mi2)=l2[1,]

smi1=mi1[intersect(rownames(mi1),rownames(mi2)),]
smi2=mi2[intersect(rownames(mi1),rownames(mi2)),]
smi=cbind(smi1,smi2)


#reduce the raw matrices to samples present in microRNA matrix 
# Mm=M[intersect(colnames(mi),rownames(M)),]
all_mouse_miRNAs <- c(colnames(mi1), colnames(mi2))
Mm=M[intersect(all_mouse_miRNAs, rownames(M)),]
Mi=t(smi[,intersect(colnames(smi),rownames(Mm))])
Mb=mb[rownames(Mm),]

#get the microRNAs of interest:
mR=read.table("/Users/jaespinosa/git/food_addiction_analysis/forJose2/myResum_miRNAs.txt",
              skip=1,
              colClasses = "character")
mymi=which(colnames(Mi)%in%mR[,1])
names(mymi)=colnames(Mi)[mymi]


gr=as.character(ma[2,])
names(gr)=rownames(M)
gr2=gr[rownames(Mm)]
col=c(1:20)
col[which(gr2=="Non-Addict")]=2
col[which(gr2=="Addict")]=1


#behaviour explained by miRNAs:
vex2=rep(0,ncol(Mi))
names(vex2)=colnames(Mi)
for (i in 1:ncol(Mi)){
  cp=rda(Mb,Mi[,i])
  vex2[i]=cp$CCA$eig/(cp$CCA$eig+sum(cp$CA$eig))
}

#the top 3 explaining behavioural variance are from the "interesting" set, there are 5 of this set in the top 10:
names(vex2)[order(vex2,decreasing=TRUE)[1:10]]
# [1] "mmu-miR-3085-3p" "mmu-miR-876-5p"  "mmu-miR-665-3p"  "mmu-miR-1197-3p"
# [5] "mmu-miR-376b-3p" "mmu-miR-381-5p"  "mmu-miR-3072-3p" "mmu-miR-124-3p" 
# [9] "mmu-miR-6516-5p" "mmu-miR-935"
intersect(names(mymi),names(vex2)[order(vex2,decreasing=TRUE)[1:10]])
# [1] "mmu-miR-876-5p"  "mmu-miR-3085-3p" "mmu-miR-665-3p"  "mmu-miR-3072-3p"
# [5] "mmu-miR-124-3p"

#take the 3rd best:
vex2["mmu-miR-665-3p"]
# mmu-miR-665-3p 
# 0.2914104 

## mmu-miR-665-3p
cp=rda(Mb,Mi[,"mmu-miR-665-3p"])

cp$CCA$eig/(cp$CCA$eig+sum(cp$CA$eig))
cp$CA$eig[1]/(cp$CCA$eig+sum(cp$CA$eig))

png("../figures/RDAmiRNA.png",width=12,height=12,unit="cm",res=300)
plot(c(cp$CCA$u,cp$CCA$v),c(cp$CA$u[,1],cp$CA$v[,1]),type="n",main="Constrained PCA miRNA",xlab="RDA1 (mmu-miR-665-3p), 29.1%",ylab="PC1, 39.3%")
text(cp$CCA$v,cp$CA$v[,1],colnames(Mb),cex=0.5,col="grey")
text(cp$CCA$u,cp$CA$u[,1],rownames(Mb),col=col+1,cex=0.7)
legend("topleft",c("addicted","control"),col=c("red","green"),pch=c("B","B"),cex=0.7)
dev.off()

## mmu-miR-29c-3p
cp_29=rda(Mb,Mi[,"mmu-miR-29c-3p"])

cp_29$CCA$eig/(cp_29$CCA$eig+sum(cp_29$CA$eig))
cp_29$CA$eig[1]/(cp_29$CCA$eig+sum(cp_29$CA$eig))

png("../figures/RDAmiRNA.png",width=12,height=12,unit="cm",res=300)
plot(c(cp_29$CCA$u,cp_29$CCA$v),c(cp_29$CA$u[,1],cp_29$CA$v[,1]),type="n",main="Constrained PCA miRNA",xlab="RDA1 (mmu-miR-29c-3p), 18.4%",ylab="PC1, 49.4%")
text(cp_29$CCA$v,cp_29$CA$v[,1],colnames(Mb),cex=0.5,col="grey")
text(cp_29$CCA$u,cp_29$CA$u[,1],rownames(Mb),col=col+1,cex=0.7)
legend("topleft",c("addicted","control"),col=c("red","green"),pch=c("B","B"),cex=0.7)
dev.off()

#our microbiota ratio
colnames(Mm)[c(10,53)]
# [1] "Bacteroidales S24-7 group_unidentified"
# [2] "Prevotellaceae_UCG001"  

#the constrained PCA we had before on the full samples:
cpc=rda(MB,log(M[,10]/M[,53]))



## spearman correlation 
spe=rep(0,ncol(Mi))
names(spe)=colnames(Mi)
for (i in 1:ncol(Mi)){
  spe[i]=cor(log(Mi[,i]+0.5),log(Mm[,10]/Mm[,53]),method="spearman")
}

## top 10:
names(spe)[order(spe,decreasing=TRUE)[1:10]]
# [1] "mmu-miR-764-3p"  "mmu-miR-540-3p"  "mmu-miR-770-5p"  "mmu-miR-3552"   
# [5] "mmu-miR-544-3p"  "mmu-miR-3085-3p" "mmu-miR-337-3p"  "mmu-miR-377-5p" 
# [9] "mmu-miR-6976-3p" "mmu-miR-137-5p" 
intersect(names(mymi),names(spe)[order(spe,decreasing=TRUE)[1:10]])
# [1] "mmu-miR-3085-3p" "mmu-miR-544-3p" 
spe[order(spe,decreasing=TRUE)]
spe[order(spe,decreasing=TRUE)]['mmu-miR-29c-3p']

#this one is also decent (ranked 11):
spe["mmu-miR-665-3p"]
# mmu-miR-665-3p 
# 0.4466165 

#the three have "significant" rank correlation with our ratio
# (but wouldn't survive full multiple testing correction)
cor.test(log(Mi[,29]),log(Mm[,10]/Mm[,53]),method="spearman")
cor_665 = cor.test(log(Mi[,29]),log(Mm[,10]/Mm[,53]),method="spearman")
# p-value = 0.04985
cor_665$estimate
# rho = 0.4466165

## Obtain index of the column
grep('mmu-miR-29c-3p', colnames(Mi))
colnames(Mi)[55]
cor.test(log(Mi[,55]),log(Mm[,10]/Mm[,53]), method="spearman")
cor_29c = cor.test(log(Mi[,55]),log(Mm[,10]/Mm[,53]), method="spearman")
# p-value = 0.40
cor_29c$estimate
# rho = 0.1954887

png("../figures/corrMratioVmiRNA.png",width=12,height=12,unit="cm",res=300)
plot((Mi[,"mmu-miR-665-3p"]),(Mm[,10]/Mm[,53]),type="n",main="Microbiota ratio vs miRNA",xlab="mmu-miR-665-3p abundance",ylab="Muribaculaceae/Prevotellaceae",log="xy")
text((Mi[,"mmu-miR-665-3p"]),(Mm[,10]/Mm[,53]),rownames(Mb),col=col+1,cex=0.7)
dev.off()

#to get the regression line for the correlation plot:
res=lm(log(Mm[,10]/Mm[,53])~log(Mi[,"mmu-miR-665-3p"]))

library(vioplot)
# setwd("/Users/jaespinosa/git/food_addiction_analysis/forJose2")
png("../figures/panels_665.png",width=25,height=25,unit="cm",res=300)
par(mfrow=c(2,2))

vioplot(M[which(gr=="Addict"),10]/M[which(gr=="Addict"),53],M[which(gr=="Non-Addict"),10]/M[which(gr=="Non-Addict"),53],col=c("red","green"),ylog=TRUE,ylab="Muribaculaceae/Prevotellaceae",names=c("addicted","control"),main="Microbiota ratio in groups")

plot(log(Mi[,"mmu-miR-665-3p"]),log(Mm[,10]/Mm[,53]),type="n",main="Microbiota ratio vs miRNA",xlab="log(mmu-miR-665-3p abundance)",ylab="log(Muribaculaceae/Prevotellaceae)")
points(log(Mi[,"mmu-miR-665-3p"]),log(Mm[,10]/Mm[,53]),col=col+1,pch=19)
legend("topleft",c("control","addicted"),col=c("green","red"),pch=19)
text(6,4.2, paste0("rho=", round(cor_665$estimate,3)),  cex=0.8)
text(6,3.9, paste0("pvalue=", round(cor_665$p.value,4)),  cex=0.8)
abline(res)

plot(c(cpc$CCA$u,cpc$CCA$v),c(cpc$CA$u[,1],cpc$CA$v[,1]),type="n",main="Constrained PCA microbiota ratio",xlab="RDA1 (log(Muribaculaceae/Prevotellaceae)), 24.6%",ylab="PC1, 42.6%")
text(cpc$CCA$v,cpc$CA$v[,1],colnames(MB),cex=0.7)
text(cpc$CCA$u,cpc$CA$u[,1],rownames(MB),col=c(rep(3,11),rep(2,13)),cex=0.7,pch=20)

plot(c(cp$CCA$u,cp$CCA$v),c(cp$CA$u[,1],cp$CA$v[,1]),type="n",main="Constrained PCA miRNA",xlab="RDA1 (mmu-miR-665-3p), 29.1%",ylab="PC1, 39.3%")
text(cp$CCA$v,cp$CA$v[,1],colnames(Mb),cex=0.7)
text(cp$CCA$u,cp$CA$u[,1],rownames(Mb),col=col+1,cex=0.7,pch=20)

dev.off()

#to get the regression line for the correlation plot:
res2=lm(log(Mm[,10]/Mm[,53])~log(Mi[,"mmu-miR-29c-3p"]))

## mmu-miR-29c-3p
png("../figures/panels_29c.png",width=25,height=25,unit="cm",res=300)
par(mfrow=c(2,2))

vioplot(M[which(gr=="Addict"),10]/M[which(gr=="Addict"),53],
        M[which(gr=="Non-Addict"),10]/M[which(gr=="Non-Addict"),53],
        # cex.lab=2,
        cex.names=1,
        # cex.axis=1.2, 
        cex.main=1.2,
        col=c("red","green"),
        ylog=TRUE, ylab="Muribaculaceae/Prevotellaceae",
        names=c("Addicted","Control"),
        main="Microbiota ratio in groups")

plot(log(Mi[,"mmu-miR-29c-3p"]),
     log(Mm[,10]/Mm[,53]),type="n",
     main="Microbiota ratio vs miRNA",
     xlab="log(mmu-miR-29c-3p abundance)",
     ylab="log(Muribaculaceae/Prevotellaceae)")
points(log(Mi[,"mmu-miR-29c-3p"]),log(Mm[,10]/Mm[,53]),col=col+1, pch=19)
text(bip$sites[,1],bip$sites[,2],rownames(Mb),col=col+1,cex=0.7)
legend("topleft",c("control","addicted"),col=c("green","red"), pch=19)
# text(bip$sites[,1],bip$sites[,2], rownames(Mb), col=col+1, cex=0.7)
text(8.45,4, paste0("rho=", round(cor_29c$estimate,3)),  cex=0.8)
text(8.45,3.4, paste0("pvalue=", round(cor_29c$p.value,3)),  cex=0.8)
abline(res2)

plot(c(cpc$CCA$u,cpc$CCA$v),
     c(cpc$CA$u[,1],cpc$CA$v[,1]),
     type="n",main="Constrained PCA microbiota ratio",
     xlab="RDA1 (log(Muribaculaceae/Prevotellaceae)), 24.6%",ylab="PC1, 42.6%")
text(cpc$CCA$v,cpc$CA$v[,1],colnames(MB),cex=0.7)
text(cpc$CCA$u,cpc$CA$u[,1],rownames(MB),col=c(rep(3,11),rep(2,13)),cex=0.7,pch=20)

plot(c(cp_29$CCA$u, cp_29$CCA$v), c(cp_29$CA$u[,1],cp_29$CA$v[,1]),type="n",
     main="Constrained PCA miRNA",
     xlab="RDA1 (mmu-miR-29c-3p), 18.4%",
     ylab="PC1, 49.4%")
text(cp_29$CCA$v,cp_29$CA$v[,1],colnames(Mb),cex=0.7)
text(cp_29$CCA$u,cp_29$CA$u[,1],rownames(Mb),col=col+1,cex=0.7,pch=20)

dev.off()

#combine ratio and microRNA in a single plot


cons=t(rbind(log(Mm[,10]/Mm[,53]),Mi[,"mmu-miR-665-3p"]))
colnames(cons)=c("ratio","miRNA")

cc=rda(Mb,cons)
cc$CCA$eig/(cc$CCA$eig+sum(cc$CA$eig))
cc$CA$eig[1]/(sum(cc$CCA$eig)+sum(cc$CA$eig))
bip=plot(cc)

png("../figures/biplot.png",width=12,height=12,unit="cm",res=300)
plot(bip$sites[,1],bip$sites[,2],type="n",xlim=c(-30,30),asp=1,xlab="RDA1 (38.2%)",ylab="RDA2 (4.8%)",main="Constrained PCA biplot")
text(bip$sites[,1],bip$sites[,2],rownames(Mb),col=col+1,cex=0.7)
text(bip$species[,1],bip$species[,2],colnames(Mb),cex=0.5)
arrows(0,0,bip$biplot[,1]*attr(bip$biplot,"arrow.mul"),bip$biplot[,2]*attr(bip$biplot,"arrow.mul"),length=0.1)
text(bip$biplot[,1]*attr(bip$biplot,"arrow.mul"),bip$biplot[,2]*attr(bip$biplot,"arrow.mul"),colnames(cons),pos=3,cex=0.7)
dev.off()

## regression only addicts
plot(log(Mi[,"mmu-miR-665-3p"]),log(Mm[,10]/Mm[,53]),type="n",main="Microbiota ratio vs miRNA",xlab="log(mmu-miR-665-3p abundance)",ylab="log(Muribaculaceae/Prevotellaceae)")
points(log(Mi[,"mmu-miR-665-3p"]),log(Mm[,10]/Mm[,53]),col=col+1,pch=19)
legend("topleft",c("control","addicted"),col=c("green","red"),pch=19)
abline(res)

cor.test(log(Mi[,29]),log(Mm[,10]/Mm[,53]),method="spearman")
cor_665 = cor.test(log(Mi[,29]),log(Mm[,10]/Mm[,53]),method="spearman")
# p-value = 0.04985
cor_665$estimate
# rho = 0.4466165
rownames(Mi)
length(Mi[,1])

group = c(1:20)
group[which(gr2=="Non-Addict")] = "Non-Addict"
group[which(gr2=="Addict")] = "Addict"

df_miRNAs <- cbind(log(Mi[,29]), log(Mi[,55]), log(Mm[,10]/Mm[,53]))
df_miRNAs <- as.data.frame(df_miRNAs)
colnames(df_miRNAs) <- c("miR-665-3p", "miR-29c-3p", "ratio")
df_miRNAs$group <- group

df_miRNAs_addicts <- subset (df_miRNAs, group=="Addict")
cor.test(df_miRNAs_addicts$`miR-29c-3p`, df_miRNAs_addicts$ratio, method = 'spearman')

res3=lm(df_miRNAs_addicts$ratio~df_miRNAs_addicts$`miR-29c-3p`)

plot(df_miRNAs_addicts$`miR-29c-3p`, df_miRNAs_addicts$ratio, type="n",main="Microbiota ratio vs miRNA",xlab="log(mmu-miR-29c-3p abundance)",ylab="log(Muribaculaceae/Prevotellaceae)")
points(df_miRNAs_addicts$`miR-29c-3p`,df_miRNAs_addicts$ratio,col="red",pch=19)
legend("topleft",c("addicted"),col=c("red"),pch=19)
abline(res3)
