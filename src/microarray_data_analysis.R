#install.packages(c("limma","GEOquery","Biobase","pheatmap",
#                   "ggplot2","gplots","reshape2","plyr"))
#install.packages("Biobase")
#install.packages("GEOquery")
#install.packages("limma")
#install.packages("gplots")
#install.packages("ggplot2")
#install.packages("plyr")
#install.packages("pheatmap")
library(limma)
library(Biobase)
library(GEOquery)
library(pheatmap)
library(gplots)
library(ggplot2)
library(plyr)
package.version("limma")

setwd("/home/alireza/Program Data/Computational Biology")
gset = getGEO("GSE52509",GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data/")
gset=getGEO(filename=file.choose(),GSEMatrix = TRUE, AnnotGPL = TRUE)

class(gset)

gset=gset[[1]]
print(gset)
class(gset)
colnames(gset)
sml=c(rep("smoke_4",3),rep("control_4",3),rep("smoke_6",3),rep("control_6",3))
length(sml)
sml <- factor(sml)
levels(sml)
class(sml)

ex<-exprs(gset)
class(ex)
dim(ex)
hist(ex[,1])
#how to know is normal
max(ex)
min(ex)
# if it is not normal run
ex<-log2(ex)
hist(ex[,1])
exprs(gset)<-ex
pdf("results/boxplots.pdf",width=6)
boxplot(ex,outline=FALSE)
dev.off()
#if data was not normal
#ex<-normalizeQuantiles(ex)
#exprs(gse)<-ex
x<-normalizeQuantiles(ex)
pdf("results/boxplot_qnormal.pdf")
boxplot(x)
dev.off()

#Plot correlation HeatMap

library(pheatmap)
pdf("results/corr_heatmap.pdf")
pheatmap(cor(ex))
dev.off()

#corr_ex<-cor(ex)
#corr_ex[1:4,1:4]
pdf("results/corr_heatmap1.pdf",width = 10,height = 10)
pheatmap(cor(ex),labels_row = sml,labels_col = sml)
dev.off()
# Principal component analysis
pca<-prcomp(ex)
pdf("results/PCA.pdf")
plot(pca)
dev.off()
names(pca)
pca$sdev
colnames(pca$x)
pdf("results/PCAxx.pdf")
plot(pca$x[,1:2])
dev.off()

ex_scale=t(scale(t(ex),scale=F))
mean(ex_scale[1,])
pca2<-prcomp(ex_scale)
pdf("results/PCA_scales.pdf")
plot(pca2)
plot(pca2$x[,1:2])
dev.off()

plot(pca2$r)
pc.sample<-data.frame(pca2$r[,1:3],Group=sml)
head(pc.sample)
pdf("results/pca_samples.pdf")
ggplot(pc.sample,aes(PC1,PC2,color=Group))+geom_point(size=3)+theme_bw()
dev.off()

pca2<-prcomp(t(ex_scale))
pc.sample<-data.frame(pca2$x[,1:3],Group=sml)
head(pc.sample)
pdf("results/pca_genes.pdf")
ggplot(pc.sample,aes(PC1,PC2,color=Group))+geom_point(size=3)+theme_bw()
dev.off()

sml <- factor(sml)
levels(sml)
sml
gset$description <- sml
design <- model.matrix(~ description + 0, gset) #112
colnames(design) <- levels(sml)
head(design)
design

gset

fit <- lmFit(gset, design)
#cont.matrix <- makeContrasts(mutant-wild, levels=design)
cont.matrix <- makeContrasts(smoke_4-control_4, levels=design)
#Differential express genes: ژن هایی تمایز بیان دارند.

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01) #119
# آماره بیزی را حساب میکند
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
#برای فیت۲۵۰ ۲ تای اولی را که بهترین 
#adjusted p-value
#را دارند انتخاب میکند . 
#یعنی ژن هایی که بین اسموک۴ و 
#کنترل۴ تمایز بیان دارن هریک یک بی ولیو ای داند و ۲۵۰ تای اول را انتخاب میکند
# ولی بی ولیوی خام به درد ما نمیخورند و باید ادجاست شوند
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
colnames(tT)
head(ex)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC",
                          "Gene.symbol","Gene.title"))
#سرستون های دلخواه را انتخاب میکنیم
write.table(tT, "results/ptA_ptB.txt", row.names=F, sep="\t",quote = F)

ms.up=subset(tT,select=c("ID","adj.P.Val","logFC","Gene.symbol"))
ms.up=subset(tT,tT$logFC>1 & tT$adj.P.Val<0.05,select=c("ID","adj.P.Val","logFC",
                                                        "Gene.symbol"))
#آنهایی که بیانشان افزایش بیدا کرده
#up regulate

write.table(ms.up,file= "results/ptA_ptB.txt", row.names=F, sep="\t",quote = F)

ms.up.genenames<- sub("///.*","",ms.up$Gene.symbol)
ms.up.genenames<- ms.up.genenames[ms.up.genenames!= ""]
ms.up.genenames<- strsplit2(ms.up.genenames,"///")
ms.up.genenames<- unique(ms.up.genenames)

ms.up.genenames
ms.up$Gene.symbol

write.table(ms.up.genenames,"results/up_ptA_ptB.txt", quote = F,row.names = F, col.names = F)
