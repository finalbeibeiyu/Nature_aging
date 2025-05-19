getwd()
setwd("/mnt/data/home/tycloud/bone")
options(stringsAsFactors = F)
#加载R包
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
#读取数据
#1、批量读取单细胞的数据
dir_name=c('Old','FBMSC-EVs','ABPC-EVs')
datalist=list()
for (i in 1:length(dir_name)){
  dir.10x = paste0('bone_marrow/',dir_name[i])
  my.data <- Read10X(data.dir = dir.10x) 
  datalist[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i], 
                                   #每个基因至少在3个细胞中表达，每一个细胞至少有250个基因表达
                                   min.cells = 3, min.features = 200)
}
#修改名称
names(datalist)=dir_name
#2、细胞质控####
# 批量计算线粒体和rRNA占比
##############过滤红细胞基因
## 将红细胞所有marker基因输入到HB_genes向量

mt_genes <- rownames(datalist$Old)[grep("^MT[^(p)]", rownames(datalist$Old),ignore.case = T)]
mt_genes

for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentGroupFeatureSet(sce,pattern = "^MT")# 计算线粒体占比 ^MT- human
  datalist[[i]] <- sce
  rm(sce)}


#质控前的
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')


#样本合并
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
#统计每一个样本的个数
raw_count <- table(sce@meta.data$orig.ident)
table(sce@meta.data$percent.mt)


pearplot_befor1<-VlnPlot(sce,
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                         pt.size = 0.01, 
                         ncol = 4,
                         cols=c('#AB3282', '#23452F', '#BD956A', '#CCC9E6', '#585658', '#58A4C3'))
pearplot_befor1
ggsave(filename = 'QC_before1.pdf',plot = pearplot_befor1,he=7,wi=15)
rm(sce)


#过滤



datalist1 <- lapply(X = datalist, FUN = function(x) {
  x<-subset(x,subset = nFeature_RNA > 300
            & nFeature_RNA < 10000  
            & nCount_RNA < quantile(nCount_RNA, 0.98)
            &percent.mt < 5)
})

#合并数据
sce <- merge(datalist1[[1]],y=datalist1[2:length(datalist1)])


clean_count <- table(sce@meta.data$orig.ident)
table(sce@meta.data$orig.ident)


#过滤前后样本细胞数据的统计
summary_cells <- as.data.frame(cbind(raw_count,clean_count))
counts <- rbind(as.data.frame(cbind(summary_cells[,1],rep("raw",each = length(summary_cells[,1])))),
                as.data.frame(cbind(summary_cells[,2],rep("clean",each = length(summary_cells[,2])))))
counts$Group <- rep(rownames(summary_cells),times =2)
colnames(counts)<- c("count","Stat","Group")
counts[,1] <- as.numeric(counts[,1])
library(ggsce)
counts$Stat <- factor(counts$Stat, levels=c("raw", "clean"), ordered=TRUE)
fit_cell_count <- ggplot(data =counts, mapping = aes(x = Group, y=count))+ 
  geom_bar(aes(fill = Stat),stat = 'identity', position = 'dodge') + 
  scale_fill_jama() +
  theme(text=element_text(size=25),
        legend.title=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black",size = 1),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.background = (element_rect(colour= "white",fill = "white")))

fit_cell_count
ggsave(filename = 'fit_cell_count.pdf',plot = fit_cell_count,width = 8,height = 5)
#质控后的小提琴图
pearplot_after1 <- VlnPlot(sce,
                           features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                           pt.size = 0.01,
                           ncol = 4,
                           cols=c('#AB3282', '#23452F', '#BD956A', '#CCC9E6', '#585658', '#58A4C3'))
pearplot_after1
ggsave(filename = 'QC_after1.pdf',plot = pearplot_after1,he=6,wi=12)



#保存datalist文件

#3、数据预处理####
#合并数据

sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", 
                            nfeatures = 2000,
                            mean.cutoff=c(0.0125,3),
                            dispersion.cutoff =c(1.5,Inf))
### 可视化前20个高变基因
top20 <- head(VariableFeatures(sce), 20)
plot1 <- VariableFeaturePlot(sce)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE, size=3.0)

feat_20 <- CombinePlots(plots = list(plot1, plot2),legend="bottom")
feat_20
ggsave(filename = 'feat_20.pdf',plot = feat_20,he=10,wi=15)
#ScaleData
scale.genes <-  rownames(sce)
sce <- ScaleData(sce, features = scale.genes)
#样本的分组
meta1<-data.frame(matrix(nrow=length(sce@meta.data$orig.ident), ncol=2)) 
colnames(meta1)=c('Sample','Group')
meta1$Sample=sce@meta.data$orig.ident
table(sce$orig.ident)
### Group Tumor 为原发性肿瘤；Normal：正常

sce <- AddMetaData(sce, meta1$Group,col.name = "Group")

table(sce@meta.data$Group)


sce$Group <- factor(sce$Group,levels = c("Old",
                                         "FBMSC",
                                         "ABPC"))

sce <- RunPCA(sce, features = VariableFeatures(sce)) 
dimplot1 <- DimPlot(sce, reduction = "pca") 
elbowplot1 <- ElbowPlot(sce, ndims=50, reduction="pca") 
sc_pca <- dimplot1+elbowplot1
sc_pca
ggsave(filename = 'sc_pca.pdf',plot = sc_pca,he=10,wi=15)
#可视化前2个PC的top20个基因
VizDimLoadings(sce, dims = 1:10, nfeatures = 20, reduction = "pca")

before_pca <- DimPlot(sce, reduction = "pca", group.by = "Sample",cols=my36colors)
ggsave('before_pca.pdf',before_pca,he=4,wi=5)
ElbowPlot(sce, ndims=50)

# 不一定要聚类，但modelHomotyp、
#####################4.去双细胞##################


library(scDblFinder)
sce1 <- as.SingleCellExperiment(sce)
sce1 <- scDblFinder(sce1,dbr=8*1e-6)
table(sce1$scDblFinder.score)
#查看双细胞得分
table(sce1$scDblFinder.class)


pbmc <- as.Seurat(sce1)
dob <- DimPlot(pbmc, reduction = "tsne", group.by = "scDblFinder.class",cols = c('#E39A35', '#E63863'))
ggsave('4_双细胞去除.pdf',dob,he=4,wi=5)
table(pbmc$scDblFinder.class)
#先把sce对象转换回Seurat对象
sce1=subset(x=pbmc,subset=scDblFinder.class=='singlet')
#singlet doublet 
#28558     852 

sce$scDblFinder.class <- sce1$scDblFinder.class
sce$scDblFinder.score <- sce1$scDblFinder.score
table(sce$scDblFinder.class)
sce=subset(x=sce,subset=scDblFinder.class=='singlet')


dim.usage<- 15
library(harmony)
sce <- RunHarmony(sce, group.by.vars = "Group")
sce <- FindNeighbors(sce, dims = 1:dim.usage, reduction = "harmony") 
table(sce$Group)
ElbowPlot(sce,reduction = 'harmony')


library(clustree)

sce <- FindClusters(sce,resolution = c(0.2,0.3,0.4,0.5,0.6))

q <- clustree(sce@meta.data, prefix = "RNA_snn_res.")+ ggsce::scale_color_jco()
ggsave('树状图.pdf',q,he=10,wi=15)


sce <- RunUMAP(sce, dims = 1:15,reduction = 'harmony')
sce <- RunTSNE(sce, 
               dims=1:15, 
               reduction="harmony",
               perplexity=30,
               max_iter=1000)
# 去批次成功
after_pca <- DimPlot(sce, reduction = "pca", group.by = "Group",cols=my36colors)
ggsave('after_pca.pdf',after_pca,he=4,wi=5)
table(sce$Group)
#颜色
Idents(sce) <- sce$RNA_snn_res.0.3


############### 继续做柱状图###################

#可视化
library(SCP)
library(BiocParallel)
library(stats)
library(ggsce)
color_cluster <- pal_jama()(10)


f <- CellDimPlot(sce, 
                 group.by = "RNA_snn_res.0.3", 
                 reduction = "umap", 
                 label =T,
                 #hex = TRUE, 
                 hex.bins = 30, 
                 palcolor=my36colors)

ggsave('RNA0.3UMAP图.pdf',f,he=6,wi=8)




############细胞周期基因去除#########
Idents(sce) <- sce$seurat_clusters
str(cc.genes)
library(scuttle)
#install.packGroups(scuttle)
mouse_cell_cycle_genes <- readRDS("mouse_cell_cycle_genes.rds")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sce <- CellCycleScoring(sce, 
                        s.features =s.genes, 
                        g2m.features = g2m.genes)

sce<- RunPCA(sce, features = c(s.genes, g2m.genes))

DimPlot(sce, reduction = "pca", group.by = "Phase", shape.by = "Group",
        cols=my36colors,raster = FALSE)

pdf('细胞周期评分.pdf',he=5,wi=5)

plot(sce$S.Score,sce$G2M.Score,
     col=factor(sce$Phase),
     main="CellCycleScoring")
legend("topleft",inset=.05,
       title = "cell cycle",  
       c("G1","S","G2M"), pch = c(1),col=c("black","green","red"))
dev.off()


pdf('细胞周期基因.pdf',he=6,wi=7)
RidgePlot(sce, 
          features = c("PCNA", "TOP2A", "MCM6", "MKI67"), 
          cols = my36colors,
          ncol = 2)
dev.off()

###消除影响
sce<- ScaleData(sce,
                vars.to.regress = c("S.Score", "G2M.Score"),
                features = rownames(sce))

p4=VlnPlot(sce, features = c("S.Score", "G2M.Score"), group.by = "Group", 
           ncol = 2, pt.size = 0.1)


load("sce20231228.RData")



Idents(sce) <- sce$RNA_snn_res.0.3
#marker基因的筛选################
#需要修改#############
#寻找差异基因时的差异倍数
Logfc = 0.5
#差异基因时最小的表达比例
Minpct = 0.15
DefaultAssay(sce) <- "RNA"
sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, 
                              min.pct = Minpct,only.pos = T)
sce.markers["pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
length(unique(sce.markers$gene))
head(sce.markers)

write.table(sce.markers,'sce_marker_gene_0.4.txt',quote = F,row.names = F,sep='\t')

DimPlot(sce,
        group.by ="RNA_snn_res.0.4", 
        reduction="umap",
        #reduction="tsne",
        label = "T", 
        pt.size = 0.2,
        label.size = 5)


##############################注释###################  
library(SingleR)
library(celldex)
library(BiocParallel)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scRNAtoolVis)
load("sce_before_anno.RData")

####################anno##################
hpca.se=get(load("ref_Human_all.RData"))
#hpca.se <- HumanPrimaryCellAtlasData()
table(sce$Group)

#获取基因的表达谱的count数据
testdata <- GetAssayData(sce, slot="data")
#获取聚类的亚群
clusters <- sce@meta.data$RNA_snn_res.0.3
pred.sce <- SingleR(test =  testdata, 
                    ref = hpca.se, 
                    labels = hpca.se$label.fine,
                    method = "clusters",
                    clusters = clusters, 
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")

plotScoreHeatmap(pred.sce)

celltype = data.frame(ClusterID=rownames(pred.sce), celltype=pred.sce$labels
                      , stringsAsFactors = F)

celltype
write.table(celltype,'celltype.txt',quote = F,sep = '\t',row.names = F)
#手工注释
celltype=read.table("celltype.txt", header=T, sep="\t", check.names=F)
celltype = data.frame(celltype, stringsAsFactors = F)

sce$cell_type=sce$RNA_snn_res.0.3
sce$seurat_clusters <- sce$cell_type
for (i in 1:nrow(celltype)){
  sce$cell_type=gsub(paste0('^',celltype[i,1],'$'),celltype[i,2],sce$cell_type)
}
unique(sce$cell_type)
table(sce1$cell_type)

length(table(sce$cell_type))

table(sce1$Group)
Idents(sce) <- sce$cell_type
###去除基质细胞################
sce1 <- sce[, sce$cell_type != "Stromal cells"]

pdf("1_TSne图注释.pdf", width=5, height=4)
CellDimPlot(sce1, 
            group.by = "cell_type", 
            reduction = "Tsne", 
            label =T,
            #hex = TRUE, 
            hex.bins = 30, 
            palcolor=color_cluster)
dev.off()

sce <- sce1
library(parallel)
mc.cores <- detectCores()
options(mc.cores = 24)
register(SnowParam(workers = 10, progressbar = TRUE))
library(SCP)
library(BiocParallel)
library(stats)

sce <-  RunDEtest(srt = sce, group_by = "cell_type", fc.threshold = 1, min.pct = 0.15,only.pos = FALSE)

ggsave("1_桑基图.pdf", height = 5, width = 6)
CellStatPlot(sce1,   
             stat.by = "cell_type", 
             group.by = "Group",
             individual = F, 
             plot_type = "trend",
             palcolor=color_cluster)


load('sce去除基质细胞.RData')

table(sce$cell_type)

myMarkers <- c("LILRA4",	"PLAC8",	"IRF8","GZMB",#Dendritic cells
               "IL7R","CD3E","CD3D","CD3G",
               "MZB1",	"JCHAIN",	"TNFRSF17",	"DERL3",#	Plasma 
               "VPREB3",	"CD79A"	,"CD79B", #B-cell 
               "HBB",	"HBA1",	"HBQ1",	"AHSP", #	Erythroid cells
               "SRGN",	"MPO",	"AZU1",	"PRTN3",	#	Granulocytic cells
               "VCAN",	"S100A12",	"LGALS3",	"ITGAM",	"CD14"	#	Monocytes
               )	#T细胞


color_cluster <- pal_jama()(10)
CellDimPlot(sce,
            group.by = "cell_type", 
           split.by = "Group", 
            reduction = "Tsne", 
            add_density = T,
            show_stat = FALSE,
            cells.highlight = TRUE, 
            theme_use = "theme_blank", 
            palcolor=color_cluster) 
ggsave("2_分组TSNE图.pdf", height = 4, width = 10)


sce$cell_type <- factor(sce$cell_type,levels = c("B cells","Dendritic cells","Erythroid cells",
                                                 "Granulocytic cells","Monocyte","Plasma cells","T cells"))

myMarkers2 <- c("LILRA4",#Dendritic cells
               "CD3E",
               "MZB1",	#	Plasma 
               	"CD79A"	, #B-cell 
               "HBA1",	 #	Erythroid cells
               	"MPO",	#	Granulocytic cells
               "VCAN"	#	Monocytes
)	#T细胞
pdf("3_点热图.pdf",width = 10, height = 5)
jjDotPlot(object = sce,
          gene  =  myMarkers2,
          anno = F,
          plot.margin = c(3,1,1,1),
          point.geom = T,
          id = 'celltype',
          tile.geom = T,
          x.text.angle=60,
          dot.col=c('White','#79AF97'))
dev.off()

pdf("3_umapmarker.pdf",width = 12, height = 3)
FeatureCornerAxes(object = sce,
                  reduction = 'tsne',
                  show.legend = T,
                  groupFacet = NULL,
                  relLength = 0.5,
                  relDist = 0.2,
                  features = myMarkers2)
dev.off()
##################semayo##################################

library(AUCell)
library(clusterProfiler)
library(ggplot2)
library(Seurat)
library(GSEABase)

library(readxl)
SenMayo <- openxlsx::read.xlsx("./SenMayo.xlsx")
unique(SenMayo$gene)
a <- SenMayo[SenMayo$gene %in% rownames(sce),]
# AUCell_buildRankings
cells_rankings <- AUCell_buildRankings(sce@assays$RNA@counts,featureType = "genes",plotStats=TRUE,splitByBlocks=TRUE)  
## 制作基因集
geneSets <- a
extraGeneSets <- c(GeneSet(sample(geneSets),setName="SenMayo_gene"))

geneSets1 <- GeneSetCollection(extraGeneSets)
names(geneSets1) 
cells_AUC <- AUCell_calcAUC(geneSets1, cells_rankings)
cells_AUC
pdf("小胶质细胞sce后cells_assignment.pdf",width = 8, height = 6)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=T, assignCells=TRUE) 
dev.off()
thresholds <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE)
getThresholdSelected(cells_assignment)
cells_assignment$ABPC_gene$aucThr$selected <- cells_assignment$ABPC_gene$aucThr$thresholds
cells_assignment$ABPC_gene$aucThr$thresholds

##set gene set of interest here for plotsceng
geneSet <- "SenMayo_gene"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
sce$AUC <- aucs
df<- data.frame(sce@meta.data, sce@reductions[["tsne"]]@cell.embeddings)
####################提取特定细胞#######################
meta1<-data.frame(matrix(nrow=length(sce@meta.data$AUC), ncol=2)) 
colnames(meta1)=c('AUC1','AUC2')
meta1$AUC1=sce@meta.data$AUC
unique(meta1$AUC1)
meta1$AUC2 <- ifelse(meta1$AUC1>0.03455305,"SenMayo","No_SenMayo")
sce <- AddMetaData(sce, meta1$AUC2,col.name = "AUC2")
table(sce$AUC2)

colnames(df)

class_avg <- df %>%
  group_by(cell_type) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )


FeatureDimPlot(
 sce, features = paste0("AUC"), 
  reduction = "tsne")

ggsave("3_AUC评分UMAP图.pdf", height = 4, width = 5)

FeatureStatPlot(sce, 
                stat.by = "AUC", 
                group.by = "cell_type", 
                add_trend = TRUE,
                palcolor=color_cluster, 
                add_box = TRUE,
                comparisons = list(c("Monocyte","B cells"),
                                   c( "Monocyte","Dendritic cells"),
                                   c( "Monocyte","Erythroid cells"),
                                   c( "Monocyte","Granulocytic cells"),
                                   c( "Monocyte","Plasma cells"),
                                   c( "Monocyte","T cells")), 
                bg.by = "cell_type")
ggsave("3_AUC评分细胞比较.pdf", height = 4, width = 8)


FeatureStatPlot(sce, 
                stat.by = "AUC", 
                group.by = "Group", 
                #split.by = "cell_type",
                add_trend = TRUE,
                palcolor=color_cluster, 
                add_box = TRUE,
                comparisons = list(c("ABPC","FBMSC"),
                                   c( "ABPC","Old")), 
                bg.by = "Group")
ggsave("3_AUC评分分组比较.pdf", height = 4, width = 5)


sce$G2M.Score
FeatureStatPlot(sce, 
                stat.by = "S.Score", 
                group.by = "Group", 
                #split.by = "cell_type",
                add_trend = TRUE,
                palcolor=color_cluster, 
                add_box = TRUE,
                comparisons = list(c("ABPC","FBMSC"),
                                   c( "ABPC","Old")), 
                bg.by = "Group")
ggsave("3_S评分分组比较.pdf", height = 4, width = 5)

table(sce$Phase)
CellStatPlot(sce, 
             #split.by = "cell_type",
             stat.by = c("AUC2"), 
             plot_type = "trend",
             group.by = "Group",
             label = T,
             palcolor=c("#79AF97","#354E6B"))
ggsave("3_AUC评分柱状图比较.pdf", height = 4, width = 4)

table(sce$Group)
M <- sce@tools$DEtest_cell_type$AllMarkers_wilcox
write.table(M,'celltypemarker.txt',quote = F,sep = '\t',row.names = F)
#################################差异基因#################
Mono <- subset(sce,cell_type=="Monocyte")
Mono <- RunDEtest(srt = Mono, 
                 group_by = "Group",
                 group1="Old",group2="ABPC",
                 fc.threshold =1, only.pos = FALSE)
VolcanoPlot(Mono, group_by = "custom")
ggsave("4_ABPC火山图.pdf", height = 4, width = 4)


Q <- Mono@tools$DEtest_custom$AllMarkers_wilcox

library(clusterProfiler)
library(org.Hs.eg.db)
library(GseaVis)
GseaVis()
# load diff
Q$avg_log2FC
diff <- Q %>%
  arrange(desc(avg_log2FC))

genelist <- diff$avg_log2FC
names(genelist) <- diff$gene

# check
head(genelist,3)
#    Aplnr    Foxf1     Bmp5
# 13.45176 13.35322 12.02845

# load tpm

#富集分析
ego <- gseGO(geneList     = genelist,
             OrgDb        = org.Hs.eg.db,
             ont          = "BP",
             keyType      = 'SYMBOL',
             minGSSize    = 5,
             maxGSSize    = 500,
             pvalueCutoff = 1,
             verbose      = FALSE)
b <- ego@result
write.table(b,"ABPCBP.txt",sep = "\t",quote = F,row.names = F)


#########多个GSEA

ego
terms <- c('GO:0006955',
           'GO:0006954',
           'GO:0045333',
           'GO:0008219',"GO:0012501","GO:0000723",
           "GO:0006281","GO:0051276")
pdf("top6单独NES比例图.pdf",width = 8, height = 10)
lapply(terms, function(x){
  gseaNb(object = ego,
         geneSetID = x,
         subPlot=2,
         addPval = T,
         pvalX = 0.5,pvalY = 0.7,
         pCol = 'black',
         pHjust = 0)
}) -> gseaList

# combine
cowplot::plot_grid(plotlist = gseaList,ncol = 2,align = 'hv')
dev.off()

pdf("top6单独NES比例图新样式.pdf",width = 8, height = 10)

######多样式
lapply(terms, function(x){
  gseaNb(object = ego,
         geneSetID = x,
         newGsea = T,
         addPval = T,
         pvalX = 0.75,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0,
         subPlot = 2)
}) -> gseaList1

# combine
cowplot::plot_grid(plotlist = gseaList1,ncol = 2,align = 'hv')
dev.off()
######一张图多个通路 最多三个三个
pdf("一张图多个通路.pdf",width = 6, height = 4)

gseaNb(object = ego,
       geneSetID = c('GO:0006954',
                     'GO:0006955','GO:0008219'),
       subPlot = 2,
       termWidth = 35,
       legend.position = c(0.8,0.8))
dev.off()



##########FBMSC
Mono <- RunDEtest(srt = Mono, 
                  group_by = "Group",
                  group1="Old",group2="FBMSC",
                  fc.threshold =1, only.pos = FALSE)
VolcanoPlot(Mono, group_by = "custom")
ggsave("4_FBMSC火山图.pdf", height = 4, width = 4)


Q <- Mono@tools$DEtest_custom$AllMarkers_wilcox

library(clusterProfiler)
library(org.Hs.eg.db)
library(GseaVis)
# load diff
Q$avg_log2FC
diff <- Q %>%
  arrange(desc(avg_log2FC))

genelist <- diff$avg_log2FC
names(genelist) <- diff$gene

# check
head(genelist,3)
#    Aplnr    Foxf1     Bmp5
# 13.45176 13.35322 12.02845


#富集分析
ego <- gseGO(geneList     = genelist,
             OrgDb        = org.Hs.eg.db,
             ont          = "BP",
             keyType      = 'SYMBOL',
             minGSSize    = 5,
             maxGSSize    = 500,
             pvalueCutoff = 1,
             verbose      = FALSE)
b <- ego@result
write.table(b,"FBMSCBP.txt",sep = "\t",quote = F,row.names = F)


#########多个GSEA

ego
terms <- c('GO:0006955',
           'GO:0006954',
           'GO:0045333',
           'GO:0008219',"GO:0012501","GO:0000723",
           "GO:0006281","GO:0051276")
pdf("top6单独NES比例图.pdf",width = 8, height = 10)
lapply(terms, function(x){
  gseaNb(object = ego,
         geneSetID = x,
         subPlot=2,
         addPval = T,
         pvalX = 0.5,pvalY = 0.7,
         pCol = 'black',
         pHjust = 0)
}) -> gseaList

# combine
cowplot::plot_grid(plotlist = gseaList,ncol = 2,align = 'hv')
dev.off()

pdf("top6单独NES比例图新样式.pdf",width = 8, height = 10)

######多样式
lapply(terms, function(x){
  gseaNb(object = ego,
         geneSetID = x,
         newGsea = T,
         addPval = T,
         pvalX = 0.75,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0,
         subPlot = 2)
}) -> gseaList1

# combine
cowplot::plot_grid(plotlist = gseaList1,ncol = 2,align = 'hv')
dev.off()
######一张图多个通路 最多三个三个
pdf("4_FBMSC一张图多个通路.pdf",width = 5, height = 4)

gseaNb(object = ego,
       geneSetID = c('GO:0006979',
                     'GO:0006954','GO:0030330'),
       subPlot = 2,
       termWidth = 35,
       curveCol=color_cluster,
       rankCol=c("white","#B6D4C6","#79AF97"),
       addPval=T,
       topGeneN=5,
       legend.position = c(0.8,0.8))
dev.off()
