#library
library(limma)
library(pheatmap)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)

mydata<-read.table("brain.txt",header = T,sep="\t",row.names = 1)
mydata<-as.matrix(mydata)

colnames(mydata)

# Assuming 'geneExp' is your matrix of gene expression data
# and 'group' is your factor variable as defined earlier
# Create a factor for the groups
group <- factor(rep(c("ABPC_EVs",  "ABMC_EVs", "FBMC_EVs", "Old"), each = 3))

# Define the comparison groups
comparison_groups <- c("ABPC_EVs", "FBMC_EVs", "ABMC_EVs")

# Load the limma package
library(limma)

# Loop through each comparison group
for (group_name in comparison_groups) {
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  # Fit the linear model
  fit <- lmFit(mydata, design)
  
  # Create contrasts: current group vs Young
  contrast_formula <- paste(group_name, "-Old", sep = "")
  contrast_name <- paste(group_name, "vs Old", sep = "_")
  contrasts <- makeContrasts(contrasts = setNames(contrast_formula, contrast_name), levels = design)
  
  # Apply the contrasts
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)
  
  # The correct coefficient name is the one used in contrast formula
  coef_name <- paste(group_name, "-Old", sep = "")
  
  # Get the top table with the correct coefficient
  allgene <- topTable(fit2, coef = coef_name, adjust.method = "BH", number = Inf)
  # Ensure row names (gene names) are preserved
  allgene <- cbind(Gene = rownames(allgene), allgene)
  # Create a directory for the group if it doesn't exist
  dir_name <- paste("Results", group_name, sep = "_")
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
  
  # Filter for differentially expressed genes
  diffgene <- allgene[abs(allgene$logFC) >= 1 & allgene$P.Value < 0.05, ]
  upgene <- allgene[allgene$logFC >= 1 & allgene$P.Value < 0.05, ]
  downgene <- allgene[allgene$logFC <= -1 & allgene$P.Value < 0.05, ]
  
  # Write the results to files in the respective directory
  write.table(allgene, file = paste0(dir_name, "/allgene_", group_name, "_vs_Old.txt"), sep = "\t", quote = FALSE, row.names = F)
  write.table(diffgene, file = paste0(dir_name, "/diffgene_", group_name, "_vs_Old.txt"), sep = "\t", quote = FALSE, row.names = F)
  write.table(upgene, file = paste0(dir_name, "/upgene_", group_name, "_vs_Old.txt"), sep = "\t", quote = FALSE, row.names = F)
  write.table(downgene, file = paste0(dir_name, "/downgene_", group_name, "_vs_Old.txt"), sep = "\t", quote = FALSE, row.names = F)
}

#GO

library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

pvalueFilter=0.05       
qvalueFilter=0.05       


colorSel="pvalue"
ontology.col=c("#CD9DC5", "#A6BAB1", "#6CA8AF")

##########################UP####################
#将基因名转换为基因id   
diffgene <- read.table("up.txt",header = T,sep = "\t")
genes=unique(as.vector(diffgene[,1]))
entrezIDs=mget(genes, org.Rn.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集
kk=enrichGO(gene=gene, OrgDb=org.Rn.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter),]
#保存结果
write.table(GO, file="upGO.txt", sep="\t", quote=F, row.names = F)

#定义GO的数目
library(dplyr)
library(ggplot2)
library(RColorBrewer)

data0=read.delim("clipboard")
colnames(data0) #包含信息列查看
#指定绘图顺序（转换为因子）：
data0$Description <- factor(data0$Description,levels = rev(data0$Description))


colnames(data0)
######
#1.常规画法：
#自定义主题：
mytheme <- theme(axis.title = element_text(size = 16),
                 axis.text = element_text(size = 14),
                 plot.title = element_text(size = 17,
                                           hjust = 0.5,
                                           face = "bold"),
                 legend.title = element_text(size = 16),
                 legend.text = element_text(size = 18))

#2.论文画法：
pdf(file="1_upGO.pdf", width=4, height=8)
mytheme2 <- mytheme + theme(axis.text.y = element_blank()) #先在自定义主题中隐去y轴文本标签显示
data0$text_x <- rep(0.03,24) #新增一列重复数组，使绘图时文本标签能从固定位置开始；
p2 <- ggplot(data = data0,
             aes(x = pvalue, y = Description)) +
  geom_bar(aes(fill = pvalue), stat = "identity", width = 0.9, alpha = 1) +
   scale_fill_gradient(high = "#79AF97",low = "#B6D4C6")+
  labs(x = "-Log10 P value", y = "pathway", title = "Genes specifically elevated in ABPC-EVs") +
  geom_text(aes(x = text_x, #用新增的重复数组控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust=0，左对齐
  theme_classic() + mytheme2
p2
dev.off()
#79AF97
############################down#############################
#定义GO的数目
library(dplyr)
library(ggplot2)
library(RColorBrewer)

data0=read.delim("clipboard")
colnames(data0) #包含信息列查看
#指定绘图顺序（转换为因子）：
data0$pathway <- factor(data0$Description,levels = rev(data0$Description))


colnames(data0)
######
#1.常规画法：
#自定义主题：
mytheme <- theme(axis.title = element_text(size = 16),
                 axis.text = element_text(size = 14),
                 plot.title = element_text(size = 17,
                                           hjust = 0.5,
                                           face = "bold"),
                 legend.title = element_text(size = 16),
                 legend.text = element_text(size = 14))
#Top20富集数目条形图：
p <- ggplot(data = data0,
            aes(x = Gene, y = Description, fill = -LogP))+
  scale_fill_distiller(palette = "RdPu",direction = 1) +
  geom_bar(stat = "identity", width = 0.8) +
  labs(x = "Number of Gene", y = "pathway", title = "GO enrichment barplot") +
  theme_bw() + mytheme
p


#2.论文画法：
pdf(file="1_downGO.pdf", width=6, height=6)
mytheme2 <- mytheme + theme(axis.text.y = element_blank()) #先在自定义主题中隐去y轴文本标签显示
data0$text_x <- rep(0.03,9) #新增一列重复数组，使绘图时文本标签能从固定位置开始；
p2 <- ggplot(data = data0,
             aes(x = Enrichment, y = Description)) +
  geom_bar(aes(fill = Gene), stat = "identity", width = 0.4, alpha = 1) +
  scale_fill_gradient(high = "#1C84BC",low = "#C2E3F6")+
  labs(x = "Enrichment score", y = "pathway", title = "Genes specifically elevated in ABPC-EVs") +
  geom_text(aes(x = text_x, #用新增的重复数组控制文本标签起始位置
                label = pathway),
            hjust = 0)+ #hjust=0，左对齐
  theme_classic() + mytheme2
p2
dev.off()
