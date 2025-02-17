####常识####
B<- read.csv("CAF.xlsx",row.names = 1,header = TRUE)
A<- read.table("HMY表达矩阵TCGA+GEO.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE,row.names = 1)
write.table(A, file = "GEO.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
A <- `上调基因6字名字`#带中文的赋值
A=read_dta("Work_Retirement_and_Pension.dta")
write.table(A, file = "Work_Retirement_and_Pension.xlsx",sep = "\t",row.names = T,col.names = T,quote = T)
B=colnames(A)
BiocManager::install("Rtools")



adev.off()#去图
rm(list = ls())#去变量
A <- A[,-ncol(A)]   #去除最后一列
A<- exp[,(1:6)]   #去除第2列
exp<-exp[,-(103:306)]   #去除第41-72行
A <- A[!duplicated(A$Gene),]#去重复
exp<-na.omit(expr)#删除包含空值的行
rownames(rt)=rt[,1]#行名用第一行
deg2 <- deg[deg$change != "stable",]#提取非stable数据
splots = list()#新增一个空表

A$group=data1$Y

rt <- as.data.frame(rt)#变表格
exp <- exp %>% t() %>% as.data.frame()#倒置，变表格
A <- t(data)#颠倒横纵
rt <- as.matrix(rowdata$gene_name)#变非表格
dat$label <- as.factor(dat$label)#变factor
dat$label <- as.numeric(dat$label)#变numeric
class(gset)  #查看数据类型


A$Gene <- as.character(rownames(A))#新增Gene Symbol
A$Gene <- substr(A$Gene,1,12)#保留列名前15位
A <- A[!duplicated(A$Gene),]#去重复
rownames(A) <- A$Gene #将行名变为Gene Symbol


countsintegers=round(counts_01A, digits = 0)#取整
DEG <- DEG[order(DEG$P.Value), ]#对差异基因的p值进行从小到大的排序
DEG <- DEG[DEG$P.Value<0.05, ]#仅保留P<0.05的

#包的安装与升级
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #中科大镜像
packageVersion("Rtools")#查看版本
remove.packages("promises")#删除
update.packages("promises")#升级
remotes::install_github('promises')#删了重下
install.packages("fromto")#普通安装
BiocManager::install("promises")#下载
install.packages("~/celldex_1.0.0.tar.gz",repos = NULL,type ="source")#特殊安装
.libPaths()#安装包路径
find.package("xCell")#寻找安装位置



library(fromto)
install.packages("installr")#更新R
require(installr)
updateR()
traceback()#寻找错误在哪

#提取基因子集
#总数据，#（11920*371）
exp<-read.table("TCGA肝癌mRNA_counts表达矩阵.txt", row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F,header = T)

exp<-read.table("TCGA肝癌表达矩阵counts01A(倒置).txt", row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F,header = T)
#子集
A<- read.table("交集基因.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
genes<- colnames(A)[2:11]#批量计算那些列
genes=c("OS","time","FKBP9","SEPT6")
genes=c("ITGAX","NFAM1","ACTA2","CCR7","CD8A","METTL3","IGFBP3","SPRED1")
genes=c("FKBP1A","FKBP1B","PRDX5","FKBP3","FKBP4","FKBP5","NINL","FKBP7",
        "FKBP8","FKBP9","FKBP10","FKBP11","FKBP14","FKBP15")
genes=c("FKBP1A","FKBP1B","FKBP2","FKBP3","FKBP4","FKBP5","FKBP6","FKBP7",
        "FKBP8","FKBP9","FKBP10","FKBP11","FKBP14","FKBP15")
genes=c("HDAC2","BRCA1","RCOR1","SNAI1","XRCC6","PHF21A","KDM1A","MECP2",
        "HDAC1","TRIM28")
									
exp1 <- A[genes,]
write.table(exp1, file = "top10.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


#取交集
A=t(A)
B=t(exp)
comgene <- intersect(rownames(B),rownames(A))
A <- A[comgene,]
B <- B[comgene,]
AB <- cbind(B,A)
write.table(B, file = "TCIA基因.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)




####——TCGA——####
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) #对应清华源
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #对应中科大源
#if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
#BiocManager::install("TCGAbiolinks") 
library(TCGAbiolinks)
library(SummarizedExperiment)
##数据处理与下载
getProjectSummary('TCGA-UCEC')
TCGAbiolinks:::getProjectSummary('TCGA-UCEC')
#下载数据，最好登陆TCGA核对一下信息
query<-GDCquery('TCGA-UCEC',data.category = 'Transcriptome Profiling',
                data.type = 'Gene Expression Quantification',
                workflow.type = 'STAR - Counts')
#如果是常见的RNAseq，STAR - Counts要写成HTSeq - Counts



##1表达矩阵####
GDCdownload(query, files.per.chunk = 50) #每次下载50个文件
#GDCprepare(query,save = T,save.filename = "胶质瘤.rdata")#保存一下
#load("肝癌.rdata")
exps<-GDCprepare(query = query)
exp.matrix<-assay(exps)
exp<-na.omit(exp)#删除包含空值的行
#write.table(exp, file = "TCGA胰腺癌原始表达矩阵.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


#####mRNA####
load("胰腺癌.rdata")
se <- data#exps
se
names(assays(se))
rowdata <- rowData(se)# 提取rowData
names(rowdata)#看看rowData包括哪些内容，可以看到里面有我们需要的gene_name和gene_type
table(rowdata$gene_type)# gene_type是基因类型，看看有哪些
head(rowdata$gene_name)#gene_name就是gene_symbol，就是我们需要的
length(rowdata$gene_name)
rowdata <- rowData(se)
#分别提取mRNA和lncRNA的SummarizedExperiment对象
#mRNA
se_mrna <- se[rowdata$gene_type == "protein_coding",]
se_mrna#还是一个SummarizedExperiment对象，神奇！
# mRNA的counts矩阵
expr_counts_mrna <- assay(se_mrna,"unstranded")
# mRNA的tpm矩阵
expr_tpm_mrna <- assay(se_mrna,"tpm_unstrand")
# mRNA的fpkm矩阵
expr_fpkm_mrna <- assay(se_mrna,"fpkm_unstrand")
# 先提取gene_name
symbol_mrna <- rowData(se_mrna)$gene_name
head(symbol_mrna)
#合并表达矩阵，需要那个就合并那个expr_counts_mrna
expr_counts_mrna_symbol <- cbind(data.frame(symbol_mrna),
                                 as.data.frame(expr_counts_mrna))
expr_fokm_mrna_symbol <- cbind(data.frame(symbol_mrna),
                                 as.data.frame(expr_fpkm_mrna))
#去重复也很简单，这里我们保留最大的那个。
suppressPackageStartupMessages(library(tidyverse))
expr_read1 <- expr_counts_mrna_symbol %>% #这里更换再保存
  as_tibble() %>% 
  mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
  arrange(desc(meanrow)) %>% 
  distinct(symbol_mrna,.keep_all=T) %>% 
  select(-meanrow) %>% 
  column_to_rownames(var = "symbol_mrna") %>% 
  as.data.frame()
expr_read1.1 <- expr_fokm_mrna_symbol %>% #这里更换再保存
  as_tibble() %>% 
  mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
  arrange(desc(meanrow)) %>% 
  distinct(symbol_mrna,.keep_all=T) %>% 
  select(-meanrow) %>% 
  column_to_rownames(var = "symbol_mrna") %>% 
  as.data.frame()
#标准化
ex <- expr_read1
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0,arr.ind = T)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
exp1 <- exprSet#如果不需要转换就显示没exprSet，是对的。
#标准化
ex1.1 <- expr_read1.1
qx <- as.numeric(quantile(ex1.1, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex1.1[which(ex1.1 <= 0,arr.ind = T)] <- NaN
exprSet1.1 <- log2(ex1.1)
print("log2 transform finished")}else{print("log2 transform not needed")}
exp1.1 <- exprSet1.1#如果不需要转换就显示没exprSet，是对的。
#删除包含空值的行
exp1<-na.omit(exp1)
exp1.1<-na.omit(exp1.1)

write.table(exp1, file = "TCGA肝癌mRNA_counts表达矩阵.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(exp1.1, file = "TCGA肝癌mRNA_fpkm表达矩阵.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#####lncRNA####
se_lnc <- se[rowdata$gene_type == "lncRNA"]
se_lnc
# lncRNA的counts矩阵
expr_counts_lnc <- assay(se_lnc,"unstranded")
# lncRNA的tpm矩阵
expr_tpm_lnc <- assay(se_lnc,"tpm_unstrand")
# lncRNA的fpkm矩阵
expr_fpkm_lnc <- assay(se_lnc,"fpkm_unstrand")
symbol_lnc <- rowData(se_lnc)$gene_name
head(symbol_lnc)
expr_counts_lncrna_symbol <- cbind(data.frame(symbol_lnc),
                                 as.data.frame(expr_counts_lnc))
expr_read2 <- expr_counts_lncrna_symbol %>% #这里更换再保存
  as_tibble() %>% 
  mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
  arrange(desc(meanrow)) %>% 
  distinct(symbol_lnc,.keep_all=T) %>% 
  select(-meanrow) %>% 
  column_to_rownames(var = "symbol_lnc") %>% 
  as.data.frame()
expr_read2<-na.omit(expr_read2)
#标准化
ex <- expr_read2
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0,arr.ind = T)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
exp2 <- exprSet#如果不需要转换就显示没exprSet，是对的。
write.table(exp2, file = "宫颈癌lncRNA_counts表达矩阵.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#mRNA+lncRNA（后期mRNA不够用再考虑）
exp1<-t(exp1)
exp2<-t(exp2)
exp3<- cbind(exp1,exp2)
exp3<-t(exp3)
write.table(exp3, file = "宫颈癌mRNA+lncRNA表达矩阵.txt",sep = "\t",row.names = T,col.names = NA,quote = F)



##2临床信息####
library(limma)
library(stringr)
#exp1 <- read.table("TCGA肝癌mRNA_counts表达矩阵.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#clinical <- read.table("TCGA肝癌临床信息.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

clinical<-GDCquery_clinic('TCGA-LIHC')
write.table(clinical, file = "TCGA肝癌临床信息.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(clinical, file = "TCGA肝癌临床信息.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
cl<-clinical
col <- data.frame(colnames(cl))
cl2 <-cl
dim(cl2)#临床信息377
dim(exp1)#表达矩阵424
#癌/癌旁分组：
group_list <- ifelse(as.numeric(str_sub(colnames(exp1),14,15)) < 10,'tumor','normal')
table(group_list)
#筛选表达矩阵中的肿瘤样本：
exp2 <- exp1[,group_list == c("tumor")]
dim(exp2)#剩374个tumor
colnames(exp2) <- substr(colnames(exp2),1,12)#保留列名前12位
#将cl样本行名的顺序按exp2列名进行匹配：
cl3 <- cl2[match(colnames(exp2),cl2$submitter_id),]
dim(cl3) #匹配完306
#检查顺序是否一致：
identical(colnames(exp2),cl3$submitter_id)
#如果匹配那就合并
exp3 <- t(exp2)#横纵颠倒
C <- cbind(cl3,exp3)
C <- C[!duplicated(C$submitter_id),]#去重复
rownames(C) <- C$submitter_id

write.table(cl3, file = "TCGA肝癌临床参数01A.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(exp3, file = "TCGA肝癌表达矩阵counts01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(exp2, file = "TCGA肝癌表达矩阵counts01A(倒置).txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#write.table(C, file = "TCGA肝癌表达矩阵+临床参数(总).txt",sep = "\t",row.names = T,col.names = NA,quote = F)



####3保留有影像的数据####
A<- read.table("TCGA肝癌表达矩阵counts01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
B <- read.table("TCGA肝癌临床信息处理.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

A <- A[!duplicated(A[,1]),]#去重复
rownames(A) <- A[,1]
A=A[,-1]
comgene <- intersect(rownames(A),rownames(B))
B <- B[comgene,]
write.table(B, file = "表达矩阵(倒置).txt",sep = "\t",row.names = T,col.names = NA,quote = F)
C <- t(C)#颠倒横纵
B <- as.data.frame(B)#变表格
B<-na.omit(B)
write.table(B, file = "TCGA表达矩阵使用.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
C <- C[!duplicated(C$submitter_id),]#去重复
rownames(C) <- C$submitter_id
comgene2 <- intersect(rownames(A),rownames(C))
C <- C[comgene2,]
write.table(C, file = "TCGA临床信息使用.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
#这里的临床信息还是需要自己手动修改的，修完之后可以和表达矩阵合并用于分析



####4差异分析####
library(limma)
library(stringr)
#加入影像组学得分分组
A<- read.table("临床信息.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
table(A$label)
group_list <- ifelse(str_detect(A$label, "1"), "tumor","normal")
group_list = factor(group_list,levels = c("normal","tumor"))
table(group_list)
#对表达矩阵进行差异分析
exp <- read.table("TCGA肝癌表达矩阵counts01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
##标记上下调基因
logFC=1
P.Value = 0.05
k1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))#exp报错subscript out of bounds
table(deg$change)

write.table(deg, file = "上调基因.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
#这里已经选出差异基因，如果数目少则直接免疫浸润、富集分析等
#如果数目多，强烈建议与别的数据集取交集，再分析


#热图###
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(ggthemes)
cg = rownames(deg)[deg$change !="stable"]
exp_diff <- exp[cg,]
group_list=factor(ifelse(substr(colnames(exp),14,16) == "01A","T","N"),levels = c("N","T"))
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(exp_diff)
pheatmap(exp_diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)

#火山图###
DEG<-deg
DEG$logP <- -log10(DEG$P.Value)
DEG$Label = ""   #新加一列label
DEG <- DEG[order(DEG$P.Value), ]   #对差异基因的p值进行从小到大的排序
DEG$Gene <- rownames(DEG)
#高表达的基因中，选择fdr值最小的5个
up.genes <- head(DEG$Gene[which(DEG$change == "up")], 5)
#低表达的基因中，选择fdr值最小的5个
down.genes <- head(DEG$Gene[which(DEG$change == "down")], 5)
#将up.genes和down.genes合并，并加入到Label中
DEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
DEG$Label[match(DEG.top5.genes, DEG$Gene)] <- DEG.top5.genes

ggscatter(DEG, x = "logFC", y = "logP",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1,
          label = DEG$Label,
          font.label = 8,
          repel = T,
          xlab = "log2FoldChange",
          ylab = "-log10 (Adjust P-value)") +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")








#——GEO——####
library(tidyverse)
library(stringr)
library(GEOquery)
library(limma) 
#BiocManager::install("limma")
####1下载数据####
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) #对应清华源
gset = getGEO('GSE193928', destdir=".", AnnotGPL = F, getGPL = F)
class(gset)
###提取子集
gset[[1]]
# 通过exprs函数表达矩阵
exp <- exprs(gset[[1]])
#通过pData函数获取分组信息
pdata <- pData(gset[[1]])
write.table(pdata, file = "GSE193928临床信息.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
#下载好之后自己分类或者利用代码分类
table(pdata$title)
group_list <- ifelse(str_detect(pdata$title, "with meta"),"tumor", "normal")#NCR=1=tumor
group_list = factor(group_list,levels = c("normal","tumor"))
table(group_list)

#对得到的表达矩阵操作
ex <- exp
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
exp <- exprSet#如果不需要转换就显示没exprSet，是对的。
#删除包含空值的行
exp<-na.omit(exp)
write.table(exp, file = "GSE205209表达矩阵标准化.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


#2差异分析####
library(limma)
library(stringr)
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
##标记上下调基因
logFC=1
P.Value = 0.05
k1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
write.table(deg, file = "GSE158408差异基因.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
#手动挑选上调基因做成上调基因txt



####3ID转换####
#官网下载GPL探针，然后与上调基因取交集
A<- read.table("GSE158408差异基因.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
B<- read.table("GPL23159.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
A <- deg
B <- GPL570
comgene <- intersect(rownames(A),rownames(B))
A <- A[comgene,]
B <- B[comgene,]
C <- cbind(A,B)
#colnames(C)
#C <- C[!duplicated(C$Symbol),]   #去重复
#rownames(C) <- C$Symbol   #将行名变为Gene Symbol
write.table(C, file = "GSE158408差异基因带名字.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)



library(ggthemes)
#火山图###
DEG<-deg#C
DEG$logP <- -log10(DEG$P.Value)
DEG$Label = ""   #新加一列label
DEG <- DEG[order(DEG$P.Value), ]   #对差异基因的p值进行从小到大的排序
DEG$Gene <- rownames(DEG)
#高表达的基因中，选择fdr值最小的5个
up.genes <- head(DEG$Gene[which(DEG$change == "up")], 5)
#低表达的基因中，选择fdr值最小的5个
down.genes <- head(DEG$Gene[which(DEG$change == "down")], 5)
#将up.genes和down.genes合并，并加入到Label中
DEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
DEG$Label[match(DEG.top5.genes, DEG$Gene)] <- DEG.top5.genes

ggscatter(DEG, x = "logFC", y = "logP",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1,
          label = DEG$Label,
          font.label = 8,
          repel = T,
          xlab = "log2FoldChange",
          ylab = "-log10 (Adjust P value)") +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")




####4取交集####
A<- read.table("TCGA表达矩阵.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
B<- read.table("肝癌临床信息01A.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#A <- `上调基因`
#B <- `放疗相关基因`
A=data
B=data3
rownames(data)=substr(rownames(data),1,12)#保留列名前12位
data3$X=substr(rownames(data3),1,12)#保留列名前12位
data3$Y<- substr(rownames(data3),14,28)#保留列名前12位

data$Y <- substr(data$T,14,28)
data3 <- data3[!duplicated(data3$submitter_id),]#去重复
rownames(data3)=data3$submitter_id
B=rownames(data3)
B=data3


comgene <- intersect(rownames(B),rownames(A))
B <- B[comgene,]
A <- A[comgene,]
AB <- cbind(B,A)
write.table(A, file = "TCGA表达矩阵.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)




#—TCGA&GEO—####
#TCGA+GEO 
library(limma)#先自行处理好两家的表达矩阵
C<- read.table("GSE76427表达矩阵带名字.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp1<- read.table("TCGA表达矩阵使用.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
colnames(C)#把多余的列删掉
C<- C[,-(168:196)] #去除第41-72行
sameGene=intersect(rownames(C),rownames(exp1))
TCGA=exp1[sameGene,,drop=F]
GEO=C[sameGene,,drop=F]
##TCGA和GEO组间矫正
TCGA=normalizeBetweenArrays(TCGA)
GEO=normalizeBetweenArrays(GEO)
###组间校正后数据合并
data0=cbind(TCGA,GEO)

#GEO+GEO



#ComBat去批次效应####
BiocManager::install("factoextra")
library(sva)
library(FactoMineR)
library(factoextra)
TCGA_smaple=as.data.frame(colnames(TCGA))
colnames(TCGA_smaple)="group"
TCGA_smaple$batch="TCGA"
GEO_sample=as.data.frame(colnames(GEO))
colnames(GEO_sample)="group"
GEO_sample$batch="GEO"
sample=rbind(TCGA_smaple,GEO_sample)

mod = model.matrix(~1,data = sample)
batch=sample$batch
data0_combat=ComBat(dat=data0,batch = batch,mod=mod,par.prior = T)

write.table(data0_combat, file = "TCGA+GEO表达矩阵.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(data0_combat, file = "TCGA+GEO表达矩阵.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


#绘制去批次前箱线图及PCA图
#boxplot(log2(data0+1),outline=T, notch=T,las=2)#样本数量太多时间长
##PCA图
batch=sample$batch
dat0=as.data.frame(t(data0))
dat.pca0 <- PCA(dat0, graph = FALSE)
pca_plot0 <- fviz_pca_ind(dat.pca0,geom.ind = "point",
                          col.ind = batch,
                          palette = c("#00AFBB", "#E7B800"),
                          addEllipses = TRUE, 
                          legend.title = "Groups")
pca_plot0
#ggsave(plot = pca_plot0,filename ="PCA_FPKM_to_GEO.pdf")

#绘制去批次后箱线图及PCA图
#boxplot(log2(data0_combat+1),outline=T, notch=T,las=2)#箱线图
##绘制PCA图
dat1=as.data.frame(t(data0_combat))
dat.pca0_combat <- PCA(dat1, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca0_combat,geom.ind = "point",
                         col.ind = batch,
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, 
                         legend.title = "Groups")
pca_plot
#ggsave(plot = pca_plot,filename ="PCA_FPKM_to_GEO_ComBat.pdf")


#差异分析####
library(stringr)
#手动为GEO表头添加-01A以区分癌与癌旁
exp <- read.table("TCGA+GEO表达矩阵.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group_list <- ifelse(str_detect(colnames(exp), "01A"), "tumor","normal")
group_list = factor(group_list,levels = c("normal","tumor"))
table(group_list)

design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
##标记上下调基因
logFC=1
P.Value = 0.05
k1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))#exp报错subscript out of bounds
table(deg$change)
write.table(deg, file = "上调基因.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(deg, file = "上调基因.txt",sep = "\t",row.names = T,col.names = NA,quote = F)



####铁死亡/铜死亡####
A<- read.table("铁死亡数据集.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
A<- read.table("铜死亡数据集.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
A<- read.table("GeneCardsm6A数据集.csv",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
A<- read.table("GeneCards代谢相关数据库.csv",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
B<- read.table("上调基因.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

comgene <- intersect(rownames(A),rownames(B))
A <- A[comgene,]
B <- B[comgene,]
AB <- cbind(A,B)
write.table(AB, file = "铜死亡交集基因.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
#下一步sting网站找hub基因，然后进行基因筛选



####批量COX####
rm(list=ls())
R=read.table("生存分析基因+临床.xlsx",header=T,sep="\t",row.names = 1)
colnames(R)
can <- colnames(R)[3:9]#批量计算那些列
can
uni_sur <- sapply(can, function(x) as.formula(paste('Surv(time, status)~', x)))
uni_cox <- lapply(uni_sur, function(x){coxph(x, data =  R)})
uni_results <- lapply(uni_cox, function(x) {#x <- uni_cox$Gender从每个模型的摘要中提取相关信息
  x <- summary(x)# 获取Wald检验的p值并截断为两位小数
  p.value <- signif(x$wald["pvalue"], digits = 2)# 获取模型系数的估计值（通常表示风险比HR）并截断为两位小数
  HR <- signif(x$coef[2], digits = 2)# 获取HR的95%置信区间的下界和上界并截断为两位小数
  HR.confint.lower <- signif(x$conf.int[,"lower .95"], digits = 2)
  HR.confint.upper <- signif(x$conf.int[,"upper .95"], digits = 2)# 构建HR和其95%置信区间的字符串表示
  HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")# 创建一个包含p值和HR（以及其95%置信区间）的结果向量
  res <- c(p.value, HR)# 设置结果向量的名称，便于后续引用
  names(res) <- c("p.value", "HR (95% CI)")# 返回结果向量
  return(res)
})
result <- as.data.frame(t(as.data.frame(uni_results, check.names = FALSE)))
result
#保存一下模型的原始数据，用于后续的分析
write.table(result, file = "批量单因素COX.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)



####批量K-M####
library(survminer)
library(survival)
######最佳截断值####
data<- read.table("胰岛分析.xlsx",header=T,sep="\t",row.names = 1)
colnames(data)
var =colnames(data[,4:7])#提取前11个列名
res.cut <- surv_cutpoint(data,time="time",event ="status",variables =var)#surv_cutpoint求取
summary(res.cut)
cut <- surv_categorize(res.cut)

plot(res.cut, "PPY")#截断值图像
#write.table(cut, file = "截断值.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)

#单个KM曲线
fit <- survfit(Surv(time, status)~PPY, data = cut)#KM曲线
ggsurvplot(fit, data = cut) #得到基础的生存曲线啦
ggsurvplot(fit, pval = TRUE, palette = "jco", 
           data = cut, legend = c(0.8, 0.8),
           ylab="OS (%)",xlab = "Months", #更改横纵坐标
           conf.int = TRUE, #给生存曲线添加上置信区间
           legend.labs = c("High","Low"), #在图中添加图例
           title ="CD68",risk.table = F)

#KM拼在一起
data=cut
splots = list()
genes<-var 
cl4=data[,(1:2)]#提取前2列作为临床信息
exp3=subset(data,select=genes)#提取子集作为表达矩阵
exprSet=t(exp3)
for(i in 1:length(genes)){
  g = genes[i]
  cl4$genes = ifelse(exprSet[g,]== "high",'high','low')#千万别搞错
  sfit1 = survfit(Surv(time,status) ~genes, data = cl4)
  splots[[i]] =  ggsurvplot(sfit1, pval = TRUE, palette = "jco", 
                            data = cl4, legend = c(0.8, 0.8),
                            ylab="OS (%)",xlab = "Months", #更改横纵坐标
                            conf.int = TRUE, #给生存曲线添加上置信区间
                            legend.labs = c("High","Low"), #在图中添加图例
                            title =genes[i],risk.table = F)}

arrange_ggsurvplots(splots,nrow = 2, ncol = 2)#nrow一列几个，ncol一行几个


#####中位数截断####
data=read.table("胰岛分析.xlsx",header=T,sep="\t",row.names = 1)  
cl4=data[,1:2]
exp3=data[,3:6]
genes =colnames(data[,3:6])#提取前11个列名

exp2=t(exp3)
exp3=as.data.frame(exp2)
exp3<- apply(exp2[1:4,],2,as.numeric) # 将2-10列的数据变成数值型
median(exp3[3,]) # 除此之外还可以是mean
cl4$genes <- ifelse(exp3[1,] > median(exp3[1,]),"High","Low")

#拼在一起
splots = list()
exprSet=t(exp3)
for(i in 1:length(genes)){
  g = genes[i]
  cl4$genes = ifelse(exp2[i,] > median(exp2[i,]),"High","Low")
  sfit1 = survfit(Surv(time,status) ~genes, data = cl4)
  splots[[i]] =  ggsurvplot(sfit1, pval = TRUE, palette = "jco", 
                            data = cl4, legend = c(0.8, 0.8),
                            ylab="OS (%)",xlab = "Months", #更改横纵坐标
                            conf.int = TRUE, #给生存曲线添加上置信区间
                            legend.labs = c("High","Low"), #在图中添加图例
                            title =genes[i],risk.table = F)}

arrange_ggsurvplots(splots,nrow = 2, ncol = 2)#一行几个，一列几个

#####X-tlie截断值####
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE) 
cl4=data[,1:2]
exp3=data[,3:6]
genes =colnames(data[,3:6])#提取前11个列名

exp2=t(exp3)
exp3=as.data.frame(exp2)
exp3<- apply(exp2[1:4,],2,as.numeric) # 将2-10列的数据变成数值型
median(exp3[3,]) # 除此之外还可以是mean
cl4$genes <- ifelse(exp3[1,] > 0,"High","Low")

#拼在一起
splots = list()
exprSet=t(exp3)
for(i in 1:length(genes)){
  g = genes[i]
  cl4$genes = ifelse(exp2[i,] > 0,"High","Low")
  sfit1 = survfit(Surv(time,status) ~genes, data = cl4)
  splots[[i]] =  ggsurvplot(sfit1, pval = TRUE, palette = "jco", 
                            data = cl4, legend = c(0.8, 0.8),
                            ylab="DFS (%)",xlab = "Months", #更改横纵坐标
                            conf.int = TRUE, #给生存曲线添加上置信区间
                            legend.labs = c("High","Low"), #在图中添加图例
                            title =genes[i],risk.table = F)}

arrange_ggsurvplots(splots,nrow = 2, ncol = 2)#一行几个，一列几个




####LASSO-COX回归模型的构建####
library("glmnet")
library("survival")
#####提取基因子集####
exp1 <- read.table("35例临床信息.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#总数据，#（11920*371）
exp<-read.table("CPTAC蛋白组.xlsx", row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F,header = T)
exp=t(exp)
comgene <- intersect(rownames(exp),rownames(exp1))
exp <- exp[comgene,]
write.table(exp, file = "转录组学生存分析.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)

#手动编辑生存信息,status+time
A<- read.table("肝癌临床信息01A.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
B=t(exp)
comgene <- intersect(rownames(A),rownames(B))
B <- B[comgene,]
A <- A[comgene,]
AB <- cbind(A,B)
write.table(AB, file = "蛋白组学生存分析.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)



#####开始分析####
rt <- read.table("HMY表达矩阵TCGA+GEO.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
rt <- as.data.frame(data)#变表格
rt <- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
set.seed(133)
#colnames(rt)
rt$time=rt$time+0.01
x<- as.matrix(rt[,c(5:ncol(rt))])#是否是3要看前面有几列无用数据
y=data.matrix(Surv(rt$time,rt$status))
fit=glmnet(x, y, family = "cox", maxit = 1000,alpha = 1)
plot(fit, xvar = "lambda", label = TRUE,lwd=1.4)
cvfit=cv.glmnet(x,y,family="cox",maxit=1000,type.measure="deviance",nfolds=5)
plot(cvfit)
lambda.min <-fit$lambda.min
lambda.1se <-fit$lambda.1se
###4. 输出预测模型的相关系数与riskScore
coef=coef(fit, s = cvfit$lambda.min)#lambda.1se
index=which(coef != 0)
actCoef=coef[index]#coef
lassoGene=row.names(coef)[index]#diffvariables
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
geneCoef #查看模型的相关系数
write.table(geneCoef, file = "LASSO系数2.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


###4.2 计算riskScore
FinalGeneExp = rt[,lassoGene]
myFun = function(x){crossprod(as.numeric(x),actCoef)}
riskScore = apply(FinalGeneExp,1,myFun)
outCol = c("time", "status", lassoGene)
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
dat = cbind(rt[,outCol], riskScore=as.vector(riskScore), risk)
write.table(dat, file = "LASSO转录组学2.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)

#4.3特征展示
lasso.result.se<-cbind(lassoGene,actCoef)
lasso.result.se<-as.data.frame(lasso.result.se)
lasso.result.se$actCoef<-as.numeric(lasso.result.se$actCoef)
lasso.result.se
ggplot(aes(x=reorder(lassoGene,actCoef),y=actCoef,fill=lassoGene),data=lasso.result.se)+
  geom_col()+
  coord_flip()+
  theme_bw()+
  labs(x="")+
  ggtitle("LASSO identified variables")+
  scale_fill_brewer(palette = "Set3")+
  theme(legend.position = "")

###5. 绘制散点分布图
library(ggpubr)  
p <- ggboxplot(dat, x = "status", y = "riskScore",
               color = "status", palette = "jco",
               add = "jitter")
p <- p + stat_compare_means()   #  Add p-value
p   #得出预测结果

###6. 判断预测结果的准确性
library(ROCR)   #使用ROCR包绘制预测模型的ROC曲线
library(glmnet)
library(caret)
library(pROC)
library(ggplot2)
library(timeROC)#timeROC
ROC=timeROC(T=dat$time, delta=dat$status,
            marker=dat$riskScore, cause=1,
            weighting='aalen',
            times=c(12,36,50), 
            ROC=TRUE)
ROC
plot(ROC,time=12, col="red", lwd=3, title = "")  
plot(ROC,time=36, col="blue", add=TRUE, lwd=3)
plot(ROC,time=50, col="green", add=TRUE, lwd=3)

legend("bottomright",c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],3)), 
                       paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],3)),
                       paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],3))),
       col=c("red", "blue", "green"),lty=1, lwd=3,bty = "n",cex=1.2) #cex为字体大小




#画风险分布图
#y=生存时间
rt <- dat
color=as.vector(rt$status)
color[color==1]="indianred1"
color[color==0]="lightseagreen"
plot(rt$time, pch=19,xlab="Patients (increasing risk socre)", ylab="Survival time (years)",col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("indianred1","lightseagreen"),cex=1.2)
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
abline(v=lowLength,lty=2)
#y=riskscore
rt <- rt[order(rt[,6]),]#riskscore的那一列
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
plot(line, type="p", pch=20,xlab="Patients (increasing risk socre)", ylab="Risk score",
col=c(rep("lightseagreen",lowLength),rep("indianred1",highLength)) )
abline(h=median(rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High risk", "Low risk"),bty="n",pch=19,col=c("indianred1","lightseagreen"),cex=1.2)
    



#——后基因时代——####
#火山图美化####
library(limma)
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(ggthemes)
library(ggplot2)
library(dplyr)
exp <- read.table("43表达矩阵counts01A(倒置).xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- read.table("宫颈癌表达矩阵counts01A(倒置).txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
B <- read.table("肝癌临床信息01A.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE,row.names = 1)

exp<-na.omit(exp1)#删除包含空值的行
A <- t(exp)#颠倒横纵
colnames(B)
#根据临床资料分组
group_list <- ifelse(str_detect(A$group, "Brucine"), "tumor","normal")
#group_list <- ifelse(str_detect(rownames(A), "Brucine"), "tumor","normal")#大NC的对应第一个tumor，别搞反了
group_list = factor(group_list,levels = c("normal","tumor"))
table(group_list)


#或者根据基因分组
gene="CDH1"
med=median(as.numeric(A[,gene]))#取中位数14.947
group_list = factor(ifelse(A[,gene]>11.94,"high","low"),levels = c("low","high"))#顺序相反，别搞反了


design=model.matrix(~group_list)
fit=lmFit(exp,design)#exp可能需要颠倒，正确的是(13019*306)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
##标记上下调基因
logFC_t1 = 1
P.Value_t1 = 0.05
logFC_t2 = 2
P.Value_t2 = 0.01
dat = deg
dat$symbol=rownames(deg)
k1 = with(dat,logFC > logFC_t2 & P.Value<P.Value_t2);table(k1)
k2 = with(dat,logFC < -logFC_t2 & P.Value<P.Value_t2);table(k2)
k3 = with(dat,logFC > logFC_t1 & P.Value < P.Value_t1 );table(k3)
k4 = with(dat,logFC < -logFC_t1 & P.Value <P.Value_t1 );table(k4)
deg$change = ifelse(k3,"up",ifelse(k4,"down","stable"))
table(deg$change)
# 设置不同颜色和大小
my_color = case_when(k1~"red",k2~"green",k3~"#f7c0be",k4~"#abd486",TRUE~"#4d4d4d")
my_size = case_when(k1|k2~2.5,k3|k4~2,TRUE~1.5)
p = ggplot(dat,aes(x = logFC,y = -log10(P.Value))) +
  geom_point(alpha=0.5, size=my_size, color=my_color) +
  geom_vline(xintercept = c(-logFC_t1,logFC_t1,-logFC_t2,logFC_t2),lty= 4,lwd=0.8,alpha = c(0.5,0.5,1,1)) +
  geom_hline(yintercept = c(-log10(P.Value_t1),-log10(P.Value_t2)),lty= 4,lwd=0.8,alpha = c(0.5,1)) +
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(-7, 7), expand = c(0,0))
p

DEG<-deg
DEG$logFC <- abs(DEG$logFC)#取绝对值
DEG$Label = ""   #新加一列label
DEG <- DEG[order(-DEG$logFC), ]#对差异基因的p值进行从大到小的排序
DEG$Gene <- rownames(DEG)
#高表达的基因中，选择fdr值最小的5个
up.genes <- head(DEG$Gene[which(DEG$change == "up")], 5)
#低表达的基因中，选择fdr值最小的5个
down.genes <- head(DEG$Gene[which(DEG$change == "down")], 5)
#将up.genes和down.genes合并，并加入到Label中
DEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
DEG$Label[match(DEG.top5.genes, DEG$Gene)] <- DEG.top5.genes
g7=DEG.top5.genes
# 加标签
for_label <- dat %>%filter(symbol %in% g7)
p + geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data =for_label,
    color="black")


write.table(deg, file = "差异基因.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)



####富集分析#####
library(tidyverse)
#library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
#基因名字为首列
DEG <- read.table("差异基因.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
deg <- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE,row.names = 1)
DEG <- deg[deg$P.Value<0.05, ]#仅保留P<0.05的,要不然结果出奇的一致


DEG <- DEG %>% rownames_to_column("Gene")
genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))


#####1GO富集####
ego <- enrichGO(gene = DEG$ENTREZID,keyType = "ENTREZID",
                OrgDb = org.Hs.eg.db, ont = "all",
                pAdjustMethod = "BH",minGSSize = 1,
                pvalueCutoff =0.05,qvalueCutoff =0.05,readable = TRUE)
ego_res <- ego@result
##柱状图+##气泡图
barplot(ego, drop = TRUE, showCategory =6,split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale='free')
#dotplot(ego,showCategory = 6,split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale='free')


#####2KEGG富集####
kk <- enrichKEGG(gene = DEG$ENTREZID,organism= 'hsa',
                 pvalueCutoff = 0.05,qvalueCutoff =0.05)
kk_res <- kk@result

barplot(kk, showCategory = 20,color = "pvalue")#柱状图
#dotplot(kk, showCategory = 20)#气泡图


#####3GSEA富集####
library(enrichplot)
msigdb_GMTs <- "msigdb_v7.0_GMTs"
msigdb <- "H.all.v7.0.entrez.gmt"#H.all.v7.0.entrez.gmt（全人类）
#c2.all.v7.0.entrez.gmt或c5.all.v7.0.entrez.gmt
#一般选择C2模块；如果进行GO富集分析选择C5模块如果进行motif分析选择C3模块
#读取上面指定的gmt文件,需要自己下载
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))

#读取差异基因文件
DEG<- read.table("差异基因.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
DEG <- deg[deg$P.Value<0.05, ]#仅保留P<0.05的,要不然结果出奇的一致
#   DEG <-exp

DEG <- DEG %>% rownames_to_column("Gene")
genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))
geneList = DEG$logFC# geneList =DEG$log2.Fold_change.
names(geneList) = as.character(DEG[,'ENTREZID'])
geneList = sort(geneList, decreasing = TRUE)
set.seed(123)
#GSEA<-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
GSEA<- GSEA(geneList,                 # 排序后的gene
       TERM2GENE = kegmt, # 基因集
       pvalueCutoff = 0.05,      # P值阈值
       minGSSize = 20,           # 最小基因数量
       maxGSSize = 1000,         # 最大基因数量
       eps = 0,                  # P值边界
       pAdjustMethod = "BH")     # 校正P值的计算方法


gseaplot2(GSEA,14,color="red",pvalue_table=TRUE)#第一个
gseaplot2(GSEA,geneSetID= c(2:4),subplots = 1:3,pvalue_table =T)#前10个

dotplot(GSEA, showCategory = 10, color = "p.adjust")# 展示富集到的通路，我们这里选择展示前15个
# 将通路分为激活和抑制两个部分进行展示
dotplot(GSEA, showCategory = 10, split = ".sign") + facet_grid(~.sign) +
  theme(plot.title = element_text(size = 10, color = "black", hjust = 0.5),
        axis.title = element_text(size = 10,color = "black"), 
        axis.text = element_text(size = 10,color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 1 ),
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))




#仙桃可视化####
rm(list = ls())
library(GOplot)
df1 <- read.table("仙桃GOKEGG联合logFC.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
df2<- read.table("差异基因.xlsx",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
mydata <- data.frame(Category=df1$ONTOLOGY,
           ID=df1$ID,Term =df1$Description,
           Genes=df1$geneID,
#一定要把geneID列的/改成，
           adj_pval=df1$p.adjust)
genelist <- data.frame(ID=df2$ID,
                       logFC=df2$logFC)
cir <- circle_dat(mydata,genelist)
head(cir)
#write.csv(cir,"GOplot结果(OS).csv")

#GOCircle圈图
#圈图展示数据前十条通路
GOCircle(cir,rad1 = 2.5,#内圈直径
         rad2 = 3.5,#外圈直径
         label.size = 4,nsub = 10)#显示通路数目，也可以指定特定的ID(IDs <- c('GO:0007159','GO:0050863',...);nsub=IDs)
#GOCircle(cir, lfc.col = c('purple', 'orange'),label.size = 4)
#GOCircle(cir, zsc.col = c('yellow', 'black', 'cyan'),label.size = 4)
#GOCircle(cir, label.size = 4, label.fontface = 'italic')

#GOChord-GO弦图展示基因和通路之间的关系
termNum=8#限定基因数目，可以选择你感兴趣的通路
termNum=ifelse(nrow(mydata)<termNum,nrow(mydata),termNum)
geneNum=500#限定基因数目，可以选择你感兴趣的基因
#构建矩阵
chord <- chord_dat(cir, genelist[1:geneNum,],mydata$Term[1:termNum])

GOChord(chord, 
        space = 0.02, #基因之间的间距
        gene.order = 'logFC', #排序基因,"logFC","alphabetical","none"
        gene.space = 0.25, #基因离圆圈距离
        gene.size = 3.6,#基因字体大小
        border.size=0.1,#线的大小
        process.label=8)#GO名称大小




#维恩图####
library (VennDiagram) 
library(openxlsx)
#分组数据导入#
set1<-read.xlsx('C:/Users/Administrator/Desktop/Venn.xlsx',sheet= "Sheet1",sep=',')
set2<-read.xlsx('C:/Users/Administrator/Desktop/Venn.xlsx',sheet= "Sheet2",sep=',')
set3<-read.xlsx('C:/Users/Administrator/Desktop/Venn.xlsx',sheet= "Sheet3",sep=',')
set4<-read.xlsx('C:/Users/Administrator/Desktop/Venn.xlsx',sheet= "Sheet4",sep=',')
#数据转置，如果不转后头函数venn.diagram对矩阵数据不识别#
set1=t(set1)
set2=t(set2)
set3=t(set3)
set4=t(set4)
#二元#
venn1 = venn.diagram(x=list(set1,set2),
             scaled = F, # 根据比例显示大小
             alpha= 0.5, #透明度
             lwd=1,lty=1,col=c('black','black'), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             cex = 3, # 数字大小
             fontface = "bold",  # 字体粗细；加粗bold
             fill=c('blue','red'), # 填充色 配色https://www.58pic.com/
             category.names = c("DEG", "Estrogen-related genes") , #标签名
             cat.dist = 0.02, # 标签距离圆圈的远近
             cat.pos = -180, # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             cat.cex = 1.5, #标签字体大小
             cat.fontface = "bold",  # 标签字体加粗
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             cat.default.pos = "outer",  #标签位置, outer内;text 外
             output=TRUE,
             #filename='D:/Desktop/维恩图.png',# 文件保存
             filename = NULL,#自己手动保存
             imagetype="png",  # 类型（tiff png svg）
             resolution = 1000,  # 分辨率
             compression = "lzw"# 压缩算法
)
dev.off()#去图
grid.draw(venn1)


#三元#
venn2 = venn.diagram(x=list(set1,set2,set3),
             scaled = F, # 根据比例显示大小
             alpha= 0.6, #透明度
             lwd=2,lty=1,col=c('black'), #圆圈线条lwd粗细、lty形状、颜色；1 实线, 2 虚线, blank无线条
             label.col ='black' , # 数字颜
             cex = 1.6, # 数字大小
             fontface = "bold",  # 字体粗细；加粗bold
             fill=c("#FF0000", "#007FFF", "#00FF00"), # 填充色 配色https://www.58pic.com/
             category.names = c("GeneCard", "TRAP","JASPAR") , #标签名
             cat.dist = 0.05, # 标签距离圆圈的远近
             cat.pos = c(-20, 20, 180), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             cat.cex = 2, #标签字体大小
             cat.fontface = "bold",  # 标签字体加粗
             cat.col='black' , #根据相应颜色改变标签颜色
             cat.default.pos = "outer",  # 标签位置, outer外;text内
             filename = NULL)
dev.off()#去图
grid.draw(venn2)

#四元#
venn3 =venn.diagram(x=list(set1,set2,set3,set4),
             scaled = T, # 根据比例显示大小
             alpha= 0.6, #透明度
             lwd=1,lty=1,col=c("black"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             cex = 1.8, # 数字大小
             fontface = "bold",  # 字体粗细；加粗bold
             fill=c('#DC143C',"#7FFFAA",'#1E90FF',"#9400D3"), # 填充色 配色https://www.58pic.com/
             category.names = c("IGF2BP3-RIP", "sh-IGF2BP3","MeRIP","Survival") , #标签名
             cat.dist = c(0.2, 0.2, 0.1, 0.1), # 标签距离圆圈的远近
             cat.pos = c(-10, 10, -10, 10), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             cat.cex = 1.5, #标签字体大小
             cat.fontface = "bold",  # 标签字体加粗
             cat.col=c('black'),   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             filename = NULL)
dev.off()#去图
grid.draw(venn3)



#箱型图####
library(ggplot2)
library(tidyr)
library(ggsci)
library(ggpubr)
library(ggrepel)
library(tidyverse)
#data=melt(b,id.vars=c("Type"))
#colnames(data)=c("Type","Gene","Expression")
#这两句可以把TCGA数据变箱型图格式

#####1肿瘤瘤周####
#全部患者
exp <- read.table("TCGA肝癌mRNA_counts表达矩阵.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#提取hub基因
data<- read.table("箱型图基因.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
colnames(data)
markers<- colnames(data)[3:3]#批量计算那些列
exp1 <- exp[markers,]#提取基因子集
exp1 = exp1 %>% t() %>% as.data.frame()#倒置，变表格
exp1<-exp1%>%mutate(Group=factor(ifelse(str_detect(rownames(exp1), "01A"), "Tumor","Normal")))
exp1 <- exp1 %>% rownames_to_column("sample")#新增一列group
b <- gather(exp1,key=Genes,value = exp,-c(Group,sample))
#方法1
ggboxplot(b, x = "Genes", y = "exp",fill = "Group", palette = "lancet")+
  stat_compare_means(aes(group = Group),method = "wilcox.test",
  label = "p.signif",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
  symbols = c("***", "**", "*", "ns")))+theme(text = element_text(size=10),
  axis.text.x = element_text(angle=45, hjust=1)) 
#方法2
ggplot(b, aes(x=Genes, y=exp, fill=Group, palette = "jco")) + #画图
  geom_boxplot()+
  stat_compare_means(label = "p.signif", method = "wilcox.test")+ #选择检验算法
  ggsci::scale_fill_lancet()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

write.table(b, file = "箱型图数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
#b<- read.table("箱型图数据.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)




#####2TNM分期####
#仅分析肿瘤患者
exp <- read.table("TCGA+GEO表达矩阵.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
data<- read.table("箱型图基因.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)#这里把TNM分期总结好
colnames(exp) <- substr(colnames(exp),1,12)#保留列名前12位
rownames(data) <- substr(rownames(data),1,12)#保留列名前12位
colnames(data)#注意要哪些基因
markers<- colnames(data)[3:3]#批量计算那些列
A <- exp[markers,]#提取基因子集
A=t(A)
comgene <- intersect(rownames(A),rownames(data))#筛选肿瘤患者
A<- A[comgene,]#A为表达矩阵
B<- data[comgene,]#B为临床信息
identical(rownames(A),rownames(B))
A=as.data.frame(A)#倒置，变表格
A$Group<-B$FIGO#新增一列,注意这里一定不能是数字
exp1 <- A %>% rownames_to_column("sample")
b <- gather(exp1,key=Genes,value = exp,-c(Group,sample))
#方法1
ggboxplot(b, x = "Genes", y = "exp",fill = "Group", palette = "lancet")+#jco
  stat_compare_means(aes(group = Group),method = "wilcox.test",#多组kruska-wallis
  label = "p.signif",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
  symbols = c("***", "**", "*", "ns")))+theme(text = element_text(size=10),
  axis.text.x = element_text(angle=45, hjust=1)) 




####——免疫——####
#开始前先求出fpkm-01A的数据：
counts <- read.table("TCGA肝癌表达矩阵counts01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#fpkm <- read.table("TCGA肝癌mRNA_fpkm表达矩阵.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group_list <- ifelse(as.numeric(str_sub(colnames(fpkm),14,15)) < 10,'tumor','normal')
table(group_list)
fpkm01A <- fpkm[,group_list == c("tumor")]#筛选表达矩阵中的肿瘤样本：
fpkm01A= data.frame(t(fpkm01A))
fpkm01A$name=substr(rownames(fpkm01A),1,12)#保留列名前12位
fpkm01A<- fpkm01A[!duplicated(fpkm01A$name),]#去重复
rownames(fpkm01A) <- fpkm01A$name
fpkm01A<- fpkm01A[,-ncol(fpkm01A)]   #去除最后一列
write.table(fpkm01A, file = "TCGA肝癌表达矩阵fpkm01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
A=t(expr)
write.table(A, file = "TCGA肝癌表达矩阵fpkm01A(倒置).txt",sep = "\t",row.names = T,col.names = NA,quote = F)


####1-ESTIMATE免疫评分与肿瘤纯度#####
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
library(utils) 
library(estimate)
library(tidyverse)
library(utils)
#读取肿瘤患者01A表达谱(19938*185)必须是这种，可能需要倒置
expr <- read.table("宫颈癌表达矩阵counts01A(倒置).txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#计算免疫评分
filterCommonGenes(input.f = "宫颈癌表达矩阵counts01A(倒置).txt",#文件名必须和上面一样
                  output.f = "肝癌免疫评分.gct",   #输出文件名
                  id = "GeneSymbol")   #行名为gene symbol
estimateScore("肝癌免疫评分.gct",   #刚才的输出文件名
              "肝癌免疫评分.txt",   #新的输出文件名（即估计的结果文件）
              platform="affymetrix") #默认平台

#3. 输出每个样品的打分
result <- read.table("肝癌免疫评分.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
result <- result[,-1]   
colnames(result) <- result[1,]   
result <- as.data.frame(t(result[-1,]))
rownames(result) <- colnames(expr)
write.table(result, file = "宫颈癌免疫评分.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F) # 保存并覆盖得分
#得到结果后自己用Prism作图



####2-Cibersort免疫####
library(e1071)
library(parallel)
library(preprocessCore)
source("CIBERSORT.R")#注释文件，自己下载 
sig_matrix <- "LM22.txt"#注释文件，自己下载 
mixture_file = '宫颈癌表达矩阵counts01A(倒置).txt' #(19938*185)
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)
res_cibersort <- res_cibersort[,1:22]#耗时   
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]   #去除丰度全为0的细胞
#save(res_cibersort,file = "res_cibersort.Rdata")   #保存中间文件
#write.table(ciber.res,"ciber.res.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#可视化（彩虹图）
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7) #创建彩虹色板（带70%透明度）
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(ciber.res)),
        border = NA, # 柱子无边框写
        names.arg = rep("",nrow(ciber.res)), # 无横坐标样本名
        yaxt = "n", # 先不绘制y轴
        ylab = "Relative percentage", # 修改y轴名称
        col = mycol) # 采用彩虹色板
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # 补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-1,par("usr")[4],legend = colnames(ciber.res),xpd = T,
       fill=mycol,cex=0.8,y.intersp=1.25,x.intersp=1,border=1,bty = "n")


#####分组比较####
library(ggsci)
library(tidyr)
library(ggpubr)
library(tidyverse)
#方法1，自己分组
a <-as.data.frame(ciber.res)
exp <- read.table("肝癌临床信息01A.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
table(exp$A)
exp <- exp%>% mutate(group=#新增一列
     factor(ifelse(str_detect(exp$B, "1"),"High","Low"),levels = c("Low","High")))
#将cl样本行名的顺序按exp2列名进行匹配：
exp <- exp[match(rownames(a),rownames(exp)),]
identical(rownames(a),rownames(exp))

#方法2，根据基因分组
exp <- read.table("宫颈癌表达矩阵counts01A(倒置).txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- exp %>% t() %>% as.data.frame()
gene="CD274"
med=median(as.numeric(exp[,gene]))#取中位数
exp <- exp %>% mutate(Group=factor(ifelse(exp[,gene]<med,"High","Low"),levels = c("Low","High")))
a <-as.data.frame(ciber.res)#备份
identical(rownames(a),rownames(exp))#.和-的差别

#重点是顺序一致
a$Group <- exp$Group
a <- a %>% rownames_to_column("sample")
b <- gather(a,key=CIBERSORT,value = Fraction,-c(Group,sample))
ggboxplot(b, x = "CIBERSORT", y = "Fraction",
  fill = "Group", palette = "jco")+ #lancet or jco
  stat_compare_means(aes(Group = Group),
  method = "wilcox.test",
  label = "p.signif",
  symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
  symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
  axis.text.x = element_text(angle=45, hjust=1)) 


#####基因与cibersort相关性热图####
library(corrplot)
library(tidyverse)
exp = read.table("宫颈癌表达矩阵counts01A(倒置).txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <-c("CD4","CD8A","ITGAX")
exp <- exp[,gene]
ciber = read.table("ciber.res.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
identical(rownames(ciber),rownames(exp))

cor<-sapply(ciber,function(x,y) cor(x,y,method="spearman"),exp)
rownames(cor)<-colnames(exp)
cor_res <- cor.mtest(cor,计算p值conf.level = 0.95)#置信区间
corrplot(cor,
         method = "color",#相关性矩阵展示的图形
         col=colorRampPalette(c("#01468b","white","#ee0000"))(100),
         addCoef.col = "black",#为相关系数添加颜色
         tl.col="black",#设置文本标签的颜色
         number.cex = 0.5,
         tl.cex = 0.7,
         cl.align = "l")


####基因与cibersort相关性散点图
exp = read.table("宫颈癌表达矩阵counts01A(倒置).txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- c("G6PD","STMN1","TMMR")
exp <- exp[gene,]
#exp <- exp %>% t() %>% as.data.frame()
ciber = read.table("ciber.res.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ciber <- ciber %>% t() %>% as.data.frame()
rownames(ciber) <- gsub(" ",".",rownames(ciber))
identical(colnames(ciber),colnames(exp))
exp_ciber <- rbind(exp,ciber)
exp_ciber <- exp_ciber %>% t() %>% as.data.frame()
update.packages("rlang")
library(ggstatsplot)
library(tidyverse)#这里只能一个基因一个基因做
ggscatterstats(data = exp_ciber, #要分析的数据
               y = NRP1, #设置Y轴
               x = B.cells.naive,#设置X轴
               type = "nonparametric", 
               margins = "both",#是否显示 边缘，默认为true                                      
               xfill = "#01468b", #x轴边缘图形的颜色
               yfill = "#ee0000", #y轴边缘图形的颜色
               marginal.type = "densigram")#在图片坐标轴边缘添加图形类型



####基因与基因相关性散点图
exp = read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#install.packages("ggstatsplot")
library(ggstatsplot)
library(tidyverse)
exp <- exp %>% t() %>% as.data.frame()
ggscatterstats(data = exp, #要分析的数据
               y = CDKN3, #设置Y轴
               x = PDCD1,#设置X轴
               type = "nonparametric", 
               margins = "both",#是否显示 边缘，默认为true                                      
               xfill = "#01468b", #x轴边缘图形的颜色
               yfill = "#ee0000", #y轴边缘图形的颜色
               marginal.type = "densigram")#在图片坐标轴边缘添加图形类型


####基因与ESTIMATE相关性散点图
exp = read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
est = read.table("ESTIMATE_result.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- exp %>% t() %>% as.data.frame()

identical(rownames(est),rownames(exp))
exp_est <- cbind(exp,est)
library(ggstatsplot)
library(tidyverse)
ggscatterstats(data = exp_est, #要分析的数据
               y = FGR, #设置Y轴
               x = ImmuneScore,#设置X轴
               type = "nonparametric", 
               margins = "both",#是否显示 边缘，默认为true                                      
               xfill = "#01468b", #x轴边缘图形的颜色
               yfill = "#ee0000", #y轴边缘图形的颜色
               marginal.type = "densigram")#在图片坐标轴边缘添加图形类型





##3-ssGSEA免疫####
library(tidyverse)
library(data.table)
library(GSVA)
library(ggsci)
library(tidyr)
library(ggpubr)
#1准备细胞marker
cellMarker <- data.table::fread("cellMarker.csv",data.table = F)
colnames(cellMarker)[2]<-"celltype"#将cellMarker文件列名的第2个修改为celltype
type <- split(cellMarker,cellMarker$celltype)#将cellMarker文件以celltype为分组拆分成list数据格式
#处理data.tables列表通常比使用group by参数按组对单个data.table进行操作要慢得多
cellMarker <- lapply(type, function(x){dd = x$Metagene
  unique(dd)})#将list中每个celltype中的基因进行合并



##2表达量矩阵的准备(19938*185)
exp1<- read.table("宫颈癌表达矩阵counts01A(倒置).txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp1 <- as.matrix(exp1) #将expr转换为矩阵格式
#gsva_data <- gsva(expr,cellMarker, method = "ssgsea")#4.3版本之前能用
gsvapar <- gsvaParam(exp1, cellMarker, maxDiff=TRUE) 
ES <- gsva(gsvapar)

#####分组比较####
#方法1，自己分组
exp <- read.table("肝癌临床信息01A.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
table(exp$A)
exp <- exp%>% mutate(group=factor(ifelse(str_detect(exp$B, "1"),"high","low"),levels = c("low","high")))
exp <- exp[match(rownames(b),rownames(exp)),]
identical(rownames(b),rownames(exp))

#方法2，根据基因分组
exp <- exp1 %>% t() %>% as.data.frame(stringsAsFactors = F)
gene="CD274"
med=median(as.numeric(exp[,gene]))#取中位数
exp <- exp %>% mutate(group=factor(ifelse(exp[,gene]>med,"High","Low"),levels = c("Low","High")))
b <- ES %>% t() %>% as.data.frame(stringsAsFactors = F)
identical(rownames(b),rownames(exp))


#4使用ssGSEA量化免疫浸润
b$group <- exp$group#B
b <- b %>% rownames_to_column("sample")
b <- gather(b,key=ssGSEA,value = Expression,-c(group,sample))
ggboxplot(b, x = "ssGSEA", y = "Expression",
          fill = "group", palette = "jco")+
          stat_compare_means(aes(group = group),
          method = "wilcox.test",
          label = "p.signif",
          symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
          symbols = c("***", "**", "*", "ns")))+
          theme(text = element_text(size=10),
          axis.text.x = element_text(angle=45, hjust=1)) 

#write.table(a,"ssGSEA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)




####4-X cell免疫####
#devtools::install_github('dviraran/xCell') #包的安装,致谢作者
library(xCell)
library(ggpubr)
library(tidyverse)
library(ggsci)
library(tidyr)
library(ggpubr)
library(GSVA)
#读取表达矩阵(19938*185)
exp <- read.table("宫颈癌表达矩阵counts01A(倒置).txt", sep = "\t",row.names = 1,check.names = F,header = T)
celltypeuse<-xCell.data$signatures
rs<-xCellAnalysis(exp,parallel.sz=10) #3.5版本之前可用


#####分组比较
#方法1，自己分组
exp <- read.table("肝癌临床信息01A.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
table(exp$B)
exp <- exp%>% mutate(group=factor(ifelse(str_detect(exp$B, "1"),"High","Low"),levels = c("Low","High")))
exp <- exp[match(rownames(rs),rownames(exp)),]
b <- rs %>% t() %>% as.data.frame(stringsAsFactors = F)
identical(rownames(b),rownames(exp))

#方法2，根据基因分组
exp <- exp %>% t() %>% as.data.frame()
gene="PDCD1"
med=median(as.numeric(exp[,gene]))#取中位数
exp <- exp %>% mutate(group=factor(ifelse(exp[,gene]>med,"High","Low"),levels = c("Low","High")))
b <- rs %>% t() %>% as.data.frame(stringsAsFactors = F)
identical(rownames(b),rownames(exp))

b$group <- exp$group#B#group
b <- b %>% rownames_to_column("sample")
b <- gather(b,key=xCell,value = Expression,-c(group,sample))
ggboxplot(b, x = "xCell", y = "Expression",fill = "group", palette = "jco")+
  stat_compare_means(aes(group = group),method = "wilcox.test",label = "p.signif",
  symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
  symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45, hjust=1)) 


#####新的方法#####
#如果上面失败了采用这个方法
exp1 <- as.matrix(exp) #将expr转换为矩阵格式#将cellMarker
Es <- ssgseaParam(exp1, celltypeuse, assay = NA_character_,
                  annotation = NA_character_,minSize = 1,maxSize = Inf,alpha = 0.25,normalize = TRUE)
rs<- gsva(Es)
#1先把原始数据方法网页计算http://xcell.ucsf.edu/(13000*180格式)
#2然后复制到Xcell里,字母顺序排序，最后三个是免疫、基质和微环境评分
b <- read.table("xCell.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#自己分组,NC组为HIGH
b <- b %>% mutate(group=factor(ifelse(str_detect(b$group, "High"),"high","low"),levels = c("low","high")))
b <- b %>% rownames_to_column("sample")
b <- gather(b,key=xCell,value = Expression,-c(group,sample))
ggboxplot(b, x = "xCell", y = "Expression",fill = "group", palette = "jco")+
  stat_compare_means(aes(group = group),method = "wilcox.test",label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45, hjust=1)) 




####5-免疫检查点####
library(tidyverse)
library(ggsci)
library(tidyr)
library(ggpubr)
#1分组
#方法1，自己分组
exp <- read.table("TCGA肝癌表达矩阵counts01A(倒置).txt", sep = "\t",row.names = 1,check.names = F,header = T)
exp1 <- read.table("肝癌临床信息01A.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- exp %>% t() %>% as.data.frame()
table(exp1$A)
exp1 <- exp1%>% mutate(group=factor(ifelse(str_detect(exp1$B, "1"),"high","low"),levels = c("low","high")))
exp <- exp[match(rownames(exp1),rownames(exp)),]
identical(rownames(exp1),rownames(exp))
exp$group<-exp1$group

#方法2，根据基因分组(19938*185)
exp <- read.table("宫颈癌表达矩阵counts01A(倒置).txt", sep = "\t",row.names = 1,check.names = F,header = T)
gene="CD274"
med=median(as.numeric(exp[gene,]))
exp <- exp %>% t() %>% as.data.frame()
exp <- exp %>% mutate(group=factor(ifelse(exp[,gene]<med,"Low","High"),levels = c("Low","High")))
class(exp$group)

#2画图,自行导入“免疫检查点基因”文件取交集
Immucheckpoints<- read.table("免疫检查点基因.xlsx", sep = "\t",check.names = F,header = T)
comgene <- intersect(colnames(Immucheckpoints),colnames(exp))
a <- exp[,comgene]
a <- a %>% rownames_to_column("sample")
b <- gather(a,key=ICB,value = Expression,-c(group,sample))
ggboxplot(b, x = "ICB", y = "Expression",fill = "group", palette = "jco")+
  stat_compare_means(aes(group = group),method = "wilcox.test",
  label = "p.signif",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
  symbols = c("***", "**", "*", "ns")))+theme(text = element_text(size=10),
  axis.text.x = element_text(angle=45, hjust=1)) 



####6-TMB肿瘤突变负荷####
library(TCGAbiolinks)
library(maftools)
library(tidyverse)
library(readxl)
library(readr)
library(OncogenicPathways)
#0，TCGA突变数据下载，网址：https://portal.gdc.cancer.gov/
query <- GDCquery(
  project = "TCGA-CESC", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open")
GDCdownload(query)
GDCprepare(query, save = T,save.filename = "TCGA-LIHC_SNP.Rdata")

#1，下次直接从这开始，读取数据
load(file = "./TCGA-LIHC_SNP.Rdata")

maf.coad <- data
class(maf.coad)
dim(maf.coad)
maf.coad[1:10,1:10]
maf <- read.maf(maf.coad)

#2，画图，总的
coad.tmb <- tmb(maf, captureSize = 50, logScale = T)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
#绘图，分开的
oncoplot(maf = maf,
         top = 10, #显示前30个的突变基因信息
         fontSize = 0.8,   #设置字体大小
         showTumorSampleBarcodes = F)   #不显示病人信息

#3，Oncostrip绘制 某些关键基因
oncostrip(maf = maf, fontSize = 0.8, showTumorSampleBarcodes = F,
      genes = c('TTN','TP53', 'MUC16'))#改你想要的基因就行

#4，titv,函数绘制SNP分类
maf.titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = maf.titv)# plot titv summary

#5，lollipop plot棒糖图显示蛋白质结构突变点
lollipopPlot(maf = maf, gene = 'TP53', AACol = 'RefSeq', showMutationRate = TRUE)

#6，rainfallPlot线性基因组尺度上绘制变异间距离图来可视化超突变的基因组区域
rainfallPlot(maf = maf, detectChangePoints = TRUE, pointSize = 0.6)

#7，可以比较给定MAF中的突变负荷与TCGA数据库中不同癌症队列的中位突变载量来进行比较较
maf.mutload = tcgaCompare(maf = maf, cohortName = 'Example-Maf')

#8，plotVaf可以将不同的等位基因频率绘制为箱式图
plotVaf(maf = maf)




####——聚类——####
####1ConsensusClusterPlus聚类分析####
#install.packages("BiocManager")
#BiocManager::install('ConsensusClusterPlus')
library(tidyverse)
library(ConsensusClusterPlus)
exp <- read.table("TCGA肝癌表达矩阵fpkm01A(倒置).txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
d=as.matrix(exp)
gene <- c("G6PD","STMN1","HMMR") #"ALKBH5","YTHDF1"，必须是多个基因聚类
d <- d[gene,]
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:3],] #3是基因个数，需要改
d = sweep(d,1, apply(d,1,median,na.rm=T))

title=("JULEI") ##文件夹输出图片的位置
set.seed(1) #我发现设不设置种子都一样
results = ConsensusClusterPlus(d,maxK=9,reps=50,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",seed=1,plot="pdf")
results[[2]][["consensusMatrix"]][1:5,1:5]
results[[2]][["consensusTree"]]
results[[2]][["consensusClass"]][1:5]
icl = calcICL(results,title=title,plot="pdf") ##画另一组图片

group<-results[[2]][["consensusClass"]]
group<-as.data.frame(group)
group$group <- factor(group$group,levels=c(1,2))
save(group,file = "group_AY.Rda")
load("group_GENE23456FINAL.Rda")

exp_gene <- exp[gene,]

# 绘制ConsensusClusterPlus后的热图
library(pheatmap)
group <- group %>% 
  rownames_to_column("sample")
annotation <- group %>% arrange(group) %>% column_to_rownames("sample")
a <- group %>% arrange(group) %>% mutate(sample=substring(.$sample,1,12))
b <- t(exp_gene) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  mutate(sample=substring(.$sample,1,12))
c <- inner_join(a,b,"sample") %>% .[,-2] %>% column_to_rownames("sample") %>% t(.)
pheatmap(c,annotation = annotation,
         cluster_cols = F,fontsize=5,fontsize_row=5,
         scale="row",show_colnames=F,
         fontsize_col=3)
pheatmap(c,annotation = annotation,
         annotation_colors = list(group = c("1" ="#01468b","2"= "#ee0000")),
         cluster_cols = F,fontsize=5,fontsize_row=5,
         scale="row",show_colnames=F,cluster_row = F,
         fontsize_col=3)
dev.off()

#小提琴图###
setwd("cibersort_AY")
library(tidyverse)
a <- read.table("CIBERSORT-Results.txt", sep = "\t",row.names = 1,check.names = F,header = T)
a <- a[,1:22]
identical(rownames(a),rownames(group))
b <- group
class(b$group)
a$group <- b$group
a <- a %>% rownames_to_column("sample")
library(ggsci)
library(tidyr)
library(ggpubr)
#install.packages("ggsci")
#install.packages("tidyr")
#install.packages("ggpubr")

b <- gather(a,key=CIBERSORT,value = Proportion,-c(group,sample))

ggboxplot(b, x = "CIBERSORT", y = "Proportion",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 
dev.off()





####2WGCNA共表达聚类####
#加权基因共表达网络分析，鉴定表达模式相似的基因集合（module）。
#解析基因集合与样品表型之间的联系，绘制基因集合中基因之间的调控网络并鉴定关键调控基因。
#install.packages("WGCNA")
#install.packages("BiocManager")
#BiocManager::install("preprocessCore")
#BiocManager::install("impute")
#if(!require(DESeq2))BiocManager::install('DESeq2')
library("tidyverse")
library("WGCNA")            
library(DESeq2)
exp <- read.table("宫颈癌表达矩阵counts01A(倒置).txt",sep = "\t",check.names = F,
        row.names = 1,stringsAsFactors = F,header = T)#（11920*371）
#exp1=t(exp)
#可以用总的数据去聚类，或者用差异基因聚类也行
#比如我先把差异基因读进来，再取交集
A <- read.table("Radscore差异基因.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
comgene <- intersect(rownames(exp),rownames(A))
exp <- exp[comgene,]
#exp =exp [,1:43]
#开始聚类
datExpr0 <- exp %>% t() %>% as.data.frame()#保证样本在行，基因在列（371*11920）
##开始WGCNA
gsg = goodSamplesGenes(datExpr0, verbose = 3)#检查缺失值
gsg$allOK
###如果没有达标就需要筛选
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]}

# 样品聚类
sampleTree = hclust(dist(datExpr0), method = "average")
# 画图
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 200, col = "red")###剪切线
###删除剪切线以下的样品
clust = cutreeStatic(sampleTree, cutHeight = 200, minSize = 10)
table(clust)
keepSamples = (clust==1)#保留1簇
datExpr0 = datExpr0[keepSamples, ]
# 重新聚类
sampleTree2 = hclust(dist(datExpr0), method = "average")
plot(sampleTree2)
# 记录基因和样本数，方便后续可视化
nGenes = ncol(datExpr0)#基因数
nSamples = nrow(datExpr0)#样本数
#save(datExpr0, nGenes, nSamples,file = "Step01-WGCNA_input.Rda")

# 构建网络，识别模块，power值散点图
enableWGCNAThreads()   #多线程工作
powers = c(1:20)       #幂指数范围1:20
sft=pickSoftThreshold(datExpr0,powerVector=powers,verbose=5)#耗时
par(mfrow = c(1,2))
###拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col="red");
abline(h=0.90,col="red") #可以修改
###平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels=powers, cex=1,col="red")

##邻接矩阵转换
sft #查看最佳power值
softPower =sft$powerEstimate #最佳power值
softPower
adjacency = adjacency(datExpr0, power = softPower)#耗时

##TOM矩阵
TOM = TOMsimilarity(adjacency)#耗时
dissTOM = 1-TOM
#save(TOM,file = "TOM.Rda")

# 基因聚类
geneTree = hclust(as.dist(dissTOM), method = "average");
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# 动态剪切模块识别
minModuleSize = 50      #模块基因数目，这里数字越大模块越少
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# 相似模块聚类
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

###相似模块合并
MEDissThres = 0.1 #剪切高度可修改
abline(h=MEDissThres, col = "red")
#如果不需要删，聚类就结束了
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
plotDendroAndColors(geneTree, mergedColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs


#####整理临床信息####
clinical <- read.table("肝癌临床信息01A.xlsx",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#这里用临床信息/免疫评分什么参数都可以，但是列名字必须一样TCGA-5R-AACC
clinical = datExpr0[,0]#新增一个空表
clinical$Radscore=datExpr0$CD4
clinical$Pathscore=datExpr0$CD274
identical(rownames(clinical),rownames(datExpr0))
head(clinical)# 查看临床信息
# 对表达矩阵进行预处理
datTraits = as.data.frame(do.call(cbind,lapply(clinical, as.numeric)))
rownames(datTraits) = rownames(clinical)
# 对样本进行聚类
sampleTree2 = hclust(dist(datExpr0), method = "average")
# 将临床信息转换为颜色，白色表示低，红色表示高，灰色表示缺失
traitColors = numbers2colors(datTraits, signed = FALSE)
# 样本聚类图与样本性状热图
plotDendroAndColors(sampleTree2, 
                    traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")

####网络的分析
MEs=orderMEs(MEs)# 对模块特征矩阵进行排序
#计算模型特征矩阵和样本信息矩阵的相关度。
moduleTraitCor=cor(MEs, datTraits, use="p")
#write.table(file="Step04-modPhysiological.cor.xls",moduleTraitCor,sep="\t",quote=F)
moduleTraitPvalue=corPvalueStudent(moduleTraitCor, nSamples)
#write.table(file="Step04-modPhysiological.p.xls",moduleTraitPvalue,sep="\t",quote=F)

#使用labeledHeatmap()将上述相关矩阵和p值可视化。
textMatrix=paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
# 基因模块与临床信息相关性图
labeledHeatmap(Matrix=moduleTraitCor,#模块和表型的相关性矩阵，这个参数最重要，其他可以不变
               xLabels=colnames(datTraits),
               yLabels=names(MEs),
               ySymbols=names(MEs),
               colorLabels=F,
               colors=blueWhiteRed(100),#颜色种类
               textMatrix=textMatrix,
               setStdMargins=T,
               cex.text=0.45,#里面字体大小
               cex.lab=0.7,#外面字体大小
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))

# 不同模块与基因性状的具体分析
##矩阵一
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
####看一下目的基因和哪个模块相关性最高
a <- geneModuleMembership
a <- a %>% rownames_to_column()

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

##矩阵二
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

##批量输出性状和模块散点图
for (trait in traitNames){
  traitColumn=match(trait,traitNames)  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      outPdf=paste(trait, "_", module,".pdf",sep="")
      pdf(file=outPdf,width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      abline(v=0.8,h=0.5,col="red")
      dev.off()
    }
  }
}

#10. 输出每个模块的基因
for (mod in 1:nrow(table(moduleColors)))
{  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0("基因",modules,".xlsx"),sep="\t",row.names=F,col.names=F,quote=F)
}




####—其他RNA—####
dev.off()
rm(list = ls())
library(tidyverse)
library(GEOquery)
library(stringr)
library(limma) 
####1、CircalRNA####
###下载数据，如果文件夹中有会直接读入
gset = getGEO('GSE184882', destdir=".", AnnotGPL = F, getGPL = F)
class(gset)
###提取子集
gset[[1]]
#通过pData函数获取分组信息
pdata <- pData(gset[[1]])
# 通过exprs函数表达矩阵
exp <- exprs(gset[[1]])
#设置参考水平
table(pdata$source_name_ch1)#注意是不是source_name_ch1
group_list <- ifelse(str_detect(pdata$source_name_ch1, "non-tumor"), "normal",
                     "tumor")
#pdata<- pdata[-(1:20),]   #去除第41-72行
#exp<- as.data.frame(exp)#变表格
#exp<- exp[,-(1:20)]
table(pdata$characteristics_ch1.4)
group_list <- ifelse(str_detect(pdata$characteristics_ch1.4, "Primary tumor"), "tumor",
                     "normal")
#group_list <- ifelse(str_detect(pdata$characteristics_ch1.1, "non-tumor"), "normal","tumor")
#因子型
group_list = factor(group_list,
                    levels = c("normal","tumor"))
table(group_list)

#画图
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
###数据校正（差异分析前进行校正是有必要的，但会减少上调基因）
exp=normalizeBetweenArrays(exp)
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
range(exp)
write.table(exp, file = "表达矩阵.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


#标准化
ex <- expB
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0,arr.ind = T)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
exp <- exprSet#如果不需要转换就显示没exprSet，是对的。
write.table(exp, file = "表达矩阵标准化.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


#PCA
library(FactoMineR)
library(factoextra)
dat=as.data.frame(t(exp))
dat.pca <- PCA(dat, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = group_list, # color by groups
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups")
pca_plot
#本来应该先进行ID转换的，但CircalRNA比较特殊
#所以先进行差差异分析

#差异分析
library(limma)
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
##标记上下调基因
logFC=1
P.Value = 0.05
k1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
write.table(deg, file = "上调基因.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
#手动挑选上调基因做成上调基因txt



####ID转换####
rm(list = ls())
#官网下载GPL探针，然后与上调基因取交集
#导入上调基因时注意head选yes + row name选fisrt coulm
A <- `上调基因`
B <- GPL10553
comgene <- intersect(rownames(A),rownames(B))
A <- A[comgene,]
B <- B[comgene,]
C <- cbind(B,A)
write.table(C, file = "上调基因6字名字.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
#导入上调基因6字名字和7字探针
rm(list = ls())
A <- `上调基因6字名字`
B <- `七字探针`
comgene <- intersect(rownames(A),rownames(B))
A <- A[comgene,]
B <- B[comgene,]
C <- cbind(B,A)
write.table(C, file = "上调基因7字名字.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
#取交集
#导入上调基因7字名字和全人类基因
rm(list = ls())
A <- `上调基因7字名字`
B <- `全人类基因`
comgene <- intersect(rownames(A),rownames(B))
A <- A[comgene,]
B <- B[comgene,]
C <- cbind(A,B)
write.table(C, file = "交集基因.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)





####2、LncRNA####
dev.off()
rm(list = ls())
#采集GEO数据
###加载R包
library(tidyverse)
library(GEOquery)
library(stringr)
library(limma) 
###下载数据，如果文件夹中有会直接读入
gset = getGEO('GSE99417', destdir=".", AnnotGPL = F, getGPL = F)
class(gset)
###提取子集
gset[[1]]
#通过pData函数获取分组信息
pdata <- pData(gset[[1]])
# 通过exprs函数表达矩阵
exp <- exprs(gset[[1]])
#设置参考水平
table(pdata$title)#注意是不是source_name_ch1
group_list <- ifelse(str_detect(pdata$title, "Gastric cancer tissues"), "tumor",
                     "normal")

table(pdata$characteristics_ch1)
group_list <- ifelse(str_detect(pdata$characteristics_ch1, "Tumor"), "tumor",
                     "normal")
#group_list <- ifelse(str_detect(pdata$characteristics_ch1.1, "non-tumor"), "normal","tumor")
#因子型
group_list = factor(group_list,
                    levels = c("normal","tumor"))
table(group_list)
write.table(exp, file = "表达矩阵.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)

#标准化表达矩阵
ex <- expB
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0,arr.ind = T)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
exp <- exprSet#如果不需要转换就显示没exprSet，是对的。
#write.table(exp, file = "表达矩阵标准化.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


#画图
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
###数据校正（差异分析前进行校正是有必要的，但会减少上调基因）
exp=normalizeBetweenArrays(exp)
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
range(exp)

#PCA
library(FactoMineR)
library(factoextra)
dat=as.data.frame(t(exp))
dat.pca <- PCA(dat, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = group_list, # color by groups
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups")
pca_plot
#本来应该先进行ID转换的，但CircalRNA比较特殊
#所以先进行差差异分析

#差异分析
library(limma)
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
##标记上下调基因
logFC=1
P.Value = 0.05
k1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
write.table(deg, file = "上调基因.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
#手动挑选上调基因做成上调基因txt


####ID转换####
rm(list = ls())
#官网下载GPL探针，然后与上调基因取交集
#导入上调基因时注意head选yes + row name选fisrt coulm
A <- `上调基因`
B <- GPL16956
comgene <- intersect(rownames(A),rownames(B))
A <- A[comgene,]
B <- B[comgene,]
C <- cbind(B,A)
write.table(C, file = "上调基因Ensembl.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
#GPL转换为Ensembl-gene-ID,如果没转换成，建议从DAVID转化
#然后从新整合到一起
rm(list = ls())
A <- `上调基因Ensembl`
B <- `上调基因EnsemblDAVID`
comgene <- intersect(rownames(A),rownames(B))
A <- A[comgene,]
B <- B[comgene,]
C <- cbind(A,B)
write.table(C, file = "上调基因Ensembl.DAVID.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
#导入上调基因7字名字和全人类基因
rm(list = ls())
A <- `m6A修饰的lncRNA`
B <- `上调基因汇总`
rownames(B) <- substr(rownames(B),1,15)#保留列名前15位
comgene <- intersect(rownames(A),rownames(B))
A <- A[comgene,]
B <- B[comgene,]
C <- cbind(A,B)
write.table(C, file = "交集基因铁死亡.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


####3、原始数据找表达矩阵####
#先用常规方式试一下能不能提出表达矩阵，同时也赋值pd
gset = getGEO('GSE61797', destdir=".", AnnotGPL = F, getGPL = F)
class(gset)
gset[[1]]#0个特征
pd <- pData(gset[[1]])
exp <- exprs(gset[[1]])#0
####Illumina####
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("lumi", quietly = TRUE)) BiocManager::install("lumi")
if (!requireNamespace("methylumi", quietly = TRUE)) BiocManager::install("methylumi")
BiocManager::install("GenomicFeatures")
BiocManager::install("methylumi")
BiocManager::install("runway")
#remotes::install_github('ML4LHS/runway')

library(runway)
library(methylumi)
library(lumi)
library(limma)
fileName <- 'GSM1513948_Control.GeneExpression.txt' # 自己手动下载，放到路径中,记得加.gz
x.lumi <- lumiR.batch(fileName) # subscript out of bounds下标出界,代表不正常
pData(phenoData(x.lumi))
lumi.N.Q <- lumiExpresso(x.lumi)
### retrieve normalized data
dataMatrix <- exprs(lumi.N.Q) #dataMatrix即表达矩阵
exp = dataMatrix#即可无缝对接pipeline 01
boxplot(exp) #看一下表达量数据分布是否一致
write.table(dataMatrix, file = "表达矩阵dataMatrix.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


#不正常的
fileName <- 'GSM1513948_Control.GeneExpression.txt' 
a=read.table(fileName,header = T,sep = '\t')
b=a #备份 a=b
colnames(a)=c('PROBE_ID',paste(rownames(pd),
                               rep(c('AVG_Signal','Detection Pval'),9)#数字需要改
                               ,sep='.'))
colnames(a)
write.table(a,file = 'raw data.txt',sep = '\t',quote = F)
x.lumi <- lumiR.batch('raw data.txt') ##, sampleInfoFile='sampleInfo.txt')
pData(phenoData(x.lumi))   
lumi.N.Q <- lumiExpresso(x.lumi)
dat <- exprs(lumi.N.Q)
dim(dat)#看一下dat这个矩阵的维度
dat[1:4,1:4] #查看dat这个矩阵的1至4行和1至4列，逗号前为行，逗号后为列
boxplot(dat[ ,1:4] ,las=2) 
write.table(dat, file = "表达矩阵GSE174237_Raw_lncRNA_counts_matrix.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)



####Agilent####
rm(list = ls())
library(stringr)
library(AnnoProbe)
library(GEOquery)
library(limma)
gse = "GSE213499"
geoChina(gse)#这个不行就用下面的
getGEO("GSE213499", destdir=".", AnnotGPL = F, getGPL = F)
f='GSE213499_eSet.Rdata'
if(!file.exists(f)){
  gset <- getGEO('GSE213499', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}

load("GSE213499_eSet.Rdata")
pd <- pData(gset[[1]])

raw_dir = "GSE213499_RAW/"#这个是原始数据文件的名字加/
raw_datas = paste0(raw_dir,"/",dir(raw_dir))

#调整pd与rawdata的顺序一致
raw_order = str_extract(raw_datas,"GSM\\d*")
pd = pd[match(raw_order,rownames(pd)),]

table(pd$source_name_ch1)#注意是不是source_name_ch1
group_list <- ifelse(str_detect(pd$source_name_ch1, "non-gastric"), "normal","tumor")

table(pd$characteristics_ch1)#注意是不是source_name_ch1
group_list <- ifelse(str_detect(pd$characteristics_ch1, "cancer"), "tumor",
                     "normal")


group_list <- factor(group_list, levels=c("normal","tumor"))

x <- read.maimages(raw_datas,
                   source="agilent", 
                   green.only=TRUE,
                   other.columns = "gIsWellAboveBG")
dim(x)
y <- backgroundCorrect(x, method="normexp")
y <- normalizeBetweenArrays(y, method="quantile")
class(y)

Control <- y$genes$ControlType==1L;table(Control)
## Control
## FALSE  TRUE 
## 43529  1486
NoSymbol <- is.na(y$genes$ProbeName);table(NoSymbol)
## NoSymbol
## FALSE 
## 45015
IsExpr <- rowSums(y$other$gIsWellAboveBG > 0) >= 16;table(IsExpr)
## IsExpr
## FALSE  TRUE 
## 13088 31927
Isdup <- duplicated(y$genes$ProbeName);table(Isdup)
## Isdup
## FALSE  TRUE 
## 30328 14687
yfilt <- y[!Control & !NoSymbol & !Isdup, ]#见机调整
yfilt <- y[!Control & !NoSymbol & IsExpr & !Isdup, ]
dim(yfilt)
exp = yfilt@.Data[[1]]
boxplot(exp)
exp[1:2,1:2]
#获取样本名
colnames(exp) = str_extract(colnames(exp),"GSM\\d*")
exp[1:2,1:2]
#获取基因名
anno = yfilt$genes
nrow(anno);nrow(exp)
## [1] 20650
## [1] 20650
rownames(exp)=rownames(anno)
ids = unique(anno$ProbeName)
exp = exp[!duplicated(anno$ProbeName),]
rownames(exp) = ids
exp[1:4,1:4]
write.table(exp, file = "Agilent表达矩阵.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)

#差异分析
design <- model.matrix(~group_list)
fit <- lmFit(exp,design)
fit <- eBayes(fit,trend=TRUE,robust=TRUE)
summary(decideTests(fit))
##        (Intercept) group_listtumor
## Down             0            2102
## NotSig           0           16928
## Up           20650            1620
deg = topTable(fit,coef=2,n=dim(y)[1])
boxplot(exp[rownames(deg)[1],]~group_list)



####affymetrix####
library(affy)
#perform mas5 normalization
affy_data = ReadAffy(celfile.path=dir_cels)
eset.mas5 = mas5(affy_data)
exprSet.nologs = exprs(eset.mas5)
exprSet = log(exprSet.nologs, 2)  #transform to Log_2 if needed

library(affy)
data <- ReadAffy(celfile.path=dir_cels) 
eset <- rma(data)
write.exprs(eset,file="data.txt")



#新版
studyID='GSE42872';
studyID_probe=paste0(studyID,'_probe')
studyID_gene=paste0(studyID,'_gene')
R_history_data <- paste0(studyID,'.Rdata')
if ( file.exists(R_history_data)){
  load( R_history_data )
}else{
  library(oligo)
  ## the directory data should have raw CEL files 
  geneCELs=list.celfiles('data',listGzipped=T,full.name=T) 
  affyGeneFS <- read.celfiles(geneCELs)   
  geneCore <- rma(affyGeneFS, target = "core") 
  genePS <- rma(affyGeneFS, target = "probeset") 
  featureData(genePS) <- getNetAffx(genePS, "probeset")
  featureData(geneCore) <- getNetAffx(geneCore, "transcript")
  
}

exprSet=exprs( genePS )
pdata=pData( genePS )








####ID转换大盘点####
#第一种：bitr
#涵盖种类包括ACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT, ENSEMBLTRANS, 
#ENTREZID, ENZYME, EVIDENCE, EVIDENCEALL, GENENAME, GENETYPE, 
#GO, GOALL, IPI, MAP, OMIM, ONTOLOGY, ONTOLOGYALL, PATH, PFAM, 
#PMID, PROSITE, REFSEQ, SYMBOL, UCSCKG, UNIPROT.
library(clusterProfiler)
library(org.Hs.eg.db)
x = bitr(A,fromType = "SYMBOL",#C为文件名，从SYMBOL转换为ENSEMBL
         toType = "ENSEMBL",OrgDb = "org.Hs.eg.db")
head(x)

#第二种：mygene
xli <-  c('DDX26B','CCDC83',  'MAST3',
          'RPL11', 'ZDHHC20',  'LUC7L3',  
          'SNORD49A',  'CTSH', 'ACOT8')
tmp = queryMany(xli, scopes="symbol", 
                fields=c("uniprot", "ensembl.gene", "reporter"), 
                species="human")
colnames(tmp)

#第三种：gtf
library(rtracklayer)
g = import("gencode.v22.annotation.gtf.gz")
colnames(as.data.frame(g))
g = as.data.frame(g)
g = g[,c("gene_name","gene_id","gene_type")]
head(g)
b = g[g$gene_name %in% xli,]
head(b)

#第四种：从GEO，soft文件提取
rm(list = ls())
options(stringsAsFactors = F)
library(GEOquery)
gse <- getGEO(filename = "GSE50710_family.soft.gz",destdir = ".")#手动下
str(gse)
length(gse)
id_probe <- gse@gpls$GPL13825@dataTable@table#记得改GPL13825
dim(id_probe)
head(id_probe)
colnames(id_probe) <- id_probe$ID#将ID变为行名
id_probe$Gene <- as.character(id_probe$ID)   #新增Gene Symbol
id_probe <- id_probe[,-1] 
#手动导入基因
A <- `上调基因`
B <- id_probe
comgene <- intersect(rownames(A),rownames(B))
A <- A[comgene,]
B <- B[comgene,]
C <- cbind(B,A)
write.table(C, file = "上调基因Ensembl2.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)



#第五种：R包转换id
index = gset[[1]]@annotation
#if(!require("hgu133plus2.db"))
#  BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
ls("package:hgu133plus2.db")
ids <- toTable(hgu133plus2SYMBOL)
head(ids)
#length(unique(ids$symbol))
#table(sort(table(ids$symbol)))
#id转换
library(tidyverse)
exp <-A#上调基因
exp <- exp %>% mutate(probe_id=rownames(exp))
exp <- exp %>% inner_join(ids,by="probe_id") 
exp <- exp[!duplicated(exp$symbol),]
rownames(exp) <- exp$symbol
write.table(exp, file = "上调基因Ensembl2.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)




#高通量测序的表达矩阵#不可用
rm(list = ls())
options(stringsAsFactors = F)#在调用as.data.frame的时，将stringsAsFactors设置为FALSE可以避免character类型自动转化为factor类型
# 注意查看下载文件的大小，检查数据 
f='GSE174237_eSet.Rdata'
library(GEOquery)
if(!file.exists(f)){
  gset <- getGEO('GSE103611', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE174237_eSet.Rdata')  ## 载入数据
class(gset)  #查看数据类型
length(gset)  #
class(gset[[1]])
gset




#ensembl_id转换方法2
rm(list = ls())  ## 魔幻操作，一键清空~
#手动载入文件
A <- `上调基因Ensembl`
library(stringr) 
ids=data.frame(Ensembl.Gene=str_split(rownames(A),
                                      '[.]',simplify = T)[,1],
               median=A$AveExpr
)
head(ids)
head(A$Ensembl.Gene)

library(org.Hs.eg.db)
g2s=unique(toTable(org.Hs.egSYMBOL))
head(g2s)
g2e=unique(toTable(org.Hs.egENSEMBL)) 
head(g2e)
s2e=merge(g2e,g2s,by='gene_id')
head(s2e)
table(ids$ensembl_id %in% s2e$ensembl)

ids=ids[ids$ensembl_id %in% s2e$ensembl,]
#取出在gencode数据库的gtf注释中能找到的ensembl_id
ids$symbol=s2e[match(ids$ensembl_id,s2e$ensembl),2]
#match返回其第二个参数中第一个参数匹配的位置
# 把s2e的ensembl按照ids$ensembl的顺序一个个取出来，从而得到ids$symbol这一列
length(unique(ids$symbol))
head(ids) 
ids=ids[order(ids$symbol,ids$median,decreasing = T),]#把ids$symbol按照ids$median排序
ids=ids[!duplicated(ids$symbol),]#取出不重复的ids$symbol
dim(ids) 
exprSet= exprSet[rownames(ids),]#取出表达矩阵中ids有的行
rownames(exprSet)=ids$symbol#把ids$symbol变为exprSet的行名
exprSet[1:4,1:4]  
dim(exprSet)

#ensembl_id转换方法3

#手动载入文件
A <- `上调基因Ensembl`
library(stringr) 
ids=data.frame(Ensembl.Gene=str_split(rownames(A),
                                      '[.]',simplify = T)[,1],
               median=A$AveExpr
)
head(ids)
head(A$Ensembl.Gene)
load('human_geneInfo_genecode_v25.rda')#gencode数据库的gtf注释
head(human_geneInfo_genecode_v25)#可以看到有symbol和ensembl的对应关系








