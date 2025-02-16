####深度学习汇总####
library(VIM)
library(mice)
rm(list = ls())
####—数据预处理—####
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE,row.names = 1)
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
data<-na.omit(data)#删除包含空值的行
table(is.na(data))#md.pattern(data)
aggr(data,prop=FALSE,numbers=TRUE)#判断数据是否有缺失值
data=data[complete.cases(data),]#删除缺失值
#标准化
#scale函数:center和scale默认为真，如果只有scale为真则 ：实际值/SD
data_scale<-scale(data[5:682],center=T,scale=T)#从第二列到第2819列#（公式：(实际值-均值)/sd）
data_scale<- as.data.frame(data_scale)#变表格
#归一化
f1<-function(x){return((x-min(x)) / (max(x)-min(x)))}
data_scale2<-as.data.frame(lapply(data[5:682],f1))
data_scale3<-as.data.frame(lapply(data_scale[1:678],f1))
#data_scale仅标准化
#data_scale2仅归一化
#data_scale3标准化+归一化
aggr(data_scale3,prop=FALSE,numbers=TRUE)#判断数据是否有缺失值
write.table(data_scale3, file = "标准化后数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)



####Spearman检验####
rm(list = ls())
data<-data_scale3
data<-na.omit(data)#删除包含空值的行
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
rt <- data
norm_result <- apply(rt, 2, function(x) shapiro.test(x)$p.value)
norm_feature <- rt[which(norm_result >= 0.05)]
cor_nor <- cor(norm_feature, method = "pearson")
cor_all <- cor(rt, method = "spearman")
num_nor <- dim(cor_nor)[1]
cor_all[1:num_nor, 1:num_nor] <- cor_nor
cor_all[upper.tri(cor_all)] <- 0
diag(cor_all) <- 0
data_reduce = rt[, !apply(cor_all, 1, function(x) any(abs(x) > 0.90))]
dim(data_reduce)
write.table(data_reduce, file = "Spearman后数据瘤周.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)



####批量相关性分析####
library(tidyverse)
library(stringr)
data<- read.csv("D:/Desktop/R.csv",header = TRUE)
m1=data[,c(1:10)]#影像组学参数
m2=data[,c(11:20)]#基因
cor_2_matrix <- function(m1,m2){
  apply(m2 , 2, function(x){
    #  x = m2[2,]
    unlist(apply(m1, 2,function(y){
      # y  = m1[1,]
      cor(as.numeric(x),
          as.numeric(y))}))})}
rdf=cor_2_matrix(m1, m2) %>% as.data.frame()
rdf$gene=rownames(rdf)
rdf <- rdf %>% gather(key = 'RNA',value = 'r',-gene)

corP_2_matrix <- function(m1,m2){
  apply(m2 , 2, function(x){
    #  x = m2[2,]
    unlist(apply(m1, 2,function(y){
      # y  = m1[1,]
      cor.test(as.numeric(x),
               as.numeric(y))$p.value}))})}
pdf=corP_2_matrix(m1, m2) %>% as.data.frame()
pdf$gene=rownames(pdf) 
pdf <- pdf %>% gather(key = 'RNA',value = 'p',-gene)

Toal <- bind_cols(pdf,rdf) %>%.[-c(4,5)]
Toal2 <- filter(Toal,abs(r)>0.6&p<0.001)
colnames(Toal2) <- c('Radiomics','RNA','pvalue','corr')

write.table(Toal2, file = "组学与基因的相关性分析.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)



####相关性热图####
library(corrplot)
library(ggplot2)
library(ggcorrplot)
library(vcd)
library(psych)
library(ggrepel)
#数据导入#
data<- read.csv("D:/Desktop/R.csv",header = TRUE)
data<-as.matrix(data) #利用as.matrix()将所需数据集转换为matrix格式，才可在corrplot中跑
data=data.frame(scale(data))#数据标准化
head(data)
#相关性计算
data<-cor(data,method="spearman") #pearson，spearman和kendall
round(data, 2)#保留两位小数

#相关性热图绘制#
ggcorrplot(data, method="circle") #圆圈大小变化

#带数值
ggcorrplot(data, method = "circle", #"square", "circle"相比corrplot少了很多种，只有方形和圆形，默认方形。
           type ="upper" , #full完全(默认)，lower下三角，upper上三角
           ggtheme = ggplot2::theme_minimal,
           title = "",
           show.legend = TRUE,  #是否显示图例。
           legend.title = "Corr", #指定图例标题。
           show.diag =T ,    #FALSE显示中间
           colors = c("blue", "white", "red"), #需要长度为3的颜色向量，同时指定low,mid和high处的颜色。
           outline.color = "gray", #指定方形或圆形的边线颜色。
           hc.order = FALSE,  #是否按hclust(层次聚类顺序)排列。
           hc.method = "complete", #相当于corrplot中的hclust.method, 指定方法一样，详情见?hclust。
           lab =T , #是否添加相关系数。FALSE
           lab_col = "black", #指定相关系数的颜色，只有当lab=TRUE时有效。
           lab_size = 4, #指定相关系数大小，只有当lab=TRUE时有效。
           p.mat = NULL,  #p.mat= p_mat,insig= "pch", pch.col= "red", pch.cex= 4,
           sig.level = 0.05,
           insig = c("pch", "blank"),
           tl.cex = 12, #指定变量文本的大小，
           tl.col = "black", #指定变量文本的颜色，
           tl.srt = 45, #指定变量文本的旋转角度。
           digits = 2 #指定相关系数的显示小数位数(默认2)。
)

#corrplot包绘图#
corrplot(data)
corrplot(data, method="circle", #square方形，ellipse, 椭圆形，number数值，shade阴影，color颜色，pie饼图
         #title = "pearson",   #指定标题
         type="full",  #full完全(默认)，lower下三角，upper上三角
         #col=c("#FF6666", "white", "#0066CC"), #指定图形展示的颜色，默认以均匀的颜色展示。支持grDevices包中的调色板，也支持RColorBrewer包中调色板。
         outline = T,  #是否添加圆形、方形或椭圆形的外边框，默认为FALSE。
         diag = TRUE,  #是否展示对角线上的结果，默认为TRUE
         mar = c(0,0,0,0), #设置图形的四边间距。数字分别对应(bottom, left, top, right)。
         bg="white", #指定背景颜色
         add = FALSE, #表示是否添加到已经存在的plot中。默认FALSE生成新plot。
         is.corr = TRUE, #是否为相关系数绘图，默认为TRUE,FALSE则可将其它数字矩阵进行可视化。
         addgrid.col = "darkgray", #设置网格线颜色，当指定method参数为color或shade时， 默认的网格线颜色为白色，其它method则默认为灰色，也可以自定义颜色。
         addCoef.col = NULL, #设置相关系数值的颜色，只有当method不是number时才有效
         addCoefasPercent = FALSE, #是否将相关系数转化为百分比形式，以节省空间，默认为FALSE。
         order = "original", #指定相关系数排序的方法, 可以是original原始顺序，AOE特征向量角序，FPC第一主成分顺序，hclust层次聚类顺序，alphabet字母顺序。
         hclust.method = "complete", # 指定hclust中细分的方法，只有当指定order参数为hclust时有效。有7种可选：complete,ward,single,average,mcquitty,median,centroid。
         addrect = NULL, #是否添加矩形框，只有当指定order参数为hclust时有效， 默认不添加， 用整数指定即可添加。
         rect.col = "black", #指定矩形框的颜色。
         rect.lwd = 2, #指定矩形框的线宽。
         tl.pos = NULL,  #为"lt"(左侧和顶部),"ld"(左侧和对角线),"td"(顶部和对角线),"d"(对角线),"n"(无);当type="full"时默认"lt"。当type="lower"时默认"ld"。
         tl.cex = 1,  #设置文本标签的大小
         tl.col = "black", #设置文本标签的颜色。
         cl.pos = NULL #设置图例位置，为"r"(右边),"b"(底部),"n"(无)之一。当type="full"/"upper"时，默认"r"; 当type="lower"时，默认"b"。
         #addshade = c("negative", "positive", "all"), # 表示给增加阴影，只有当method="shade"时有效"all"(对所有相关系数增加阴影)。
         #shade.lwd = 1,  #指定阴影线宽。
         #shade.col = "white",  #指定阴影线的颜色。
         #p.mat= res1$p,sig.level= 0.01,insig= "pch", pch.col= "blue", pch.cex= 3,#只有指定矩阵的P值，sig.level，pch等参数才有效。。
)

#显示数字与图形混合
corrplot(data, method="circle", #square方形，ellipse, 椭圆形，number数值，shade阴影，color颜色，pie饼图
         title = "pearson",   #指定标题
         type="full", #full完全(默认)，lower下三角，upper上三角
         #col=c("#FF6666", "white", "#0066CC"), #指定图形展示的颜色，默认以均匀的颜色展示。支持grDevices包中的调色板，也支持RColorBrewer包中调色板。
         outline = F,  #是否添加圆形、方形或椭圆形的外边框，默认为FALSE。
         diag = TRUE,  #是否展示对角线上的结果，默认为TRUE
         mar = c(0,0,0,0), #设置图形的四边间距。数字分别对应(bottom, left, top, right)。
         bg="white", #指定背景颜色
         add = FALSE, #表示是否添加到已经存在的plot中。默认FALSE生成新plot。
         is.corr = TRUE, #是否为相关系数绘图，默认为TRUE,FALSE则可将其它数字矩阵进行可视化。
         addgrid.col = "darkgray", #设置网格线颜色，当指定method参数为color或shade时， 默认的网格线颜色为白色，其它method则默认为灰色，也可以自定义颜色。
         addCoef.col = NULL, #设置相关系数值的颜色，只有当method不是number时才有效
         addCoefasPercent = FALSE, #是否将相关系数转化为百分比形式，以节省空间，默认为FALSE。
         order = "original", #指定相关系数排序的方法, 可以是original原始顺序，AOE特征向量角序，FPC第一主成分顺序，hclust层次聚类顺序，alphabet字母顺序。
         hclust.method = "complete", # 指定hclust中细分的方法，只有当指定order参数为hclust时有效。有7种可选：complete,ward,single,average,mcquitty,median,centroid。
         addrect = NULL, #是否添加矩形框，只有当指定order参数为hclust时有效， 默认不添加， 用整数指定即可添加。
         rect.col = "black", #指定矩形框的颜色。
         rect.lwd = 2, #指定矩形框的线宽。
         tl.pos = NULL,  #指定文本标签(变量名称)相对绘图区域的位置，为"lt"(左侧和顶部),"ld"(左侧和对角线),"td"(顶部和对角线),"d"(对角线),"n"(无)。
         tl.cex = 1,  #设置文本标签的大小
         tl.col = "black", #设置文本标签的颜色。
         cl.pos = NULL #设置图例位置，为"r"(右边),"b"(底部),"n"(无)之一。当type="full"/"upper"时，默认"r"; 当type="lower"时，默认"b"。
         #addshade = c("negative", "positive", "all"), ##为"negative"(对负相关系数增加阴影135度)；"positive"(对正相关系数增加阴影45度)。
         #shade.lwd = 1,  #指定阴影线宽。
         #shade.col = "white",  #指定阴影线的颜色。
         #p.mat= res1$p,sig.level= 0.01,insig= "pch", pch.col= "blue", pch.cex= 3只有当insig = "pch"时，pch.col和pch.cex参数才有效。
)

corrplot(data, title = "",        
         method = "number", #square方形，ellipse, 椭圆形，number数值，shade阴影，color颜色，pie饼图       
         outline = F, #是否添加圆形、方形或椭圆形的外边框，默认为FALSE。
         add = TRUE, #表示是否添加到已经存在的plot中。默认FALSE生成新plot。
         type = "full", #full完全(默认)，lower下三角，upper上三角       
         order="original",
         col="black", #指定图形展示的颜色，默认以均匀的颜色展示。支持grDevices包中的调色板，也支持RColorBrewer包中调色板。
         diag=FALSE, #是否展示对角线上的结果，默认为TRUE
         tl.pos="n",  #指定文本标签(变量名称)相对绘图区域的位置，为"lt"(左侧和顶部),"ld"(左侧和对角线),"td"(顶部和对角线),"d"(对角线),"n"(无)
         cl.pos=NULL #设置图例位置，为"r"(右边),"b"(底部),"n"(无)之一。
)


#椭圆加数值#
corrplot(data, method = "ellipse", order = "original",         
         addCoef.col = "black",#设置相关系数值的颜色，只有当method不是number时才有效
         type="full", #full完全(默认)，lower下三角，upper上三角
         title = "椭圆与黑色系数值",
         add = FALSE, #表示是否添加到已经存在的plot中。默认FALSE生成新plot。
         diag = TRUE, #是否展示对角线上的结果，默认为TRUE
         tl.cex = 1,  #设置文本标签的大小
         tl.col = "black", #设置文本标签的颜色。
         cl.pos = NULL, #设置图例位置，为"r"(右边),"b"(底部),"n"(无)之一。当type="full"/"upper"时，默认"r"; 当type="lower"时，默认"b"。
         mar = c(1,1,1,1)) #设置图形的四边间距。数字分别对应(bottom, left, top, right)。

#百分比表示#
corrplot(data, method = "ellipse", order = "original",         
         addCoef.col = "black",#设置相关系数值的颜色，只有当method不是number时才有效
         addCoefasPercent = TRUE, #是否将相关系数转化为百分比形式，以节省空间，默认为FALSE。
         type="full", #full完全(默认)，lower下三角，upper上三角
         title = "椭圆与黑色百分比",
         add = FALSE, #表示是否添加到已经存在的plot中。默认FALSE生成新plot。
         diag = TRUE, #是否展示对角线上的结果，默认为TRUE
         tl.cex = 1,  #设置文本标签的大小
         tl.col = "black", #设置文本标签的颜色。
         cl.pos = NULL, #设置图例位置，为"r"(右边),"b"(底部),"n"(无)之一。当type="full"/"upper"时，默认"r"; 当type="lower"时，默认"b"。
         mar = c(1,1,1,1)) #设置图形的四边间距。数字分别对应(bottom, left, top, right)。




####批量单因素COX####
rm(list=ls())
set.seed(123)#
data<- read.csv("D:/Desktop/R.csv",header = TRUE)
credit_df <- data
index<-sample(1:nrow(credit_df),nrow(credit_df)*0.67,replace=F)#0.7代表37分
training_dataset<-credit_df[index,]
test_dataset<-credit_df[-index,]
dim(training_dataset)#训练集training_dataset
dim(test_dataset)#验证集test_dataset

R<- training_dataset
can <- colnames(R)[5:53]#批量计算那些列
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
write.table(result, file = "批量单因素COX肿瘤.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)



#多个ROC####
library(pROC)#不同数据集
ROC1 <- roc(a$label,a$"0", ci = T)
ROC2 <- roc(b$label,b$"0", ci = T)
plot(smooth(ROC1), col = "red", lwd = 2)
plot(smooth(ROC2), col = "blue", lwd = 2, add = TRUE)
ROC1
ROC2
#两条线
legend("bottomright", cex = 1.1, title="NaiveBayes",
       legend = c("Training Set (AUC: 0.901)", "Validation Set (AUC: 0.603)"),
       col = c("red", "blue"), lty = 1, lwd = 2)  #定义标签
#三条线
plot(smooth(ROC3), col = "darkgreen", lwd = 2, add = TRUE)
legend("bottomright", cex = 1.1,
       legend = c("s100b (AUC: 0.731)", "ndka (AUC: 0.612)", "wfns (AUC: 0.824)"),
       col = c("red", "blue", "darkgreen"), lty = 1, lwd = 2)  #定义标签



####——二分类——####
###1分类Lasso####
library(Matrix)
library(foreach)
library(glmnet)
library(foreign)
library(pROC)
library(ggplot2)
#二分类
set.seed(120)  
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
credit_df <- data
index<-sample(1:nrow(credit_df),nrow(credit_df)*0.7,replace=F)#0.7代表37分
train_data<-credit_df[index,]
test_data<-credit_df[-index,]
dim(train_data)
dim(test_data)

x <- as.matrix(train_data[,c(2:108)])
y <- as.matrix(train_data[,1])
#glmnet()中alpha=0表示岭回归，alpha=1表示lasso回归
lasso <- glmnet(x,y,alpha=1,family = "binomial",nlambda = 100)
#print(lasso)
plot(lasso, xvar="lambda", label=TRUE,lwd=1.4)#,lty=1
#predict(lasso, newx=x[2:5,], type = "response")
cvfit=cv.glmnet(x,y)
plot(cvfit)
lambda.min<-cvfit$lambda.min
lambda.1se<-cvfit$lambda.1se
coef_best<-coef(cvfit,s ="lambda.min",exact = F)#lamda.min精准，1se简洁
#coef_best
index<-which(coef_best!=0)#非零系数
index=data.frame(index)
index=index[-1,]
coef<-coef_best[index]#对应回归系数
diffvariables=row.names(coef_best)[index]#非零变量
sames=intersect(diffvariables,colnames(data))

#训练集
lassoout0=train_data[ ,sames,drop=F]
lassoout=data.frame(label=train_data$label,lassoout0)
myFun=function(x){crossprod(as.numeric(x),coef)}
lassoScore=apply(lassoout0,1,myFun)
outCol1=c("label",diffvariables)#这句代码就是提取data里的表头内容
outCol1<- as.data.frame(outCol1)#变表格
outCol1 <- outCol1[-2,]   #去除第2行
outTab1=cbind(train_data[,outCol1],riskScore=as.vector(lassoScore))
write.table(outTab1, file = "安医lasso纳入参数1.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
roc = roc(outTab$label,outTab$riskScore) #注意一次只能一列
roc_auc <- auc(roc)
roc_auc
plot(roc, col="red", lwd=3, title = "")  



#验证集
lassoout1=test_data[ ,sames,drop=F]
lassoout2=data.frame(label=test_data$label,lassoout1)
lassoScore=apply(lassoout1,1,myFun)
outCol2=c("label",diffvariables)
outCol2<- as.data.frame(outCol2)#变表格
outCol2 <- outCol2[-2,]   #去除第2行
outTab2=cbind(test_data[,outCol2],riskScore=as.vector(lassoScore))
roc = roc(outTab2$label,outTab2$riskScore) #注意一次只能一列
roc_auc <- auc(roc)
roc_auc



####2决策树(训练+验证)####
rm(list=ls())
library(pROC)
library(rpart.plot)
library(treeheatr)
library(caret)
library(rpart)
library(tidyverse)
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
#设置随机种子，保证结果可复现
set.seed(100)  
credit_df <- data
index<-sample(1:nrow(credit_df),nrow(credit_df)*0.5,replace=F)#0.7代表37分
train_data<-credit_df[index,]
test_data<-credit_df[-index,]
dim(train_data)
dim(test_data)
#训练集
tree = rpart(label~.,data =train_data,method = 'class')
Yhat = predict(tree,train_data,type = 'prob') #预测概率值  
a <- as.data.frame(Yhat)#变表格
a$label<- as.character(train_data$label) #新增一列label
a <- a[,-2] #去除第2列
modelroc1 = roc(a$label,a$"0") #注意一次只能一列
roc_auc <- auc(modelroc1)
roc_auc

#决策树画图
par(xpd = NA)
plot(tree)
text(tree)
predicted<- tree %>% predict(train_data, type = "class")
head(predicted)
#准确率
mean(predicted == train_data$label)

#保存一下模型的原始数据，用于后续的分析
write.table(a, file = "训练集决策树ROC数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(train_data, file = "训练集决策树原始数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)

#验证集
Yhat = predict(tree,test_data,type = 'prob') #预测概率值  
#modelroc2 = roc(test_data$label,Yhat[,1])
b<- as.data.frame(Yhat)#变表格
b$label<- as.character(test_data$label) #新增一列label
b <- b[,-2] #去除第2列
modelroc2 = roc(b$label,b$"0") #注意一次只能一列
#plot(modelroc2)
roc_auc <- auc(modelroc2)
roc_auc

#保存一下模型的原始数据，用于后续的分析
write.table(b, file = "验证集决策树ROC数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(test_data, file = "验证集决策树原始数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)






####3随机森林####
rm(list=ls())
library(dplyr) #数据处理使用
library(data.table) #数据读取使用
library(randomForest) #RF模型使用
library(caret) # 调参和计算模型评价参数使用
library(pROC) #绘图使用
library(ggplot2) #绘图使用
library(ggpubr) #绘图使用
library(ggprism)

credit_df <- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
set.seed(123)  #设置随机种子，保证结果可复现
index<-sample(1:nrow(credit_df),nrow(credit_df)*0.7,replace=F)#0.7代表37分
train_data<-credit_df[index,]
test_data<-credit_df[-index,]
dim(train_data)
dim(test_data)
# 定义训练集特征和目标变量
X_train <- train_data[, -1]
y_train <- as.factor(train_data[, 1])
# 创建随机森林分类模型
model<- randomForest(x=X_train, y=y_train)
# 输出默认参数下的模型性能
print(model)


####进行参数调优
ctrl <- trainControl(method = "cv", number = 10) #使用五折交叉验证，也可以选择10折交叉验证。
grid <- expand.grid(mtry = c(2, 4, 10))#定义参数网格,每棵树中用于分裂的特征数量，这里只是随便给的测试，主要为了介绍如何调参，并非最优选择。
rf_model <- train(x = X_train, y = y_train,
                  method = "rf",
                  trControl = ctrl,
                  tuneGrid = grid)#使用caret包进行调参

#输出最佳模型和参数
print(rf_model)#注意The final value used for the model was mtry = ？
#手动调参。mtry = c(?)
grid <- expand.grid(mtry = c(4))#要根据上面mtry改，代表每棵树中用于分裂的特征数量
modellist <- list()# 定义模型列表，存储每一个模型评估结果
# 调整的参数是决策树的数量
for (ntree in c(100,200,300)) {
  set.seed(123)
  fit <- train(x = X_train, y = y_train, method="rf", 
               metric="Accuracy", tuneGrid=grid, 
               trControl=ctrl, ntree=ntree)
  key <- toString(ntree)
  modellist[[key]] <- fit
}
results <- resamples(modellist)#compare results
#输出最佳模型和参数
summary(results)#根据准确性确定ntree是100,200,300
#使用最佳参数训练最终模型
final_model <- randomForest(x=X_train,y=y_train,mtry = 2,ntree = 300)
# 输出最终模型
print(final_model)



#####在训练集上进行预测
train_predictions <- predict(model, data=train_data)#用那个model自己看
#用哪个模型
#train_predictions <- predict(final_model, data = train_data)#用那个model自己看
#计算模型指标
confusion_matrix <- confusionMatrix(train_predictions, y_train)
accuracy <- confusion_matrix$overall["Accuracy"]
precision <- confusion_matrix$byClass["Pos Pred Value"]
recall <- confusion_matrix$byClass["Sensitivity"]
f1_score <- confusion_matrix$byClass["F1"]
# 输出模型指标
print(confusion_matrix)
print(paste("Accuracy:", accuracy))
print(paste("Precision:", precision))
print(paste("Recall:", recall))
print(paste("F1 Score:", f1_score))
# 绘制混淆矩阵热图
# 将混淆矩阵转换为数据框
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("cluster1","cluster2")
rownames(confusion_matrix_df) <- c("cluster1","cluster2")
draw_data <- round(confusion_matrix_df / rowSums(confusion_matrix_df),2)
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)

# 绘制ROC曲线需要将预测结果以概率的形式输出
train_predictions <- predict(model, data = train_data,type = "prob")
# 计算ROC曲线的参数
roc_obj1 <- roc(response = y_train, predictor = train_predictions[, 2])
roc_auc <- auc(roc_obj1)
roc_auc 

# 将ROC对象转换为数据框
roc_data <- data.frame(1 - roc_obj1$specificities, roc_obj1$sensitivities)

#保存一下模型的原始数据，用于后续的分析
a <- as.data.frame(train_predictions)#变表格
a$label<- as.character(train_data$label) #新增一列label
a <- a[,-2] #去除第2列
write.table(a, file = "汇总训练集随机森林ROC数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


#####在测试集上进行预测
X_test <- test_data[, -1]
y_test <- as.factor(test_data[, 1])
test_predictions <- predict(model, newdata = test_data)#final_model
# 计算模型指标
confusion_matrix <- confusionMatrix(test_predictions, y_test)
accuracy <- confusion_matrix$overall["Accuracy"]
precision <- confusion_matrix$byClass["Pos Pred Value"]
recall <- confusion_matrix$byClass["Sensitivity"]
f1_score <- confusion_matrix$byClass["F1"]
# 输出模型指标
print(confusion_matrix)
print(paste("Accuracy:", accuracy))
print(paste("Precision:", precision))
print(paste("Recall:", recall))
print(paste("F1 Score:", f1_score))
# 绘制混淆矩阵热图
# 将混淆矩阵转换为数据框
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("cluster1","cluster2")
rownames(confusion_matrix_df) <- c("cluster1","cluster2")
draw_data <- round(confusion_matrix_df / rowSums(confusion_matrix_df),2)
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)

# 绘制ROC曲线需要将预测结果以概率的形式输出
test_predictions <- predict(model, newdata = test_data,type = "prob")
# 计算ROC曲线的参数
roc_obj2 <- roc(response = y_test, predictor = test_predictions[, 2])
roc_obj2
# 将ROC对象转换为数据框
roc_data <- data.frame(1 - roc_obj2$specificities, roc_obj2$sensitivities)


#保存一下模型的原始数据，用于后续的分析
b <- as.data.frame(test_predictions)#变表格
b$label<- as.character(test_data$label) #新增一列label
b <- b[,-2] #去除第2列
write.table(b, file = "汇总验证集随机森林ROC数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)

#画图
plot(roc_obj1, col = "red", lwd = 2)
plot(roc_obj2, col = "blue", lwd = 2, add = TRUE)
roc_obj1
roc_obj2
legend("bottomright", cex = 1.1, title="随机森林",
       legend = c("训练集 (AUC: 0.780)", "验证集 (AUC: 0.781)"),
       col = c("red", "blue"), lty = 1, lwd = 2)  #定义标签


#####特征排序####
DEG <-as.data.frame(model$importance)
DEG$Gene <- rownames(DEG)
DEG <- DEG[order(-DEG$MeanDecreaseGini), ]#对差异基因的p值进行从大到小的排序
plot.df <- DEG[1:10,]

ggplot(aes(x=reorder(Gene,MeanDecreaseGini),y=MeanDecreaseGini,fill=Gene),data=plot.df)+
  geom_col()+
  coord_flip()+
  theme_bw()+
  labs(x="")+
  ggtitle("Feature sorting in RF")+
  scale_fill_brewer(palette = "Set3")+
  theme(legend.position ="")



####4朴素贝叶斯####
rm(list=ls())
library(e1071)
library(pROC)
library(tidyverse)
library(naivebayes)
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
#设置随机种子，保证结果可复现19
set.seed(123)  
credit_df <- data
index<-sample(1:nrow(credit_df),nrow(credit_df)*0.7,replace=F)#0.7代表37分
train_data<-credit_df[index,]
test_data<-credit_df[-index,]
dim(train_data)
dim(test_data)
#训练集
nb1 <- naiveBayes(label~., data =train_data )#计算NB
Yhat <- predict(nb1,train_data,type = 'raw')
a <- as.data.frame(Yhat)#变表格
a$label<- as.character(train_data$label) #新增一列label
rownames(a) <- rownames(train_data)   #将行名变为Gene Symbol
modelroc1 = roc(a$label,a$"0") #注意一次只能一列
#plot(modelroc1)
roc_auc <- auc(modelroc1)
roc_auc
#ggplot画风的ROC 
(g=ggroc(list(modelroc1),legacy.axes = T,lwd=1.5))
#g+theme_minimal() + ggtitle("My ROC curve") + geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed")
#保存一下模型的原始数据，用于后续的分析
write.table(a, file = "训练集朴素贝叶斯ROC数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(train_data, file = "训练集朴素贝叶斯原始数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)



#验证集
Yhat <- predict(nb1,test_data,type = 'raw')
b<- as.data.frame(Yhat)#变表格
b$label<- as.character(test_data$label) #新增一列label
b <- b[,-2] #去除第2列
modelroc2 = roc(b$label,b$"0") #注意一次只能一列
#plot(modelroc2)
roc_auc <- auc(modelroc2)
roc_auc

#保存一下模型的原始数据，用于后续的分析
write.table(b, file = "验证集朴素贝叶斯ROC数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(test_data, file = "验证集朴素贝叶斯原始数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)





####5支持向量机####
#线性核，多项式核，高斯径向基核，sigmoid核
rm(list = ls())
library(modeldata)
library(tidyverse)
library(pROC)
library(e1071)
library(ggplot2)
library(GGally)
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
credit_df<-na.omit(data)
#经典的37分
set.seed(124)
index<-sample(1:nrow(credit_df),nrow(credit_df)*0.7,replace=F)#0.7代表37分
train_df<-credit_df[index,]
test_df<-credit_df[-index,]
dim(train_df)
dim(test_df)


svmLinear<-svm(label~.,data=train_df,
               probability=TRUE,
               kernel="linear")

svmPoly<-svm(label~.,data=train_df,
             probability=TRUE,
             kernel="polynomial")

svmRadial<-svm(label~.,data=train_df,
               probability=TRUE,
               kernel="radial")

svmSigmoid<-svm(label~.,data=train_df,
                probability=TRUE,
                kernel="sigmoid")
#定义函数
getres<-function(svmfunc,dataset){
  data_pred<-predict(svmfunc,newdata=dataset,probability=T,type="class")
  data_pred_df<-dataset%>%select(label)%>%
    bind_cols(status_pred=data_pred)%>%
    bind_cols(attr(data_pred,"probabilities"))}

Linear_train_pred_df<-getres(svmLinear,train_df)
head(Linear_train_pred_df)


#提取4种核函数分别在训练集、测试集的结果
Linear_test_pred_df<-getres(svmLinear,test_df)
Poly_train_pred_df<-getres(svmPoly,train_df)
Poly_test_pred_df<-getres(svmPoly,test_df)

Radial_train_pred_df<-getres(svmRadial,train_df)
Radial_test_pred_df<-getres(svmRadial,test_df)

Sigmoid_train_pred_df<-getres(svmSigmoid,train_df)
Sigmoid_test_pred_df<-getres(svmSigmoid,test_df)

#首先构建训练集中4个ROC对象
roc_train_linear<-roc(Linear_train_pred_df$label,
                      Linear_train_pred_df$status_pred,
                      auc=T)

roc_train_Poly<-roc(Poly_train_pred_df$label,
                    Poly_train_pred_df$status_pred,
                    auc=T)


roc_train_Radial<-roc(Radial_train_pred_df$label,
                      Radial_train_pred_df$status_pred,
                      auc=T)


roc_train_Sigmoid<-roc(Sigmoid_train_pred_df$label,
                       Sigmoid_train_pred_df$status_pred,
                       auc=T)

RColorBrewer::brewer.pal(4,"Set1")

plot.roc(Linear_train_pred_df$label,
         Linear_train_pred_df$status_pred,
         col="#1c61b6",legacy=T,lwd=2)
lines.roc(Poly_train_pred_df$label,
          Poly_train_pred_df$status_pred,col="#008600")
lines.roc(Radial_train_pred_df$label,
          Radial_train_pred_df$status_pred,col="#E41A1C")
lines.roc(Sigmoid_train_pred_df$label,
          Sigmoid_train_pred_df$status_pred,col="#984EA3")

legend("bottomright",
       legend=c(paste0("train_svmLinear  AUC:",round(roc_train_linear[["auc"]],3)),
                paste0("train_svmPoly     AUC:",round(roc_train_Poly[["auc"]],3)),
                paste0("train_svmRadial   AUC:",round(roc_train_Radial[["auc"]],3)),
                paste0("train_svmSigmoid AUC:",round(roc_train_Sigmoid[["auc"]],3))
       ),
       col=c("#1c61b6","#008600","#E41A1C","#984EA3"),
       lwd=1)
#保存一下模型的原始数据，用于后续的分析
write.table(Linear_train_pred_df, file = "训练集线性核SVM原始数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(Poly_train_pred_df, file = "训练集多项式核SVM原始数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(Radial_train_pred_df, file = "训练集高斯径向基核SVM原始数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(Sigmoid_train_pred_df, file = "训练集sigmoid核SVM原始数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)



#构建测试集中4个ROC对象
roc_test<-lapply(list(Linear_test_pred_df,Poly_test_pred_df,
                      Radial_test_pred_df,Sigmoid_test_pred_df),function(x){
                        roc_res<-roc(x$label,x$status_pred,auc=T)})

roc_test[[1]]
plot.roc(Linear_test_pred_df$label,
         Linear_test_pred_df$status_pred,
         col="#1c61b6",legacy=T)
lines.roc(Poly_test_pred_df$label,
          Poly_test_pred_df$status_pred,col="#008600")
lines.roc(Radial_test_pred_df$label,
          Radial_test_pred_df$status_pred,col="#E41A1C")
lines.roc(Sigmoid_test_pred_df$label,
          Sigmoid_test_pred_df$status_pred,col="#984EA3")
legend("bottomright",
       legend=c(paste0("test_svmLinearAUC:",round(roc_test[[1]][["auc"]],3)),
                paste0("test_svmPolyAUC:",round(roc_test[[2]][["auc"]],3)),
                paste0("test_svmRadialAUC:",round(roc_test[[3]][["auc"]],3)),
                paste0("test_svmSigmoidAUC:",round(roc_test[[4]][["auc"]],3))
       ),
       col=c("#1c61b6","#008600","#E41A1C","#984EA3"),
       lwd=2)

#保存一下模型的原始数据，用于后续的分析
write.table(Linear_train_pred_df, file = "验证集线性核SVM原始数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(Poly_train_pred_df, file = "验证集多项式核SVM原始数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(Radial_train_pred_df, file = "验证集高斯径向基核SVM原始数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(Sigmoid_train_pred_df, file = "验证集sigmoid核SVM原始数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)




####6xgboost####
rm(list = ls())
library(tidyverse)    # for general data wrangling needs
library(gbm)      # for original implementation of regular and stochastic GBMs
library(h2o)      # for a java-based implementation of GBM variants
library(xgboost)  # for fitting extreme gradient boosting
library(rsample)
library(pROC)
library(ggplot2)
library(recipes)
# for data split
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
set.seed(123)  # 设置随机种子，保证结果可复现
credit_df <- data
#经典的37分
index<-sample(1:nrow(credit_df),nrow(credit_df)*0.7,replace=F)#0.7代表37分
train_data<-credit_df[index,]
validation_data<-credit_df[-index,]
dim(train_data)
dim(validation_data)

xgb_prep <- recipe(label ~ ., data =data) %>%
  step_integer(all_nominal()) %>%#数据预处理，将分类变量转变为数值变量
  prep(training = train_data, retain = TRUE) %>%
  juice()
X <- as.matrix(xgb_prep[setdiff(names(xgb_prep), "label")])#特征变量矩阵
Y <- xgb_prep$label
#网格搜索
set.seed(123)
data_xgb <- xgb.cv(
  data = X,
  label = Y,
  nrounds = 6000,
  objective = "binary:logistic",
  early_stopping_rounds = 50, 
  nfold = 10,
  params = list(
    eta = 0.1,
    max_depth = 3,
    min_child_weight = 3,
    subsample = 0.8,
    colsample_bytree = 1.0),
  verbose = 0
)  
# minimum test CV RMSE
min(data_xgb$evaluation_log$test_logloss_mean)
# hyperparameter grid
hyper_grid <- expand.grid(
  eta = 0.01,
  max_depth = 3, 
  min_child_weight = 3,
  subsample = 0.5, 
  colsample_bytree = 0.5,
  gamma = c(0, 1, 10, 100, 1000),
  lambda = c(0, 1e-2, 0.1, 1, 100, 1000, 10000),
  alpha = c(0, 1e-2, 0.1, 1, 100, 1000, 10000),
  logloss = 0,          # a place to dump logloss results
  trees = 0          # a place to dump required number of trees
)

# grid search
for(i in seq_len(nrow(hyper_grid))) {
  set.seed(123)
  m <- xgb.cv(
    data = X,
    label = Y,
    nrounds = 400,
    objective = "binary:logistic",
    early_stopping_rounds = 50, 
    nfold = 10,
    verbose = 0,
    params = list( 
      eta = hyper_grid$eta[i], 
      max_depth = hyper_grid$max_depth[i],
      min_child_weight = hyper_grid$min_child_weight[i],
      subsample = hyper_grid$subsample[i],
      colsample_bytree = hyper_grid$colsample_bytree[i],
      gamma = hyper_grid$gamma[i], 
      lambda = hyper_grid$lambda[i], 
      alpha = hyper_grid$alpha[i]
    ) 
  )
  hyper_grid$logloss[i] <- min(m$evaluation_log$test_logloss_mean)
  hyper_grid$trees[i] <- m$best_iteration
}
hyper_grid %>%
  filter(logloss > 0) %>%
  arrange(logloss) %>%
  glimpse()
# optimal parameter list
params <- list(
  eta = 0.01,
  max_depth = 3,
  min_child_weight = 3,
  subsample = 0.5,
  colsample_bytree = 0.5
)
set.seed(123)
xgb.fit.final <- xgboost(
  params = params,
  data = X,
  label = Y,
  nrounds = 317,
  objective = "binary:logistic",
  verbose = 0
)
train_data$predictions <- predict(xgb.fit.final,as.matrix(xgb_prep[,-13]))
modelrocXGboosttrain<-roc(label~predictions,data=train_data)
modelrocXGboosttrain
round(ci(modelrocXGboosttrain),3)
#验证集
xgb_validation <- recipe(label ~ ., data =data) %>%
  step_integer(all_nominal()) %>%
  prep(training = validation_data, retain = TRUE) %>%
  juice()
X <- as.matrix(xgb_validation[setdiff(names(xgb_validation), "label")])#特征变量矩阵
Y <- xgb_validation$label
set.seed(123)
xgb.fit.final <- xgboost(
  params = params,
  data = X,
  label = Y,
  nrounds = 317,
  objective = "binary:logistic",
  verbose = 0
)
validation_data$predictions <- predict(xgb.fit.final,as.matrix(xgb_validation[,-13]))
modelrocXGboostvalidation<-roc(label~predictions,data=validation_data)
modelrocXGboostvalidation
round(ci(modelrocXGboostvalidation),3)
ggroc(modelrocXGboostvalidation,
      legacy.axes = T,
      size=1.5,
      color="#4b71af")+
  geom_abline(slope =1,intercept = 0,linewidth=1,color="grey")+
  theme_bw()+
  annotate(geom ="text",label=paste0("AUC in the validation-set:",round(modelrocXGboostvalidation$auc,3)),x=0.6,y=0.2)
#测试集
xgb_test <- recipe(label ~ ., data =data) %>%
  step_integer(all_nominal()) %>%
  prep(training = test_data, retain = TRUE) %>%
  juice()
X <- as.matrix(xgb_test[setdiff(names(xgb_test), "label")])#特征变量矩阵
Y <- xgb_test$label
set.seed(123)
xgb.fit.final <- xgboost(
  params = params,
  data = X,
  label = Y,
  nrounds = 317,
  objective = "binary:logistic",
  verbose = 0
)
test_data$predictions <- predict(xgb.fit.final,as.matrix(xgb_test[,-13]))
modelrocXGboosttest<-roc(label~predictions,data=test_data)
modelrocXGboosttest
round(ci(modelrocXGboosttest),3)
ggroc(modelrocXGboosttest,
      legacy.axes = T,
      size=1.5,
      color="#4b71af")+
  geom_abline(slope =1,intercept = 0,linewidth=1,color="grey")+
  theme_bw()+
  annotate(geom ="text",label=paste0("AUC in the test-set:",round(modelrocXGboosttest$auc,3)),x=0.6,y=0.2)
#TCIA数据集
xgb_TCIA <- recipe(label ~ ., data =data) %>%
  step_integer(all_nominal()) %>%
  prep(training = TCIA_data, retain = TRUE) %>%
  juice()
X <- as.matrix(xgb_TCIA[setdiff(names(xgb_TCIA), "label")])#特征变量矩阵
Y <- xgb_TCIA$label
set.seed(123)
xgb.fit.final <- xgboost(
  params = params,
  data = X,
  label = Y,
  nrounds = 317,
  objective = "binary:logistic",
  verbose = 0
)
TCIA_data$predictions <- predict(xgb.fit.final,as.matrix(xgb_TCIA[,-13]))
modelrocXGboostTCIA<-roc(label~predictions,data=TCIA_data)
modelrocXGboostTCIA
round(ci(modelrocXGboostTCIA),3)
ggroc(modelrocXGboostTCIA,
      legacy.axes = T,
      size=1.5,
      color="#4b71af")+
  geom_abline(slope =1,intercept = 0,linewidth=1,color="grey")+
  theme_bw()+
  annotate(geom ="text",label=paste0("AUC in the TCIA-set:",round(modelrocXGboostTCIA$auc,3)),x=0.6,y=0.2)

    


#7K临近####
library(class)
library(caret)
library(sampling)
library(ggplot2)
library(tidyverse)
library(kknn)
library(gmodels)
library(pROC)
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
credit_df <- data
set.seed(123)
index<-sample(1:nrow(credit_df),nrow(credit_df)*0.7,replace=F)#0.7代表37分
train_data<-credit_df[index,]
test_data<-credit_df[-index,]

k = ceiling(sqrt(158))#K等于样本数开根号，568个样本，
knn.pred = knn(train = train_data[, -1], test = test_data[, -1], cl = train_data$label,k = k)
CrossTable(x = test_data$label,y = knn.pred,dnn = c("Actual", "Predicted"),prop.chisq = FALSE)
## knn模型训练
control <- trainControl(method = "cv", number = 10)
grid1 <- expand.grid(.k = seq(2, 24, by = 1))
model <- train(label~ ., train_data, method = "knn", trControl = control, 
         preProcess = c("center","scale"), tuneLength = 5, tuneGrid = grid1)
model
plot(model$results$k, model$results$Accuracy, type = "l", xlab = "K", ylab = "Accuracy",
     lwd = 2)
points(model$results$k, model$results$Accuracy, col = "red", pch = 20, cex = 2)
abline(v = 4, col = "grey", lwd = 1.5)
knn.pred_new = knn(train = train_data[, -1],test = test_data[, -1], cl = train_data$label,k = 4)
CrossTable(x = test_data$label, y = knn.pred_new,dnn = c("Actual", "Predicted"),prop.chisq = FALSE)

#训练集
label <- train_data$label# 测试数据真实值
pred <- predict(model, newdata = train_data)# 测试数据预测值
A=as.data.frame(cbind(label,pred))
roc = roc(A$label,A$pred) #注意一次只能一列
auc <- auc(roc)
auc

write.table(A, file = "训练集K临近ROC数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)

#测试集
label <- test_data$label# 测试数据真实值
pred <- predict(model, newdata = test_data)# 测试数据预测值
A=as.data.frame(cbind(label,pred))
ROC2 = roc(A$label,A$pred) #注意一次只能一列
auc2 <- auc(ROC2)
auc2

write.table(A, file = "测试集K临近ROC数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)

#画图KNN
library(fmsb)
mytable <- table(knn.pred_new, test_data$label)
Kappa.test(mytable, conf.level = 0.95)
kknn.train <- train.kknn(label ~ ., train_data, kmax = 25, distance = 2, kernel = c("rectangular","triangular", "epanechnikov"))
plot(kknn.train)
kknn.pred <- predict(kknn.train, newdata = test_data[, -1])
table(kknn.pred, test_data$label)




#7K均值算法（k-Means）
library(cluster)

#8主成分分析算法（PCA）
library(stats)
princomp(train, cor = TRUE)




####——生存分类—####
#单独验证时test
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
test_df<-data
train_df<-data
  
###01生存Lasso####
rm(list = ls())
library(survival)
library(corrplot)
library(ggplot2)
library(glmnet)
library(timeROC)
set.seed(109)
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
credit_df <- data
index<-sample(1:nrow(credit_df),nrow(credit_df)*0.7,replace=F)#0.7代表37分
train_df<-credit_df[index,]
test_df<-credit_df[-index,]
dim(train_df)#训练集train_df
dim(test_df)#验证集test_df


x<- as.matrix(train_df[,c(5:ncol(train_df))])#是否是3要看前面有几列无用数据
x<- as.matrix(train_df[,c(205:600)])#是否是3要看前面有几列无用数据
y <- cbind(time=train_df$time,status=train_df$status)
lasso <- glmnet(x,y,family = "cox",nlambda=1000,alpha = 1)
#  print(lasso)
#  plot(lasso,xvar="lambda",lable=TRUE,lwd=1.4)
lassoCV <- cv.glmnet(x, y, family = "cox",type.measure = "deviance",nfolds =5)
#  print(lassoCV)
#  plot(lassoCV)
lambda.min<-lassoCV$lambda.min
lambda.1se <- lassoCV$lambda.1se
coef_best<-coef(lassoCV,s ="lambda.min")#lambda.min精准，1se简洁
#coef_best
index<-which(coef_best!=0)#非零系数
coef<-coef_best[index]#对应回归系数
diffvariables=row.names(coef_best)[index]#非零变量


sames=intersect(diffvariables,colnames(train_df))
lassoout0=train_df[ ,sames,drop=F]
lassoout=data.frame(time=train_df$time,status=train_df$status,lassoout0)
myFun=function(x){crossprod(as.numeric(x),coef)}
lassoScore=apply(lassoout0,1,myFun)
outCol=c("time","status",diffvariables)
outTab=cbind(train_df[,outCol],riskScore=as.vector(lassoScore))
write.table(outTab, file = "lasso纳入参数(训练集).xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


ROC=timeROC(T=outTab$time, delta=outTab$status,
            marker=outTab$riskScore, cause=1,
            weighting='aalen',
            times=c(12,36,60), 
            ROC=TRUE)
ROC
plot(ROC, 
     time=12, col="red", lwd=3, title = "")  
plot(ROC,
     time=36, col="blue", add=TRUE, lwd=3)
plot(ROC,
     time=60, col="green", add=TRUE, lwd=3)

legend("bottomright",c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],3)), 
                       paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],3)),
                       paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],3))),
       col=c("red", "blue", "green"),lty=1, lwd=3,bty = "n",cex=1.2) #cex为字体大小


lasso.result.se<-cbind(diffvariables,coef)
lasso.result.se<-as.data.frame(lasso.result.se)
lasso.result.se$coef<-as.numeric(lasso.result.se$coef)
lasso.result.se
ggplot(aes(x=reorder(diffvariables,coef),y=coef,fill=diffvariables),data=lasso.result.se)+
  geom_col()+
  coord_flip()+
  theme_bw()+
  labs(x="")+
  ggtitle("LASSO identified variables")+
  scale_fill_brewer(palette = "Set3")+
  theme(legend.position = "")
write.table(outTab, file = "lasso纳入参数(训练集).xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


#验证集
lassoout1=test_df[ ,sames,drop=F]
lassoScore=apply(lassoout1,1,myFun)
outCol2=c("time","status",diffvariables)
outTab2=cbind(test_df[,outCol2],riskScore=as.vector(lassoScore))
ROC=timeROC(T=outTab2$time, delta=outTab2$status,
            marker=outTab2$riskScore, cause=1,
            weighting='aalen',
            times=c(12,36,60), 
            ROC=TRUE)
ROC
plot(ROC,time=12, col="red", lwd=3, title = "")  
plot(ROC,time=36, col="blue", add=TRUE, lwd=3)
plot(ROC,time=60, col="green", add=TRUE, lwd=3)

legend("bottomright",c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],3)), 
                       paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],3)),
                       paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],3))),
       col=c("red", "blue", "green"),lty=1, lwd=3,bty = "n",cex=1.2) #cex为字体大小


#特征排序
ggplot(df1, aes(x = group,y=value,fill = variable))+
  geom_col(position = 'stack', width = 0.6)+
  geom_errorbar(aes(ymin=xx-sd,ymax=xx+sd),
                width=0.1,linewidth=0.8)+
  scale_y_continuous(expand = c(0,0),limits = c(0,70))+
  labs(x=NULL,y=NULL)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.92,0.88),
        legend.background = element_blank())+
  scale_fill_manual(values = c("#0099cc","#ff9933"))


#4.3特征排序
A = list()#新增一个空表
A$Genes=lassoGene
A$geneCoef=actCoef
A$group="A"
A=as.data.frame(A)
#绘图
ggplot(A, aes(Genes,geneCoef))+ theme_classic()+#主题
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", width = 0.15,size=1)+#误差线
  geom_bar(aes(fill=group),color="black",stat="summary",fun=mean,position="dodge",size=0.5)+#柱状图
  geom_jitter(color="black",size = 2.5,width = 0.2,alpha=0.9)+#散点图
  geom_jitter(color="black",size = 2.5,width = 0.2,alpha=0.9)



###02随机生存森林RSF####
rm(list = ls())
library(randomForest)
library(randomForestSRC)
library(timeROC)
library(ggplot2)#二分类化

set.seed(170)
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
train_df <- data
index<-sample(1:nrow(data),nrow(data)*0.55,replace=F)#0.7代表37分
train_df<-data[index,]
test_df<-data[-index,]


rfs <- rfsrc(Surv(time,status) ~.,train_df,
             ntree = 1000,
             importance =TRUE,
             seed=123)
#print(rfs)
#plot(rfs)

head(rfs$predicted) # 这个就是预测值，当概率用
head(rfs$predicted.oob)
tmp.train <- train_df
tmp.train$riskScore <- rfs$predicted
dim(tmp.train)
#ggplot(tmp.train, aes(x = riskScore))+geom_histogram(fill = "steelblue",color="white")+theme_classic()
median(tmp.train$riskScore)#方法有很多，比如中位数、平均数、百分位数等，建议都试一下，选择效果最好的那个，然后你再想办法解释为什么这样选！
risk_group <- ifelse(tmp.train$riskScore > median(tmp.train$riskScore), "high_risk","low_risk")
table(risk_group)
tmp.train$riskGroup <- risk_group
dim(tmp.train)
#write.table(tmp.train, file = "复发GTV+CTV随机森林(训练组).xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


ROC <- timeROC(T = tmp.train$time,   
               delta = tmp.train$status,   
               marker = tmp.train$riskScore, # 这里就是risk score
               cause = 1,           
               weighting = "marginal",   
               times = c(18,35,60),       
               iid = TRUE)
ROC
plot(ROC,time=18, col="red", lwd=3, title = "")  
plot(ROC,time=35, col="blue", add=TRUE, lwd=3)
plot(ROC,time=60, col="green", add=TRUE, lwd=3)
legend("bottomright",c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],3)), 
                       paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],3)),
                       paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],3))),
       col=c("red", "blue", "green"),lty=1, lwd=3,bty = "n",cex=1.2) #cex为字体大小


#测试集
rfs.test <- predict(rfs, newdata = test_df)
rfs.test
tmp.test <- test_df
tmp.test$riskScore <- rfs.test$predicted
dim(tmp.test)
ROC <- timeROC(T= tmp.test$time,   
                    delta= tmp.test$status,   
                    marker=tmp.test$riskScore,   
                    cause=1,                
                    weighting="marginal",   
                    times=c(12,36,60),       
                    iid=TRUE)
ROC
plot(ROC,time=12, col="red", lwd=3, title = "")  
plot(ROC,time=36, col="blue", add=TRUE, lwd=3)
plot(ROC,time=60, col="green", add=TRUE, lwd=3)

legend("bottomright",c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],3)), 
                       paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],3)),
                       paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],3))),
       col=c("red", "blue", "green"),lty=1, lwd=3,bty = "n",cex=1.2) #cex为字体大小

write.table(tmp.test, file = "复发GTV+CTV随机森林(验证组).xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


#####特征排序####
print(sort(rfs$importance,decreasing = TRUE)[1:10])
#变量筛选
res.val <- rfs$importance[rfs$importance > 0.005]#筛选重要性大于0.005的变量
#计数，保留10个
plot.df <- data.frame(val = names(res.val),importance = res.val)
library(ggplot2)
#方法1
ggplot(aes(x=reorder(val,importance),y=importance,fill=val),data=plot.df)+
  geom_col()+
  coord_flip()+
  theme_bw()+
  labs(x="")+
  ggtitle("Feature sorting in RSF")+
  scale_fill_brewer(palette = "Set3")+
  theme(legend.position ="")
#方法2
ggplot(plot.df, aes(importance, reorder(val, importance)))+
  geom_col(fill="skyblue")+
  labs(y=NULL)+
  theme_bw()
write.csv(plot.df,file = "重要特征排序复发随机森林.csv")

##这几个变量就是我们最终的变量，接下来你可以使用这几个变量重新建立模型##。
#计算变量重要性的可信区间
rfs.smp <- subsample(rfs, B=100, #重抽样次数
                     subratio = 0.2) #重抽样比例
plot.subsample(rfs.smp)





###03GBM梯度提升机####
rm(list = ls())
library(survival)
library(gbm)
library(limma)
set.seed(121)#
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
credit_df <- data
index<-sample(1:nrow(credit_df),nrow(credit_df)*0.8,replace=F)#0.7代表37分
train_df<-credit_df[index,]
test_df<-credit_df[-index,]


#构建GBM预后模型
model <- gbm(Surv(time, status) ~ .,
             distribution = "coxph",
             data = train_df,
             n.trees = 3000,
             shrinkage = 0.005,
             interaction.depth = 2,
             n.minobsinnode = 5,
             cv.folds = 5
)
#计算GBM风险评分
best.iter <- gbm.perf(model, plot.it = TRUE, method = "cv")
pred.train <- predict(model, train_df, n.trees = best.iter)
pred.test <- predict(model, test_df, n.trees = best.iter)
GBMtrain_data=data.frame(train_df,riskscore=pred.train )
GBMtest_data=data.frame( test_df,riskscore=pred.test)

ROC <- timeROC(T = GBMtrain_data$time,   
               delta = GBMtrain_data$status,   
               marker = GBMtrain_data$riskscore, #这里就是risk score
               cause = 1,           
               weighting = "marginal",   
               times = c(12,36,60),       
               iid = TRUE)
ROC
plot(ROC,time=12, col="red", lwd=3, title = "")  
plot(ROC,time=36, col="blue", add=TRUE, lwd=3)
plot(ROC,time=60, col="green", add=TRUE, lwd=3)
legend("bottomright",c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],3)), 
                       paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],3)),
                       paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],3))),
       col=c("red", "blue", "green"),lty=1, lwd=3,bty = "n",cex=1.2) #cex为字体大小

     
#验证集
ROC <- timeROC(T = GBMtest_data$time,   
               delta = GBMtest_data$status,   
               marker = GBMtest_data$riskscore, # 这里就是risk score
               cause = 1,           
               weighting = "marginal",   
               times = c(12,36,60),       
               iid = TRUE)
ROC
plot(ROC,time=12, col="red", lwd=3, title = "")  
plot(ROC,time=36, col="blue", add=TRUE, lwd=3)
plot(ROC,time=60, col="green", add=TRUE, lwd=3)
legend("bottomright",c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],3)), 
                       paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],3)),
                       paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],3))),
       col=c("red", "blue", "green"),lty=1, lwd=3,bty = "n",cex=1.2) #cex为字体大小

#保存
write.csv(GBMtrain_data,"GBMtrain_data.csv")
write.csv(GBMtest_data,"GBMtest_data.csv")


#_特征排序####
rfs1 = summary.gbm(model, plotit=TRUE)
print(sort(rfs1$rel.inf,decreasing = TRUE)[1:10])
res.val <- rfs1$rel.inf[rfs1$rel.inf > 0.631]#通过修改这个值，筛选变量
#只要前十个

#将cl样本行名的顺序按exp2列名进行匹配：
features <- rfs1[match(rfs1$rel.inf,res.val),]
#删除包含空值的行
features<-na.omit(features)
plot.df <- data.frame(val = features$var, importance = res.val)
head(plot.df)
library(ggplot2)

#方法1
ggplot(aes(x=reorder(val,importance),y=importance,fill=val),data=plot.df)+
  geom_col()+
  coord_flip()+
  theme_bw()+
  labs(x="")+
  ggtitle("Top 10 feature in GBM")+
  scale_fill_brewer(palette = "Set3")+
  theme(legend.position ="")
#方法2
ggplot(plot.df, aes(importance, reorder(val, importance)))+
  geom_col(fill="skyblue")+
  labs(y=NULL)+
  theme_bw()
write.csv(plot.df,file = "重要特征排序生存GBM.csv")





###04XGboost####
rm(list = ls())
set.seed(125)
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
credit_df <- data
index<-sample(1:nrow(credit_df),nrow(credit_df)*0.67,replace=F)#0.7代表37分
data<-credit_df[index,]
data2<-credit_df[-index,]
train_df <- data
test_df <- data2

library(xgboost)
label <- ifelse(data$status== 1, data$time, -data$time)
x_train <- as.matrix(data[,c(3:ncol(data))])
x_label <- label
x_val <- xgb.DMatrix(x_train,label =x_label)
# train surv xgboost
xgb_param<-list("objective"="survival:cox",
                "eval_metric" = "cox-nloglik",
                alpha = 0.1,     # 设置L1正则化参数alpha
                lambda = 0.1     # 设置L2正则化参数lambda
)
xgb.model<-xgb.train(params = xgb_param,
                     data=x_val,
                     nrounds=300,#迭代次数，数值越大，运行时间越长
                     eta=0.7,
                     watchlist=list(val2=x_val))
#结果中每行的train-cox-nloglik代表模型效果，一般越小越好。
#但是要注意，在这里我们观察到在第n次之后该值逐渐增大，这表示模型出现了过度拟合
#因此，我们选择迭代次数n来构建模型
pred <- predict(xgb.model,x_val)
XGBOOSTtrain_data=data.frame(train_df,riskscore=pred)
ROC <- timeROC(T = XGBOOSTtrain_data$time,   
               delta = XGBOOSTtrain_data$status,   
               marker = XGBOOSTtrain_data$riskscore, # 这里就是risk score
               cause = 1,           
               weighting = "marginal",   
               times = c(12,36,60),       
               iid = TRUE)
ROC
plot(ROC,time=12, col="red", lwd=3, title = "")  
plot(ROC,time=36, col="blue", add=TRUE, lwd=3)
plot(ROC,time=60, col="green", add=TRUE, lwd=3)
legend("bottomright",c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],3)), 
                       paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],3)),
                       paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],3))),
       col=c("red", "blue", "green"),lty=1, lwd=3,bty = "n",cex=1.2) #cex为字体大小


#验证组
label2 <- ifelse(data2$status== 1, data2$time, -data2$time)
x_test <- as.matrix(data2[,c(3:ncol(data2))])
x_label2 <- label2
x_val2 <- xgb.DMatrix(x_test,
                      label =x_label2)
pred <- predict(xgb.model,x_val2)
XGBOOSTtest_data=data.frame(test_df,riskscore=pred)
ROC <- timeROC(T = XGBOOSTtest_data$time,   
               delta = XGBOOSTtest_data$status,   
               marker = XGBOOSTtest_data$riskscore, # 这里就是risk score
               cause = 1,           
               weighting = "marginal",   
               times = c(12,36,60),       
               iid = TRUE)
ROC
plot(ROC,time=12, col="red", lwd=3, title = "")  
plot(ROC,time=36, col="blue", add=TRUE, lwd=3)
plot(ROC,time=60, col="green", add=TRUE, lwd=3)

legend("bottomright",c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],3)), 
                       paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],3)),
                       paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],3))),
       col=c("red", "blue", "green"),lty=1, lwd=3,bty = "n",cex=1.2) #cex为字体大小


#_特征排序shap####
library(shapviz)
library(ggplot2)
library(xgboost)
X_small <- x_train[sample(nrow(x_train), 50L), ]#100L是人数
shp <- shapviz(xgb.model, X_pred = data.matrix(X_small), X = X_small)
sv_importance(shp, kind = "bar")+  theme_bw()
sv_importance(shp,kind = "beeswarm")+  theme_bw()





###05survivalSVM####
rm(list = ls())
library(survival)
library(survivalsvm)
library(timeROC)
set.seed(123)
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
index<-sample(1:nrow(data),nrow(data)*0.8,replace=F)#0.9代表37分
train_df<-data[index,]
test_df<-data[-index,]

n <- nrow(train_df)
train.index <- sample(1:n,n,replace = FALSE)
survsvm<- survivalsvm(Surv(time, status) ~ ., gamma.mu = 1,#参数gamma.mu用于调节正则化参数。
                       subset = train.index, data = train_df,opt.meth = "quadprog",
                       type = "regression", #模式1
                       kernel = "lin_kernel") #方法1
survsvm<- survivalsvm(Surv(time, status) ~ .,gamma.mu = 1,#参数gamma.mu用于调节正则化参数。
                       subset = train.index, data = train_df,opt.meth = "quadprog",
                       type = "vanbelle1", #模式2
                       kernel = "lin_kernel",#方法1
                       diff.meth = 'makediff1')
survsvm<- survivalsvm(Surv(time, status) ~ .,gamma.mu = c(0.01, 0.1),
                       subset = train.index, data = train_df,opt.meth = "quadprog",
                       type = "hybrid", #模式3
                       kernel = "add_kernel",#方法1
                       diff.meth = 'makediff1')
#模式type包括4个：'regression'，'vanbelle1'和'vanbelle2'，'hybrid'，
#方法kernel包括四个：linear kern ('lin_kernel')；additive kernel ('add_kernel')
#radial basis kernels ('rbf_kernel')和the polynomial kernel ('poly_kernel')

riskscore <- predict(survsvm,train_df,type = "risk")
riskscore1 <- riskscore[["predicted"]]
riskscore2 <- riskscore1[1, ]
SVMtrain_data=data.frame(train_df,riskscore=-riskscore2)
ROC <- timeROC(T = SVMtrain_data$time,   
               delta = SVMtrain_data$status,   
               marker = SVMtrain_data$riskscore,
               cause = 1,           
               weighting = "marginal",   
               times = c(12,36,60),       
               iid = TRUE)
ROC
plot(ROC,time=12, col="red", lwd=3, title = "")  
plot(ROC,time=36, col="blue", add=TRUE, lwd=3)
plot(ROC,time=60, col="green", add=TRUE, lwd=3)
legend("bottomright",c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],3)), 
                       paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],3)),
                       paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],3))),
       col=c("red", "blue", "green"),lty=1, lwd=3,bty = "n",cex=1.2) #cex为字体大小

#验证组
n <- nrow(test_df)
riskscore <- predict(survsvm,test_df,type = "risk")
riskscore1 <- riskscore[["predicted"]]
riskscore2 <- riskscore1[1, ]
SVMtrain_data=data.frame(test_df,riskscore=-riskscore2)
ROC <- timeROC(T = SVMtrain_data$time,   
               delta = SVMtrain_data$status,   
               marker = SVMtrain_data$riskscore,
               cause = 1,           
               weighting = "marginal",   
               times = c(22,36,60),       
               iid = TRUE)
ROC
plot(ROC, time=22, col="red", lwd=3, title = "")  
plot(ROC,time=36, col="blue", add=TRUE, lwd=3)
plot(ROC,time=60, col="green", add=TRUE, lwd=3)
legend("bottomright",c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],3)), 
                       paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],3)),
                       paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],3))),
       col=c("red", "blue", "green"),lty=1, lwd=3,bty = "n",cex=1.2) #cex为字体大小

write.csv(SVMtrain_data,"SVMtrain_data.csv")
write.csv(SVMtrain_data,"SVMtrain_data.csv")



###—自动机器学习—####
library(h2o)
library(pROC)
h2o.init()
credit_df<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
set.seed(1)
index<-sample(1:nrow(credit_df),nrow(credit_df)*1,replace=F)#0.7代表37分
train_df<-credit_df[index,]
test_df<-credit_df[-index,]
df=as.h2o(train_df)#训练集
df2=as.h2o(test_df)#验证集


colnames(df)
df$status <- h2o.asfactor(df$status)
y="status"#这个注意点，我改了
x <- setdiff(names(df),y)
dl <- h2o.deeplearning(y=y,   # 因变量
                       x=x,   # 预测变量
                       training_frame = df,  # 数据集
                       epochs=200,   # 迭代次数
                       variable_importances=T, # 变量重要性
                       model_id="dl"  # 模型ID
)
dl_perf <- h2o.performance(dl,train = T)
dl_perf
#AUC
h2o.auc(dl_perf)
plot(dl_perf)
plot(dl)
h2o.varimp_plot(dl)#变量重要性
#部分依赖图（PDP）
h2o.partialPlot(dl,df,"INS")#这里可以查看每个参数的PDP
#ICE图
h2o.ice_plot(dl,df,"TNM")#这里可以查看每个参数的ICE

#原始数据,训练集
predictions <- h2o.predict(dl, newdata = df)
A= as.data.frame(predictions)#变表格
A$status=train_df$status
roc = roc(A$status,A$p0) #注意一次只能一列
roc_auc <- auc(roc)
roc_auc
plot(roc, col="red", lwd=3, title = "")  

write.table(A, file = "训练集原始数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


#测试集原始数据
predictions2 <- h2o.predict(dl, newdata = df2)
A= as.data.frame(predictions2)#变表格
A$status=test_df$status
roc = roc(A$status,A$p0) #注意一次只能一列
roc_auc <- auc(roc)
roc_auc

write.table(A, file = "训练集原始数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


#真实版的验证集AUC会比较低，这里是重新建模
df2$status <- h2o.asfactor(df2$status)
dl2 <- h2o.deeplearning(y=y,   # 因变量
                       x=x,   # 预测变量
                       training_frame = df2,  # 数据集
                       epochs=200,   # 迭代次数
                       variable_importances=T, # 变量重要性
                       model_id="dl" )
dl_perf2 <- h2o.performance(dl2,train = T)
#AUC
h2o.auc(dl_perf2)
plot(dl_perf2)
plot(dl2)
#原始数据
predictions2 <- h2o.predict(dl2, newdata = df2)
A= as.data.frame(predictions2)#变表格
A$status=test_df$status
roc = roc(A$status,A$p0) #注意一次只能一列
roc_auc <- auc(roc)
roc_auc

write.table(A, file = "训练集原始数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)






#不指定模型####
library(h2o)
library(pROC)
h2o.init()
credit_df<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
index<-sample(1:nrow(credit_df),nrow(credit_df)*0.7,replace=F)#0.7代表37分
train_df<-credit_df[index,]
test_df<-credit_df[-index,]
df=as.h2o(train_df)#训练集train_df
df2=as.h2o(test_df)#验证集test_df

am <- h2o.automl(y="status",#二分类
                 training_frame =df,#max_runtime_secs：限制时间
                 max_models = 1)#max_models限制模型数量
b <- h2o.get_leaderboard(am)
b
h2o.get_best_model(am,criterion = "auc")
pred <- h2o.predict(am@leader,df)
pred
pred<- as.data.frame(pred)#变表格
roc = roc(pred$predict,pred$p1) #注意一次只能一列
roc_auc <- auc(roc)
roc_auc
plot(roc, col="red", lwd=3, title = "") 

pred2 <- h2o.predict(am@leader,df2)
pred2
pred2<- as.data.frame(pred2)#变表格
roc = roc(pred2$predict,pred2$p1) #注意一次只能一列
roc_auc <- auc(roc)
roc_auc
plot(roc, col="red", lwd=3, title = "")  


best <- h2o.get_best_model(am)
best
perf <- h2o.performance(best)
perf
h2o.auc(perf)
plot(perf,type="roc")







####—mlr机器学习—####
#BiocManager::install("survcomp")
library(survival)
library(survivalsvm)
library(dplyr)
library(Hmisc)
library(timeROC)
library(survcomp)
library(mlr3proba)#install.packages("mlr3proba", repos = "https://mlr-org.r-universe.dev")
library(mlr3verse)#remotes::install_github("mlr-org/mlr3verse")
library(mlr3extralearners)#remotes::install_github("mlr-org/mlr3extralearners")
#原始数据集
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header=TRUE,row.names=1,check.names=F)
#创建任务
task = as_task_surv(data, time = "time",event = "status")#生存分类
task
task$head()
autoplot(task,rhs="Stage")#绘制KM曲线
#数据划分
split = partition(task)
#训练模型的同时，还可以预测模型
p = lrn("surv.coxph")$#在此处修改模型
        train(task, split$train)$predict(task, split$test)
#训练集
t.train<-data[split$train,2]#time
s.train<-data[split$train,1]#status
#查看C-index指数
Cindex<-1-rcorr.cens(p$crank,Surv(t.train,s.train))[[1]]
Cindex
#AUC
ROC<-timeROC(T = t.test, delta = s.test, marker = p$crank,times= c(12,36,60), cause=1)
ROC
plot(ROC,time=12, col="red", lwd=2, title = "")  
plot(ROC,time=36, col="blue", add=TRUE, lwd=2)
plot(ROC,time=60, col="green", add=TRUE, lwd=2)
legend("bottomright",c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],3)), 
                       paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],3)),
                       paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],3))),
       col=c("red", "blue", "green"),lty=1, lwd=2,bty = "n") 
#验证集
t.test<-data[split$test,2]#time
s.test<- data[split$test,1]#status
#查看C-index指数
Cindex<-1-rcorr.cens(p$crank,Surv(t.test,s.test))[[1]]
Cindex
#AUC
ROC<-timeROC(T = t.test, delta = s.test, marker = p$crank,times= c(12,36,60), cause=1)
ROC
plot(ROC,time=12, col="red", lwd=2, title = "")  
plot(ROC,time=36, col="blue", add=TRUE, lwd=2)
plot(ROC,time=60, col="green", add=TRUE, lwd=2)
legend("bottomright",c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],3)), 
                       paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],3)),
                       paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],3))),
       col=c("red", "blue", "green"),lty=1, lwd=2,bty = "n") 
#查看次模型的其他评价指标有哪些
as.data.table(mlr_measures)[task_type == "surv", c("key", "predict_type")]
#查看某一项
p$score(msrs("surv.calib_beta"))

#我们还可以同时设定多个生存模型，然后进行基准测试
pipe = as_learner(ppl("distrcompositor",
  learner = lrn("surv.coxph"),
  estimator = "kaplan",form = "ph"))
pipe$id = "Coxnet"
#查看所有模型
as.data.table(mlr_learners)[task_type == "surv",key]
#设定学习器集
learners = c(lrns(c("surv.gbm","surv.kaplan"))
             ,pipe)#这里可以比较多个
# 基准测试
bmr = benchmark(benchmark_grid(task,learners,rsmp("cv", folds = 3)))
#设置多个评估指标
msr_txt = c("surv.rcll", "surv.cindex", "surv.dcalib")#评估指标
measures = msrs(msr_txt)
#基准测试
bmr$aggregate(measures)[, c("learner_id", ..msr_txt)]
#箱线图可视化，在mlr3，你可以实现万物皆一行
autoplot(bmr, measure = msr("surv.cindex"))



