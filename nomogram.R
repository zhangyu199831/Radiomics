#linshi
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE,row.names = 1)
complex1<- decision_curve(status~Pathscore,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex2<- decision_curve(status~CCINC,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)


List1<- list(complex1,complex2)
plot_decision_curve(List1,curve.names= c("Pathscore",'CCINC'),
                    cost.benefit.axis =T,confidence.intervals =FALSE,standardize = FALSE) 





####Spearman检验####
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE,row.names = 1)
rt <- data
norm_result <- apply(rt, 2, function(x) shapiro.test(x)$p.value)
norm_feature <- rt[which(norm_result >= 0.05)]
cor_nor <- cor(norm_feature, method = "pearson")
cor_all <- cor(rt, method = "spearman")
num_nor <- dim(cor_nor)[1]
cor_all[1:num_nor, 1:num_nor] <- cor_nor
cor_all[upper.tri(cor_all)] <- 0
diag(cor_all) <- 0
data_reduce = rt[, !apply(cor_all, 1, function(x) any(abs(x) > 0.9))]
dim(data_reduce)
write.table(data1, file = "Spearman后数据.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)



####批量相关性分析####
library(tidyverse)
library(stringr)
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
m1=data[,c(1:2)]#影像组学参数
m2=data[,c(2:48)]#基因
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
colnames(Toal) <- c('Radiomics','RNA','pvalue','corr')
colnames(Toal2) <- c('Radiomics','RNA','pvalue','corr')

write.table(Toal, file = "放疗抵抗炎症因子的相关性分析.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
#write.table(Toal2, file = "放疗抵抗炎症因子的相关性分析过滤.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)



####相关性热图####
library(corrplot)
library(ggplot2)
library(ggcorrplot)
library(vcd)
library(psych)
library(ggrepel)
#数据导入#
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
#data=data[,3:5]
data<-as.matrix(data) #利用as.matrix()将所需数据集转换为matrix格式，才可在corrplot中跑
data=data.frame(scale(data))#数据标准化
#相关性计算
data<-cor(data,method="spearman") #pearson，spearman和kendall
round(data, 2)#保留两位小数

#相关性热图绘制#
ggcorrplot(data, method="circle") #圆圈大小变化
#带数值
ggcorrplot(data, method = "square", #"square", "circle"
           type ="upper" , #full完全(默认)，lower下三角，upper上三角
           ggtheme = ggplot2::theme_minimal,
           title = "",
           show.legend = TRUE,  #是否显示图例。
           legend.title = "Corr", #指定图例标题。
           show.diag =T ,    #FALSE显示中间
           colors = c("blue", "white", "red"), #需要长度为3的颜色向量，同时指定low,mid和high处的颜色。
           outline.color = "gray", #指定方形或圆形的边线颜色。
           hc.order = F,  #是否按hclust(层次聚类顺序)排列。
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
corrplot(data, method="circle", #square方形，ellipse, 椭圆形，number数值，shade阴影，color颜色，pie饼图
         type="full",  #full完全(默认)，lower下三角，upper上三角
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
         tl.pos = NULL,  #指定文本标签(变量名称)相对绘图区域的位置，为"lt"(左侧和顶部),"ld"(左侧和对角线),"td"(顶部和对角线),"d"(对角线),"n"(无);当type="full"时默认"lt"。当type="lower"时默认"ld"。当type="upper"时默认"td"。
         tl.cex = 1.5,  #设置文本标签的大小
         tl.col = "black", #设置文本标签的颜色。
         cl.pos = NULL )



####线性相关可视化####
library("ggpubr")
x = "METTL3"
y = "EMP1"

data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
ggscatter(data, x = x, y = y, 
          add = "reg.line", conf.int = TRUE, color = "#00CED1", size = 3,
          cor.coef = TRUE, cor.method = "spearman",#pearson
          add.params = list(color = "#00CED1", fill = "#D3D3D3", size = 2), 
          xlab = x, ylab = y)
#换颜色
ggscatter(data, x = x, y = y, 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",#pearson
          add.params = list(color = "#0000FF", fill = "#007FFF", size = 1), 
          xlab = x, ylab = y)
#换颜色
ggscatter(data, x = x, y = y, 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",#pearson
          add.params = list(color = "red", fill = "#808080", size = 1), 
          xlab = x, ylab = y)+
          theme_bw()+#带边框
          theme(panel.grid=element_blank(),#去除网格
          legend.text=element_text(colour= 'black',size=10),
          axis.text= element_text(colour= 'black',size=10),
          axis.line= element_line(colour= 'black'),
          panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))

#第二种方法
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
A<-data$IGF2BP3
B<-data$METTL3
data <- data.frame(A,B)
head(data)
cor.test(A,B,data=data)
df_cor<-lm(B~A,data=data)
summary(df_cor)
#最简单展示——plot函数
plot(A, B,xlab = "株高", ylab = "根毛数量",pch = 16, frame = T)
# 添加回归线
abline(lm(B ~ A), col = "red")
#加载包
library(ggplot2)
library(ggprism)
p1<-ggplot(data,aes(x=A,y=B,color="orange"))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(size=3,shape=16)+#绘制点图并设定大小
  theme(panel.grid = element_blank())+
  labs(x="株高",y="根毛数量")+#x、y轴标题
  geom_smooth(method='lm', se=FALSE, color='turquoise4')+#添加回归线
  geom_text(aes(x=55,y=124,label="R^2=0.88\ny=0.75x+58.96"),
            color="red",family = "serif",fontface = "plain",size = 5)+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.8, # 坐标轴的粗细
              axis_text_angle = 45)+ # 可选值有 0，45，90，270
  scale_fill_prism(palette = "candy_bright")+
  theme(legend.position = 'none')#去除图例
p1
#回归诊断
par(mfrow=c(2,2))
plot(df_cor)  #绘制回归诊断图

##根据图可以看出第10、13，14个点偏离较远，需要剔除重新进行分析
data2<-data[c(-10,-13,-14),] 
cor.test(A,B,data=data2)
df_cor2<-lm(B~A,data=data2)
summary(df_cor2)

p2<-ggplot(data2,aes(x=A,y=B,color="orange"))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(size=3,shape=16)+#绘制点图并设定大小
  theme(panel.grid = element_blank())+
  labs(x="株高",y="根毛数量")+#x、y轴标题
  geom_smooth(method='lm', se=FALSE, color='turquoise4')+#添加回归线
  geom_text(aes(x=55,y=124,label="R^2=0.98\ny=0.73x+60.47"),
            color="red",family = "serif",fontface = "plain",size = 5)+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.8, # 坐标轴的粗细
              axis_text_angle = 45)+ # 可选值有 0，45，90，270
  scale_fill_prism(palette = "candy_bright")+
  theme(legend.position = 'none')#去除图例

p2


p2+geom_point(aes(x=69,y=104),shape=8,color="red",size=3)+
  geom_point(aes(x=75,y=120),shape=8,color="red",size=3)+
  geom_point(aes(x=68,y=116),shape=8,color="red",size=3)






####ICC####
rm(list = ls())
library(readxl)
library(psych)
#先改路径，放到文件所在路径
df_icc1 <- read_excel("ICC瘤周.xlsx")#第一位医生首次勾画
df_icc2 <- read_excel("ICC1瘤周.xlsx")#第一位医生二次勾画
df_icc3 <- read_excel("ICC2瘤周.xlsx")#第二位医生首次勾画

n <- dim(df_icc1)[1]
p <- dim(df_icc1)[2]
df_intra <- rbind(df_icc1,df_icc2)
df_inter <- rbind(df_icc1,df_icc3)
#组内一致性ICC
icc_all_intra <- apply(df_intra, 2, function(x) ICC(x = data.frame(x[1:n], x[(n+1):(2*n)]), lmer = F)$results[1,2])
icc_all_intra <- as.data.frame(icc_all_intra)#变表格
sum(icc_all_intra$icc > 0.75, na.rm=TRUE)
#组间一致性分析!
icc_all_inter <- apply(df_inter, 2, function(x) ICC(x = data.frame(x[1:n], x[(n+1):(2*n)]), lmer = F)$results[2,2])
icc_all_inter  <- data.frame(icc_all_inter)
sum(icc_all_inter$icc > 0.75, na.rm=TRUE)

df_icc_all_intra <- data.frame(feature = names(df_intra),icc = icc_all_intra,row.names = 1:ncol(df_intra))
df_icc_all_inter <- data.frame(feature = names(df_inter),icc = icc_all_inter,row.names = 1:ncol(df_inter))
#记得改保存位置
write.csv(df_icc_all_intra,"瘤周组内ICC.csv")
write.csv(df_icc_all_inter,"瘤周组间ICC.csv")



#多个ROC####
library(pROC)
library(survminer)
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
colnames(data)
#1、同个数据集
ROC<-roc(label~DIE+Dose+points,data=data,aur=TRUE,ci=TRUE)
ROC
#结合ggplot2绘图函数
ggroc(ROC, legacy.axes = TRUE)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
  theme_bw()+ggtitle('ROC')+ggsci::scale_color_lancet()+
  annotate("text",x=0.85,y=0.25,label=paste("DIE", round(ROC$DIE$auc,3)))+
  annotate("text",x=0.85,y=0.2,label=paste("Dose", round(ROC$Dose$auc,3)))+
  annotate("text",x=0.85,y=0.15,label=paste("points", round(ROC$points$auc,3)))



#2、不同数据集
ROC1 <- roc(a$label,a$"X0", ci = T)
ROC2 <- roc(b$label,b$"X0", ci = T)

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
#四条线
plot(smooth(ROC3), col = "darkgreen", lwd = 2, add = TRUE)
legend("bottomright", cex = 1.1,
       legend = c("s100b (AUC: 0.731)", "ndka (AUC: 0.612)", "wfns (AUC: 0.824)"),
       col = c("red", "blue", "darkgreen"), lty = 1, lwd = 2)  #定义标签

#3、不同数据集,同个excel
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
colnames(data)
set.seed(108)
index<-sample(1:nrow(data),nrow(data)*0.62,replace=F)#0.7代表37分
train_data<-data[index,]
test_data<-data[-index,]

ROC1 <- roc(train_data$label,train_data$"A", ci = T)#记得改标题
ROC2 <- roc(test_data$label,test_data$"A", ci = T)
auc1 <- auc(ROC1)
auc1
auc2 <- auc(ROC2)
auc2
plot(smooth(ROC1), col = "red", lwd = 2)
plot(smooth(ROC2), col = "blue", lwd = 2, add = TRUE)
#两条线
legend("bottomright", cex = 1.1, title="LASSO",
       legend = c("Training Set (AUC: 0.850)", "Validation Set (AUC: 0.807)"),
       col = c("red", "blue"), lty = 1, lwd = 2)  #定义标签
#六条线
ROC1 <- roc(train_data$label,train_data$"CD56", ci = T)#记得改标题
ROC2 <- roc(test_data$label,test_data$"A", ci = T)

auc1 <- auc(ROC1)
auc1
auc2 <- auc(ROC2)
auc2
plot(smooth(ROC1), col = "red", lwd = 2)
plot(smooth(ROC2), col = "blue", lwd = 2, add = TRUE)
#两条线
legend("bottomright", cex = 1.1, title="LASSO",
       legend = c("Training Set (AUC: 0.850)", "Validation Set (AUC: 0.807)"),
       col = c("red", "blue"), lty = 1, lwd = 2)  #定义标签



####ROC的参数####
library(reportROC)
credit_df<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
set.seed(121)  
index<-sample(1:nrow(credit_df),nrow(credit_df)*0.6,replace=F)#0.7代表37分
data<-data1[index,]

colnames(data)
ROC=reportROC(gold=data$label,
          predictor=data$Dose,
          important="se",
          plot=TRUE)
ROC$AUC#曲线下面积；
ROC$SEN#灵敏度；
ROC$SPE#特异度；
(as.numeric(ROC$NPV)+as.numeric(ROC$PPV))/2#准确率Accuracy
ROC$PLR#阳性似然比；
ROC$NLR#阴性似然比；
ROC$PPV#阳性预测值；
ROC$NPV#阴性预测值；



write.csv(ROC, "ROC.info.csv")



####——Logistic—####
#devtools::install_github('yikeshu0611/ggDCA')
#BiocManager::install("devtools")
rm(list=ls())
library(lattice)
library(Formula)
library(ggplot2)
library(Hmisc)
library(SparseM)
library(foreign)
library(rms)
library(caret)
library(car)
library(pROC)
library(timeROC)
library(ggDCA)
library(rmda)
library(dcurves)
library(dplyr)
library(tidyverse)
library(broom)
library(scales)
####批量单因素logistc####
rm(list=ls())
Pima.tr <- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
options(datadist='Pima.tr')

#多个单因素Logistic回归
names(Pima.tr)
vars = names(Pima.tr)[2:6]#前7个参数
vars
model = tibble(vars) %>% 
  mutate(model = map(Pima.tr[vars], 
                     ~ glm(label~ .x, data = Pima.tr, family = binomial()))) %>% 
  mutate(result = map(model, tidy),
         OR = map(model, ~ exp(coef(.x))),
         OR_ci = map(model, ~ exp(confint(.x)))) %>% 
  select(-model) %>% 
  unnest(c(result, OR, OR_ci))

model = model %>% 
  mutate(OR_ci %>% as_tibble()) %>% 
  select(-OR_ci) %>% 
  rename(LL = V1, UL = V2) %>% 
  mutate(across(term, ~ str_remove(.x, '.x'))) %>% 
  filter(if_all(term, ~ !.x=='(Intercept)')) %>% 
  mutate(`OR(95%CI)` = str_c(round(OR,2), ' (', round(LL,2), '-', round(UL,2), ')')) %>% 
  select(vars, term, `OR(95%CI)`, p.value, OR, LL, UL, ) %>% 
  mutate(p.value = pvalue(p.value))
model
#model文件就是结果
write.table(model, file = "单因素分析结果.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


#单个单因素Logistic回归
model <- glm(label ~ N, data = Pima.tr, family = binomial(link = 'logit'))
tidy(model) # 也可以summary(model)
exp(coef(model)) # 计算OR
exp(confint(model)) # OR的95%置信区间


####列线图####
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
#需要分类的话
data$Hormones<- factor(data$Hormones,labels=c("低","高")) 
data$Dose<-factor(data$Dose,labels=c("低","高")) 

dd=datadist(data)
options(datadist="dd") 
colnames(data)
ckf<-lrm(label ~ Hormones+Dose,x=T,y=T,data=data)
#查看单因素Cox分析结果，最下方可见其对应的Coef、p值
print(ckf)

nom<-nomogram(ckf,fun=plogis,lp=T, funlabel = "Risk",maxscale =100,
              fun.at=c('0.8','0.7','0.6','0.5','0.4','0.3','0.2',"0.1"))
#方法1
plot((nom),xfrac=.3)
#方法2,可以修改参数
nomogram <- nomogram(ckf,fun=function(x)1/(1+exp(-x)), ##逻辑回归计算公式
                     fun.at = c(0.01,0.1,0.5,0.9,0.99),#风险轴刻度
                     funlabel = "Prob of recurrence", #风险轴便签
                     lp=F,  ##是否显示系数轴
                     conf.int = F, ##每个得分的置信度区间，用横线表示,横线越长置信度越
                     abbrev = F)#是否用简称代表因子变量

#方法3
plot(nom,lplabel="Linear Predictor",xfrac = 0.3, # 左侧标签距离坐标轴的距离
     varname.label = TRUE,tcl = -0.2, # 刻度长短和方向 
     lmgp = 0.1, # 坐标轴标签距离坐标轴远近
     points.label ='Points',total.points.label = 'Total Points',
     cap.labels = FALSE,cex.var = 1, # 左侧标签字体大小
     cex.axis = 1, col.grid = c('pink','grey')) # 竖线颜色


#查看参数
ckf
#计算总分
library(nomogramFormula)
results<-formula_rd(nomogram=nom)
data$points<-points_cal(formula = results$formula,rd=data)
head(data)
write.table(data,"联合模型ROC数据1.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)



#森林图####
library(forestmodel) # 森林图
library(ggplot2) # 使用到里面的数据集“msleep”
library(dplyr) # 数据处理
# 建立回归模型
logit_model <- glm(label ~ Hormones+Dose, data = data, family = "binomial")
# 制作森林图
forest_model(logit_model)
#加彩
forest_model(logit_model,format_options=forest_model_format_options(colour = "steelblue", # 修改颜色
             shape = 15, # 改变形状
             text_size = 5, # 字体文本的大小
             point_size = 5, # 森林图中“方框”的大小
             banded = TRUE), # 如“FALSE”，可去掉森林图背景条带
             factor_separate_line = TRUE)




##校准曲线####
cal<-calibrate(ckf,method='boot',B=500)

plot(cal)
plot(cal,xlab="Nomogram-Predicted Probability")

#方法2
pred <- predict(ckf,data,type = "fitted")
val.prob(pred, data$label,cex = 0.8)#可以调整字体大小


#决策曲线####
complex1<- decision_curve(label~Hormones,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex2<- decision_curve(label~Dose,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex3<- decision_curve(label~Hormones+Dose,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex4<- decision_curve(label~N+AFP,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex5<- decision_curve(label~N+AFP+Pathscore,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex6<- decision_curve(label~DMIr +SCC+f + Radscore,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)

List1<- list(complex1,complex2,complex3)
par(mai=c(1,1,0.5,0.5))
plot_decision_curve(List1,curve.names= c('Hormones','Dose','Model'),
                    cost.benefit.axis =T,
                    confidence.intervals =FALSE,standardize = FALSE) 
#如果不需要图例，可以加上legend.position = "none"。如下：
plot_decision_curve(List1,curve.names= c('DMIr','SCC','f','Radscore','Model'),
                    legend.position = "none",
                    cost.benefit.axis =T,
                    confidence.intervals =FALSE,standardize = FALSE)


#临床影响曲线####
plot_clinical_impact(complex3,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits= 8,col = c('red','blue'),
                     confidence.intervals= T)



####IDI及NRI####
library(PredictABEL)
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
#方法1
logistic.model<- 
  list(Basic_model=glm(label~D,data=data,family=binomial), 
         New_model=glm(label~N+AFP+Pathscore,data=data, family = binomial))#加一个参数
reclassification(data=data,cOutcome=1,#Y在第几列这里就填几
   predrisk1 = fitted(logistic.model[["Basic_model"]]), 
   predrisk2 = fitted(logistic.model[["New_model"]]), 
   cutoff = c(0, 0.40, 1))
#1表示连续NRI分类变量；2表示NRI连续变量；3表示IDI.

#方法2
#拟合logistic回归模型
oldmodel=glm(label~N,family=binomial(),data=data,x=TRUE)
newmodel=glm(label~N+AFP+Pathscore,family=binomial(),data=data,x=TRUE)
#分别提取两个模型的预测值。
pold<-oldmodel$fitted.values
pnew<-newmodel$fitted.values
#设置数据集为矩阵格式
mydata2<-as.matrix(data)
reclassification(data=mydata2,cOutcome=2,predrisk1=pold,predrisk2=pnew,cutoff=c(0,0.6,1))


####ROC#### 
library(pROC)
library(ggplot2)

rocobj1<-roc(data$label,-data$Hormones)
rocobj2<-roc(data$label,data$Dose)
rocobj3<-roc(data$label,data$points)
rocobj4<-roc(data$label,data$points2)
rocobj5<-roc(data$label,data$points)
auc(rocobj1)
auc(rocobj2)
auc(rocobj3)
auc(rocobj4)
auc(rocobj5)
g7<-ggroc(list(Hormones=rocobj1,Dose=rocobj2,points=rocobj3))
plot(g7)
g7 + theme_classic() + ggtitle("ROC curve") +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), 
               color="grey", linetype="dashed")




####——COX——####
rm(list=ls())
library(foreign)
library(survival)
library(rms)
library(regplot)
library(survminer)
library(Formula)
library(nomogramFormula)
library(ggplot2)
library(Hmisc)
library(SparseM)
library(caret)
library(car)
library(pROC)
library(timeROC)
library(ggDCA)
library(dcurves)
library(dplyr)
library(tidyverse)
library(rmda)
#求取截断值####
#为连续变量求截断值,把连续变量放到低3-11列
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
colnames(data)
var =colnames(data[,3:7])#提取前11个列名
res.cut <- surv_cutpoint(data,time="time",event ="status",variables =var)
summary(res.cut)
cut <- surv_categorize(res.cut)

plot(res.cut, "FLT1")
write.table(cut, file = "截断值.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
#手动修订文件为0,1

fit <- survfit(Surv(time, status)~SFRP2, data = cut)#KM曲线
ggsurvplot(fit, data = cut) #得到基础的生存曲线啦
ggsurvplot(fit, pval = TRUE, palette = "jco", 
           data = cut, legend = c(0.8, 0.8),
           ylab="OS (%)",xlab = "Month", #更改横纵坐标
           conf.int = TRUE, #给生存曲线添加上置信区间
           legend.labs = c("High","Low"), #在图中添加图例
           title ="SFRP2",risk.table = F)


#KM拼在一起####
data=cut
splots = list()
colnames(data)
genes<-var 

cl4=data[,(1:2)]#提取前3列作为临床信息
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


arrange_ggsurvplots(splots,nrow = 1, ncol = 5)#一行几个，一列几个




####批量单因素COX####
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
colnames(data)
can <- colnames(data)[4:31]#批量计算那些列
can
uni_sur <- sapply(can, function(x) as.formula(paste('Surv(time, status)~', x)))
uni_cox <- lapply(uni_sur, function(x){coxph(x, data = data)})
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
write.table(result, file = "临床单因素cox结果.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)



nom <- nomogram(coxm, 
                fun = list(
                  function(x) surv(12, x), 
                  function(x) surv(36, x), 
                  function(x) surv(60, x)), 
                lp = F, 
                funlabel = c("1-year survival", "3-year survival", "5-year survival"),
                maxscale = 10, 
                fun.at = c('0.99', '0.9', '0.7', '0.5', '0.3', '0.1', '0.01'),
                varlabel = c("FIGO", "Differentiation", "Radscore")  # 使用varlabel参数 instead of termlabel
)

####列线图####
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
dd=datadist(data)
options(datadist="dd") 
colnames(data)
coxm <- cph(Surv(time,status) ~ ADC+FIGO+Radscore+EQD2,
            x=T, y=T, surv=T, 
            data=data)
print(coxm)
surv<-Survival(coxm)
nom<- nomogram(coxm,fun=list(function(x) surv(12, x), function(x) surv(36, x),
               function(x)surv(60, x)), lp=F, funlabel=c("1-year survival", "3-year survival", "5-year survival"),
               maxscale=10, fun.at=c('0.99','0.9','0.7','0.5','0.3','0.1','0.01'))

#方法1
plot((nom),xfrac=0.2) 


#方法2，分类展示
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
data$TNM<- factor(data$TNM,labels=c("I","II","III","IV")) 
data$BCLC<-factor(data$BCLC,labels=c("0","A","B","C")) 
data$Pathscore<-factor(data$Pathscore,labels=c("Low","High")) 
data$Differentiation<-factor(data$Differentiation,labels=c("G1","G2","G3")) 


data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
data$FIGO<- factor(data$FIGO,labels=c("I","II","III","IV")) 
data$SCC<-factor(data$SCC,labels=c("Low","High")) 
data$EQD2<-factor(data$EQD2,labels=c("Low","High")) 
data$Radscore<-factor(data$Radscore,labels=c("Low","High")) 


str(data)
dd=datadist(data)
options(datadist="dd") 
head(dd)
coxm <- cph(Surv(time, status) ~ SCC+FIGO+EQD2+Radscore,
         x=T, y=T, surv=T, data=data)
print(coxm)
surv<- Survival(coxm)
nom<- nomogram(coxm, fun=list(function(x) surv(12, x), function(x) surv(36, x),
               function(x)surv(60, x)), lp=F, funlabel=c("1-year OS", "3-year OS", "5-year OS"),
               maxscale=10, fun.at=c('0.99','0.9','0.7','0.5','0.3','0.1','0.01'))
plot(nom)


#方法3
plot(nom, 
     lplabel="Linear Predictor",xfrac = 0.3, # 左侧标签距离坐标轴的距离
     varname.label = TRUE,tcl = -0.2, # 刻度长短和方向 
     lmgp = 0.1, # 坐标轴标签距离坐标轴远近
     points.label ='Points',total.points.label = 'Total Points',
     cap.labels = FALSE,cex.var = 1, # 左侧标签字体大小
     cex.axis = 1, col.grid = c('pink','grey')) # 竖线颜色

#方法4
regplot(coxm,plots = c("violin", "boxes"), ##连续性变量形状，可选"no plot" "density" "boxes" "ecdf" "bars" "boxplot" "violin" "bean" "spikes"；分类变量的形状，可选"no plot" "boxes" "bars" "spikes"
        observation = F,center = T, # 对齐变量
        subticks = T,droplines = T,#是否画竖线
        title = "nomogram",points = T, # 截距项显示为0-100
        odds = T, # 是否显示OR值
        showP = T, # 是否显示变量的显著性标记
        rank = "sd", # 根据sd给变量排序
        interval="confidence", # 展示可信区间
        clickable = F)


#计算总分，并作为新变量添加进原始数据
results<-formula_rd(nomogram=nom)
data$points<-points_cal(formula = results$formula,rd=data)

#write.table(data, file = "列线图score.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


#计算c_index
f<-coxph(Surv(time,status==1)~A+B+C+D, data=data)
sum.surv<-summary(f)
c_index<-sum.surv$concordance
c_index
#95%CI
confitOR<-exp(confint(f))
confitOR
#OR值
valueOR<-exp(coef(f))
valueOR
#95%CI
CI95=paste(round(valueOR,3),"(",round(confitOR[,1],3),"-",round(confitOR[,2], 3),")",sep="")
CI95


####森林图####
library(survminer)
library(forplo)
library(survival)
#黑白森林图
colnames(data)
multiCox <- coxph(Surv(time,status) ~FIGO+Differentiation+Radscore,data)
multiCox
multicoxF<-ggforest(multiCox, data = data,  #数据集
           main = 'Hazard ratio of multi cox',  #标题
           cpositions = c(0.05, 0.15, 0.35),  #前三列距离
           fontsize = 0.9, #字体大小
           refLabel = 1, noDigits = 3 )#保留HR值以及95%CI的小数位数
multicoxF



##校准曲线####
coxm<-cph(Surv(time,status==1)~SCC+FIGO+SCC+Radscore+EQD2, x=T,y=T,data=data,surv=T)

cal1<-calibrate(coxm,cmethod='KM',method='boot',u=22,m=13,B=200)
plot(cal1,xlim = c(0,1),ylim =c(0,1))
cal2<-calibrate(coxm,cmethod='KM',method='boot',u=36,m=14,B=200)
plot(cal2,xlim = c(0,1),ylim =c(0,1))
cal3<-calibrate(coxm,cmethod='KM',method='boot',u=48,m=15,B=200)

plot(cal3,lwd=0,lty=0,xlim = c(0,1),ylim =c(0,1),##x轴和y轴范围
     xlab="Nomogram-Predicted Probability of DFS",
     ylab="Actual OS (%)",
     mgp = c(2, 1, 0))#控制坐标轴的位置

lines(cal1[,c("mean.predicted","KM")],type="b",lwd=2,col=c("red"),pch=16)
lines(cal2[,c("mean.predicted","KM")],type="b",lwd=2,col=c("blue"),pch=16)
lines(cal3[,c("mean.predicted","KM")],type="b",lwd=2,col=c("green"),pch=16)

abline(0,1,lty=3,lwd=2,col=c("black"))
legend("bottomright",c(paste0("1 year  "),paste0("3 year  "),paste0("5 year  ")),
       col=c("red", "blue", "green"),lty=1, lwd=3,bty = "n",cex=1.5) #cex为字体大小



data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
## 决策曲线####
complex1<- decision_curve(status~SCC,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex2<- decision_curve(status~FIGO,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex3<- decision_curve(status~EQD2,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex4<- decision_curve(status~Radscore,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex5<- decision_curve(status~SCC+FIGO+EQD2,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex6<- decision_curve(status~SCC+FIGO+Radscore+EQD2,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
#先看一下前四个
List2<- list(complex1,complex2,complex3,complex4,complex5,complex6)
par(mai=c(1,1,0.5,1))
plot_decision_curve(List2,curve.names= c('SCC',"FIGO","EQD2",'Radscore','SCC+FIGO+EQD2',"Model"),
                    cost.benefit.axis =FALSE,confidence.intervals =FALSE,standardize = FALSE) 


complex5<- decision_curve(status~T+TNM,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex6<- decision_curve(status~T+TNM+Radscore+A,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex7<- decision_curve(status~T+TNM+Pathscore+C,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex8<- decision_curve(status~TNM+Pathscore+Radscore+B+C,data = data,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)

#后四个
List2<- list(complex5,complex6,complex7,complex8)
par(mai=c(1,1,0.5,1))
plot_decision_curve(List2,curve.names= c("T+TNM","T+TNM+Radscore","T+TNM+Pathscore",'Combined Model'),
                    cost.benefit.axis =FALSE,
                    confidence.intervals =FALSE,standardize = FALSE) 

#加上所有
List2<- list(complex1,complex2,complex3,complex4,complex5,complex6,complex7,complex8)
par(mai=c(1,1,0.5,1))
plot_decision_curve(List2,curve.names= c('T','TNM','Radscore',"Pathscore","T+TNM",
                     "T+TNM+Radscore","T+TNM+Pathscore",'Combined Model'),
                    cost.benefit.axis =FALSE,
                    confidence.intervals =FALSE,standardize = FALSE) 


#如果不需要图例，可以加上legend.position = "none"。如下：
plot_decision_curve(List2,curve.names= c('1','2','3','4','5','6'),
                    cost.benefit.axis =FALSE,
                    confidence.intervals =FALSE,standardize = FALSE,legend.position = "none")

#临床影像曲线####
plot_clinical_impact(complex5,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits= 8,col = c('red','blue'),#不要图例
                     confidence.intervals= T,legend.position = "none")

plot_clinical_impact(complex6,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits= 8,col = c('red','blue'),
                     confidence.intervals= T)


#timeROC####
#data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
colnames(data)

ROC=timeROC(T=data$time, delta=data$status,
            marker=data$Transcore, cause=1,
            weighting='marginal',
            times=c(12,36,55), 
            ROC=TRUE)
ROC
plot(ROC,time=12, col="red", lwd=3, title = "")  
plot(ROC,time=36, col="blue", add=TRUE, lwd=3)
plot(ROC,time=55, col="green", add=TRUE, lwd=3)

legend("bottomright",c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],3)), 
                       paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],3)),
                       paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],3))),
       col=c("red", "blue", "green"),lty=1, lwd=3,bty = "n",cex=1.2) #cex为字体大小


####IDI及NRI####
library(PredictABEL)
library(survIDINRI)
b<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
colnames(b)
sur.IDI.INF <- IDI.INF(indata = b[,c("time","status")],
    covs0 = b[,c("A","C")],            
    covs1 = b[,c("A","D","C")], 
    t0=60,npert=300,alpha=0.05)#一年12，3年36,5年60
IDI.INF.OUT(sur.IDI.INF)
#M1的结果表示IDI。M2的结果表示连续NRI。M3的结果表示中值改善风险评分.
#IDI图
IDI.INF.GRAPH(sur.IDI.INF,
              main = paste0("5-year IDI",
              round(sur.IDI.INF$point$IDI,digits = 4)))

#NRI图
library(nricens)
newfit <- coxph(Surv(time,status)~A+B+C,data=b,x=T)
oldfit <- coxph(Surv(time,status)~A+C,data=b,x=T)

NRI<- nricens(mdl.std =oldfit,mdl.new = newfit,#新模型
                       t0 =  60, #预测年限,记得改
                       cut = c(0.05),#截断值
                       updown = 'diff',
                       niter = 10)#迭代次数(bootstrap)
legend("bottomrigh", cex = 1.3, title="3-year NRI",legend=NA, bty="n",)  #定义标签

caseDat <- NRI$rtab.case
controlDat <- NRI$rtab.ctrl
(5-18)/(5+5+18+87)+(13-1)/(1+1+13+28)  
NRI$nri 



####动态列线图####
#shinyapps的链接
rsconnect::setAccountInfo(name='zy1998',
                          token='B36D0AC69E3761A19990D4E5460CA4ED',
                          secret='YQAda9YoTqby7O6f5yc+a8I8ljvn+DekL1qBuBKx')
rsconnect::setAccountInfo(name='zy0301',
                          token='714634041095F02EEF975F3E88E41B61',
                          secret='<SECRET>')
#进行lrm的分析时，原装DynNom包中对lrm对象进行处理时存在一处bug
#下方为修改后的代替函数DynNom_czx_lrm及DNbuilder_czx_lrm，请直接运行加载使用。参数与原装DynNom及DNbuilder函数一致。
{
  #进行lrm的分析时，原装DynNom包中对lrm对象进行处理时存在一处bug
  #代码中term原意应为因变量及自变量的名称的character形式，但是代码内仅生成自变量的部分，缺少因变量。
  #报错显示为
  DynNom_czx_lrm<-function (model, data = NULL, clevel = 0.95, m.summary = c("raw","formatted"), covariate = c("slider", "numeric"), 
                            ptype = c("st","1-st"), DNtitle = NULL, DNxlab = NULL, DNylab = NULL, DNlimits = NULL,
                            KMtitle = NULL, KMxlab = NULL, KMylab = NULL) 
  {
    mclass <- getclass.DN(model)$model.class
    mfamily <- getclass.DN(model)$model.family
    if (mclass %in% c("coxph", "cph")) {
      Surv.in <- length(model$terms[[2]]) != 1
    }
    if (mclass %in% c("ols", "Glm", "lrm", "cph")) {
      model <- update(model, x = T, y = T)
    }
    if (!is.data.frame(data)) {
      if (any(class(try(getdata.DN(model), silent = TRUE)) == 
              "try-error")) {
        stop("Dataset needs to be provided in a data.frame format")
      }
      else {
        data <- getdata.DN(model)
      }
    }
    covariate <- match.arg(covariate)
    m.summary <- match.arg(m.summary)
    ptype <- match.arg(ptype)
    if (mclass %in% c("lm", "glm", "ols", "Glm", "lrm", "gam", 
                      "Gam")) {
      Terms.T <- all(all.vars(model$terms) %in% names(data))
    }
    if (mclass %in% c("coxph")) {
      if (Surv.in) {
        Terms.T <- all(all.vars(model$terms)[-c(1:2)] %in% 
                         names(data))
      }
      else {
        Terms.T <- all(all.vars(model$terms)[-1] %in% names(data))
      }
    }
    if (mclass %in% c("cph")) {
      Terms.T <- all(names(model$Design$units) %in% names(data))
    }
    if (!Terms.T) 
      stop("Error in model syntax: some of model's terms do not match to variables' name in dataset")
    if (!is.null(DNlimits) & !length(DNlimits) == 2) 
      stop("A vector of 2 is required as 'DNlimits'")
    if (is.null(DNtitle)) 
      DNtitle <- "Dynamic Nomogram"
    if (is.null(DNxlab)) {
      DNxlab <- ifelse((mclass %in% c("glm") & mfamily %in% 
                          c("binomial", "quasibinomial")) | mclass == "lrm" | 
                         mclass %in% c("coxph", "cph"), "Probability", "Response variable")
    }
    if (mclass %in% c("coxph", "cph")) {
      if (is.null(KMtitle)) {
        if (ptype == "st") {
          KMtitle <- "Estimated Survival Probability"
        }
        else {
          KMtitle <- "Estimated Probability"
        }
      }
      if (is.null(KMxlab)) {
        KMxlab <- "Follow Up Time"
      }
      if (is.null(KMylab)) {
        if (ptype == "st") {
          KMylab <- "S(t)"
        }
        else {
          KMylab <- "F(t)"
        }
      }
    }
    if (mclass %in% c("lm", "glm", "ols", "Glm", "lrm", "gam", 
                      "Gam")) {
      DynNom.core_czx_lrm(model, data, clevel, m.summary, covariate, 
                          DNtitle, DNxlab, DNylab, DNlimits)
    }
    if (mclass %in% c("coxph", "cph")) {
      DynNom.surv(model, data, clevel, m.summary, covariate, 
                  ptype, DNtitle, DNxlab, DNylab, KMtitle, KMxlab, 
                  KMylab)
    }
  }
  
  DynNom.core_czx_lrm<-function (model, data, clevel, m.summary, covariate, DNtitle,DNxlab, DNylab, DNlimits) 
  {
    mclass <- getclass.DN(model)$model.class
    mfamily <- getclass.DN(model)$model.family
    if (mclass %in% c("lm", "ols")) {
      mlinkF <- function(eta) eta
    }
    else {
      mlinkF <- ifelse(mclass == "lrm", function(mu) plogis(mu), 
                       model$family$linkinv)
    }
    input.data <- NULL
    old.d <- NULL
    if (mclass %in% c("ols", "Glm", "lrm")) {
      model <- update(model, x = T, y = T)
    }
    allvars <- all.vars(model$terms)
    if (mclass %in% c("ols", "lrm", "Glm")) {
      #进行lrm的分析时，此处term原意应为因变量及自变量的名称的character形式，但是此处仅生成自变量的部分，缺少因变量。
      #因此debug此处，增加一项因变量名称在前
      terms <- model$Design$assume[model$Design$assume != 
                                     "interaction"]
      terms<-c("asis",terms)
      names(terms) = c(model$terms[[2]],model$Design$name[model$Design$assume != "interaction"])
      #结束debug
      
      if (mclass %in% c("Glm")) {
        terms <- c(attr(attr(model$model, "terms"), "dataClasses")[1], 
                   terms)
      }
      else {
        terms <- c(attr(model$terms, "dataClasses")[1], 
                   terms)
      }
    }
    if (mclass %in% c("lm", "glm", "gam", "Gam")) {
      if (length(attr(model$terms, "dataClasses")) == length(allvars)) {
        terms <- attr(model$terms, "dataClasses")
        names(terms) = allvars
      }
      else {
        terms <- attr(model$terms, "dataClasses")[which(names(attr(model$terms, 
                                                                   "dataClasses")) %in% allvars)]
      }
    }
    if (terms[[1]] == "logical") 
      stop("Error in model syntax: logical form for response not supported")
    terms[terms %in% c("numeric", "asis", "polynomial", "integer", 
                       "double", "matrx") | grepl("nmatrix", terms, fixed = T) | 
            grepl("spline", terms, fixed = T)] = "numeric"
    terms[terms %in% c("factor", "ordered", "logical", "category", 
                       "scored")] = "factor"
    resp <- terms[1]
    names(resp) <- allvars[1]
    if ("(weights)" %in% names(terms)) {
      preds <- as.list(terms[-c(1, length(terms))])
    }
    else {
      preds <- as.list(terms[-1])
    }
    names(preds) <- allvars[-1]
    for (i in 1:length(preds)) {
      if (preds[[i]] == "numeric") {
        i.dat <- which(names(preds[i]) == names(data))
        preds[[i]] <- list(v.min = floor(min(na.omit(data[, 
                                                          as.numeric(i.dat)]))), v.max = ceiling(max(na.omit(data[, 
                                                                                                                  as.numeric(i.dat)]))), v.mean = zapsmall(mean(data[, 
                                                                                                                                                                     as.numeric(i.dat)], na.rm = T), digits = 4))
        next
      }
      if (preds[[i]] == "factor") {
        i.dat <- which(names(preds[i]) == names(data))
        if (mclass %in% c("ols", "Glm", "lrm", "cph")) {
          preds[[i]] <- list(v.levels = model$Design$parms[[which(names(preds[i]) == 
                                                                    names(model$Design$parms))]])
        }
        else {
          preds[[i]] <- list(v.levels = model$xlevels[[which(names(preds[i]) == 
                                                               names(model$xlevels))]])
        }
      }
    }
    if (!is.null(DNlimits) & !length(DNlimits) == 2) 
      stop("A vector of 2 is required as 'DNlimits'")
    if (is.null(DNlimits)) {
      if ((mclass %in% c("glm") & mfamily %in% c("binomial", 
                                                 "quasibinomial")) | mclass == "lrm") {
        limits0 <- c(0, 1)
      }
      else {
        if (mclass %in% c("lm", "glm", "gam", "Gam")) {
          limits0 <- c(mean(model$model[, names(resp)]) - 
                         3 * sd(model$model[, names(resp)]), mean(model$model[, 
                                                                              names(resp)]) + 3 * sd(model$model[, names(resp)]))
        }
        if (mclass %in% c("ols", "lrm", "Glm")) {
          limits0 <- c(mean(model$y) - 3 * sd(model$y), 
                       mean(model$y) + 3 * sd(model$y))
        }
      }
      if (mclass %in% c("glm", "Glm") & mfamily %in% c("poisson", 
                                                       "quasipoisson", "Gamma")) {
        limits0[1] <- 0
      }
    }
    neededVar <- c(names(resp), names(preds))
    data <- data[, neededVar]
    input.data <- data[0, ]
    runApp(list(ui = bootstrapPage(fluidPage(titlePanel(DNtitle), 
                                             sidebarLayout(sidebarPanel(uiOutput("manySliders"), 
                                                                        uiOutput("setlimits"), actionButton("add", "Predict"), 
                                                                        br(), br(), helpText("Press Quit to exit the application"), 
                                                                        actionButton("quit", "Quit")), mainPanel(tabsetPanel(id = "tabs", 
                                                                                                                             tabPanel("Graphical Summary", plotlyOutput("plot")), 
                                                                                                                             tabPanel("Numerical Summary", verbatimTextOutput("data.pred")), 
                                                                                                                             tabPanel("Model Summary", verbatimTextOutput("summary"))))))), 
                server = function(input, output) {
                  observe({
                    if (input$quit == 1) stopApp()
                  })
                  limits <- reactive({
                    if (!is.null(DNlimits)) {
                      limits <- DNlimits
                    } else {
                      if (input$limits) {
                        limits <- c(input$lxlim, input$uxlim)
                      } else {
                        limits <- limits0
                      }
                    }
                  })
                  output$manySliders <- renderUI({
                    slide.bars <- list()
                    for (j in 1:length(preds)) {
                      if (terms[j + 1] == "factor") {
                        slide.bars[[j]] <- list(selectInput(paste("pred", 
                                                                  j, sep = ""), names(preds)[j], preds[[j]]$v.levels, 
                                                            multiple = FALSE))
                      }
                      if (terms[j + 1] == "numeric") {
                        if (covariate == "slider") {
                          slide.bars[[j]] <- list(sliderInput(paste("pred", 
                                                                    j, sep = ""), names(preds)[j], min = preds[[j]]$v.min, 
                                                              max = preds[[j]]$v.max, value = preds[[j]]$v.mean))
                        }
                        if (covariate == "numeric") {
                          slide.bars[[j]] <- list(numericInput(paste("pred", 
                                                                     j, sep = ""), names(preds)[j], value = zapsmall(preds[[j]]$v.mean, 
                                                                                                                     digits = 4)))
                        }
                      }
                    }
                    do.call(tagList, slide.bars)
                  })
                  output$setlimits <- renderUI({
                    if (is.null(DNlimits)) {
                      setlim <- list(checkboxInput("limits", "Set x-axis ranges"), 
                                     conditionalPanel(condition = "input.limits == true", 
                                                      numericInput("uxlim", "x-axis upper", 
                                                                   zapsmall(limits0[2], digits = 2)), numericInput("lxlim", 
                                                                                                                   "x-axis lower", zapsmall(limits0[1], 
                                                                                                                                            digits = 2))))
                    } else {
                      setlim <- NULL
                    }
                    setlim
                  })
                  a <- 0
                  new.d <- reactive({
                    input$add
                    input.v <- vector("list", length(preds))
                    for (i in 1:length(preds)) {
                      input.v[[i]] <- isolate({
                        input[[paste("pred", i, sep = "")]]
                      })
                      names(input.v)[i] <- names(preds)[i]
                    }
                    out <- data.frame(lapply(input.v, cbind))
                    if (a == 0) {
                      input.data <<- rbind(input.data, out)
                    }
                    if (a > 0) {
                      if (!isTRUE(compare(old.d, out))) {
                        input.data <<- rbind(input.data, out)
                      }
                    }
                    a <<- a + 1
                    out
                  })
                  p1 <- NULL
                  old.d <- NULL
                  data2 <- reactive({
                    if (input$add == 0) return(NULL)
                    if (input$add > 0) {
                      if (!isTRUE(compare(old.d, new.d()))) {
                        isolate({
                          mpred <- getpred.DN(model, new.d())$pred
                          se.pred <- getpred.DN(model, new.d())$SEpred
                          if (is.na(se.pred)) {
                            lwb <- "No standard errors"
                            upb <- paste("by '", mclass, "'", sep = "")
                            pred <- mlinkF(mpred)
                            d.p <- data.frame(Prediction = zapsmall(pred, 
                                                                    digits = 3), Lower.bound = lwb, Upper.bound = upb)
                          } else {
                            lwb <- sort(mlinkF(mpred + cbind(1, 
                                                             -1) * (qnorm(1 - (1 - clevel)/2) * 
                                                                      se.pred)))[1]
                            upb <- sort(mlinkF(mpred + cbind(1, 
                                                             -1) * (qnorm(1 - (1 - clevel)/2) * 
                                                                      se.pred)))[2]
                            pred <- mlinkF(mpred)
                            d.p <- data.frame(Prediction = zapsmall(pred, 
                                                                    digits = 3), Lower.bound = zapsmall(lwb, 
                                                                                                        digits = 3), Upper.bound = zapsmall(upb, 
                                                                                                                                            digits = 3))
                          }
                          old.d <<- new.d()
                          data.p <- cbind(d.p, counter = 1, count = 0)
                          p1 <<- rbind(p1, data.p)
                          p1$counter <- seq(1, dim(p1)[1])
                          p1$count <- 0:(dim(p1)[1] - 1)%%11 + 1
                          p1
                        })
                      } else {
                        p1$count <- seq(1, dim(p1)[1])
                      }
                    }
                    rownames(p1) <- c()
                    p1
                  })
                  output$plot <- renderPlotly({
                    if (input$add == 0) return(NULL)
                    if (is.null(new.d())) return(NULL)
                    coll = c("#0E0000", "#0066CC", "#E41A1C", "#54A552", 
                                      "#FF8000", "#BA55D3", "#006400", "#994C00", 
                                      "#F781BF", "#00BFFF", "#A9A9A9")
                                      lim <- limits()
                                      yli <- c(0 - 0.5, 10 + 0.5)
                                      dat2 <- data2()
                                      if (dim(data2())[1] > 11) {
                                        input.data = input.data[-c(1:(dim(input.data)[1] - 
                                                                        11)), ]
                                        dat2 <- data2()[-c(1:(dim(data2())[1] - 11)), 
                                        ]
                                        yli <- c(dim(data2())[1] - 11.5, dim(data2())[1] - 
                                                   0.5)
                                      }
                                      in.d <- input.data
                                      xx <- matrix(paste(names(in.d), ": ", t(in.d), 
                                                         sep = ""), ncol = dim(in.d)[1])
                                      Covariates <- apply(xx, 2, paste, collapse = "<br />")
                                      p <- ggplot(data = dat2, aes(x = Prediction, 
                                                                   y = counter - 1, text = Covariates, label = Prediction, 
                                                                   label2 = Lower.bound, label3 = Upper.bound)) + 
                                        geom_point(size = 2, colour = coll[dat2$count], 
                                                   shape = 15) + ylim(yli[1], yli[2]) + coord_cartesian(xlim = lim) + 
                                        labs(title = paste(clevel * 100, "% ", "Confidence Interval for Response", 
                                                           sep = ""), x = DNxlab, y = DNylab) + theme_bw() + 
                                        theme(axis.text.y = element_blank(), text = element_text(face = "bold", 
                                                                                                 size = 10))
                                      if (is.numeric(dat2$Upper.bound)) {
                                        p <- p + geom_errorbarh(xmax = dat2$Upper.bound, 
                                                                xmin = dat2$Lower.bound, size = 1.45, height = 0.4, 
                                                                colour = coll[dat2$count])
                                      } else {
                                        message(paste("Confidence interval is not available as there is no standard errors available by '", 
                                                      mclass, "' ", sep = ""))
                                      }
                                      gp <- ggplotly(p, tooltip = c("text", "label", 
                                                                    "label2", "label3"))
                                      gp$elementId <- NULL
                                      gp
                  })
                  output$data.pred <- renderPrint({
                    if (input$add > 0) {
                      if (nrow(data2()) > 0) {
                        if (dim(input.data)[2] == 1) {
                          in.d <- data.frame(input.data)
                          names(in.d) <- names(terms)[2]
                          data.p <- cbind(in.d, data2()[1:3])
                        }
                        if (dim(input.data)[2] > 1) {
                          data.p <- cbind(input.data, data2()[1:3])
                        }
                      }
                      stargazer(data.p, summary = FALSE, type = "text")
                    }
                  })
                  output$summary <- renderPrint({
                    if (m.summary == "formatted") {
                      if (mclass == "lm") {
                        stargazer(model, type = "text", omit.stat = c("LL", 
                                                                      "ser", "f"), ci = TRUE, ci.level = clevel, 
                                  single.row = TRUE, title = paste("Linear Regression:", 
                                                                   model$call[2], sep = " "))
                      }
                      if (mclass %in% c("glm")) {
                        stargazer(model, type = "text", omit.stat = c("LL", 
                                                                      "ser", "f"), ci = TRUE, ci.level = clevel, 
                                  single.row = TRUE, title = paste(mfamily, 
                                                                   " regression (", model$family$link, 
                                                                   "): ", model$formula[2], " ", model$formula[1], 
                                                                   " ", model$formula[3], sep = ""))
                      }
                      if (mclass == "gam") {
                        Msum <- list(summary(model)$formula, summary(model)$p.table, 
                                     summary(model)$s.table)
                        invisible(lapply(1:3, function(i) {
                          cat(sep = "", names(Msum)[i], "\n")
                          print(Msum[[i]])
                        }))
                      }
                      if (mclass == "Gam") {
                        Msum <- list(model$formula, summary(model)$parametric.anova, 
                                     summary(model)$anova)
                        invisible(lapply(1:3, function(i) {
                          cat(sep = "", names(Msum)[i], "\n")
                          print(Msum[[i]])
                        }))
                      }
                      if (mclass %in% c("ols", "lrm", "Glm")) {
                        stargazer(model, type = "text", omit.stat = c("LL", 
                                                                      "ser", "f"), ci = TRUE, ci.level = clevel, 
                                  single.row = TRUE, title = paste("Linear Regression:", 
                                                                   model$call[2], sep = " "))
                      }
                    }
                    if (m.summary == "raw") {
                      if (mclass %in% c("ols", "Glm", "lrm")) {
                        print(model)
                      } else {
                        summary(model)
                      }
                    }
                  })
                }))
  }
  
  DNbuilder_czx_lrm<-function (model, data = NULL, clevel = 0.95, m.summary = c("raw","formatted"), 
                               covariate = c("slider", "numeric"), ptype = c("st","1-st"), 
                               DNtitle = NULL, DNxlab = NULL, DNylab = NULL, DNlimits = NULL,KMtitle = NULL, KMxlab = NULL, KMylab = NULL) 
  {
    mclass <- getclass.DN(model)$model.class
    mfamily <- getclass.DN(model)$model.family
    if (mclass %in% c("coxph", "cph")) {
      Surv.in <- length(model$terms[[2]]) != 1
    }
    if (mclass %in% c("ols", "Glm", "lrm", "cph")) {
      model <- update(model, x = T, y = T)
    }
    if (!is.data.frame(data)) {
      if (any(class(try(getdata.DN(model), silent = TRUE)) == 
              "try-error")) {
        stop("Dataset needs to be provided in a data.frame format")
      }
      else {
        data <- getdata.DN(model)
      }
    }
    covariate <- match.arg(covariate)
    m.summary <- match.arg(m.summary)
    ptype <- match.arg(ptype)
    if (mclass %in% c("lm", "glm", "ols", "Glm", "lrm", "gam", 
                      "Gam")) {
      Terms.T <- all(all.vars(model$terms) %in% names(data))
    }
    if (mclass %in% c("coxph")) {
      if (Surv.in) {
        Terms.T <- all(all.vars(model$terms)[-c(1:2)] %in% 
                         names(data))
      }
      else {
        Terms.T <- all(all.vars(model$terms)[-1] %in% names(data))
      }
    }
    if (mclass %in% c("cph")) {
      Terms.T <- all(names(model$Design$units) %in% names(data))
    }
    if (!Terms.T) 
      stop("Error in model syntax: some of model's terms do not match to variables' name in dataset")
    if (!is.null(DNlimits) & !length(DNlimits) == 2) 
      stop("A vector of 2 is required as 'DNlimits'")
    if (is.null(DNtitle)) 
      DNtitle <- "Dynamic Nomogram"
    if (is.null(DNxlab)) {
      if ((mclass %in% c("glm") & mfamily %in% c("binomial", 
                                                 "quasibinomial")) | mclass == "lrm") {
        DNxlab <- "Probability"
      }
      else {
        DNxlab <- ifelse(mclass %in% c("coxph", "cph"), 
                         "Survival probability", "Response variable")
      }
    }
    if (mclass %in% c("coxph", "cph")) {
      if (is.null(KMtitle)) {
        if (ptype == "st") {
          KMtitle <- "Estimated Survival Probability"
        }
        else {
          KMtitle <- "Estimated Probability"
        }
      }
      if (is.null(KMxlab)) {
        KMxlab <- "Follow Up Time"
      }
      if (is.null(KMylab)) {
        if (ptype == "st") {
          KMylab <- "S(t)"
        }
        else {
          KMylab <- "F(t)"
        }
      }
    }
    if (mclass %in% c("lm", "glm", "ols", "Glm", "lrm", "gam", 
                      "Gam")) {
      DNbuilder.core_czx_lrm(model, data, clevel, m.summary, covariate, 
                             DNtitle, DNxlab, DNylab, DNlimits)
    }
    if (mclass %in% c("coxph", "cph")) {
      DNbuilder.surv(model, data, clevel, m.summary, covariate, 
                     ptype, DNtitle, DNxlab, DNylab, KMtitle, KMxlab, 
                     KMylab)
    }
  }
  DNbuilder.core_czx_lrm<-function (model, data, clevel, m.summary, covariate, DNtitle, 
                                    DNxlab, DNylab, DNlimits) 
  {
    mclass <- getclass.DN(model)$model.class
    mfamily <- getclass.DN(model)$model.family
    if (mclass %in% c("lm", "ols")) {
      mlinkF <- function(eta) eta
    }
    else {
      mlinkF <- ifelse(mclass == "lrm", function(mu) plogis(mu), 
                       model$family$linkinv)
    }
    input.data <- NULL
    old.d <- NULL
    if (mclass %in% c("ols", "Glm", "lrm")) {
      model <- update(model, x = T, y = T)
    }
    allvars <- all.vars(model$terms)
    if (mclass %in% c("ols", "lrm", "Glm")) {
      
      #进行lrm的分析时，此处term原意应为因变量及自变量的名称的character形式，但是此处仅生成自变量的部分，缺少因变量。
      #因此debug此处，增加一项因变量名称在前
      terms <- model$Design$assume[model$Design$assume != 
                                     "interaction"]
      terms<-c("asis",terms)
      names(terms) = c(model$terms[[2]],model$Design$name[model$Design$assume != "interaction"])
      #结束debug
      
      if (mclass %in% c("Glm")) {
        terms <- c(attr(attr(model$model, "terms"), "dataClasses")[1], 
                   terms)
      }
      else {
        terms <- c(attr(model$terms, "dataClasses")[1], 
                   terms)
      }
    }
    if (mclass %in% c("lm", "glm", "gam", "Gam")) {
      if (length(attr(model$terms, "dataClasses")) == length(allvars)) {
        terms <- attr(model$terms, "dataClasses")
        names(terms) = allvars
      }
      else {
        terms <- attr(model$terms, "dataClasses")[which(names(attr(model$terms, 
                                                                   "dataClasses")) %in% allvars)]
      }
    }
    if (terms[[1]] == "logical") 
      stop("Error in model syntax: logical form for response not supported")
    terms[terms %in% c("numeric", "asis", "polynomial", "integer", 
                       "double", "matrx") | grepl("nmatrix", terms, fixed = T) | 
            grepl("spline", terms, fixed = T)] = "numeric"
    terms[terms %in% c("factor", "ordered", "logical", "category", 
                       "scored")] = "factor"
    resp <- terms[1]
    names(resp) <- allvars[1]
    if ("(weights)" %in% names(terms)) {
      preds <- as.list(terms[-c(1, length(terms))])
    }
    else {
      preds <- as.list(terms[-1])
    }
    names(preds) <- allvars[-1]
    for (i in 1:length(preds)) {
      if (preds[[i]] == "numeric") {
        i.dat <- which(names(preds[i]) == names(data))
        preds[[i]] <- list(v.min = floor(min(na.omit(data[, 
                                                          as.numeric(i.dat)]))), v.max = ceiling(max(na.omit(data[, 
                                                                                                                  as.numeric(i.dat)]))), v.mean = zapsmall(mean(data[, 
                                                                                                                                                                     as.numeric(i.dat)], na.rm = T), digits = 4))
        next
      }
      if (preds[[i]] == "factor") {
        i.dat <- which(names(preds[i]) == names(data))
        if (mclass %in% c("ols", "Glm", "lrm", "cph")) {
          preds[[i]] <- list(v.levels = model$Design$parms[[which(names(preds[i]) == 
                                                                    names(model$Design$parms))]])
        }
        else {
          preds[[i]] <- list(v.levels = model$xlevels[[which(names(preds[i]) == 
                                                               names(model$xlevels))]])
        }
      }
    }
    if (!is.null(DNlimits) & !length(DNlimits) == 2) 
      stop("A vector of 2 is required as 'DNlimits'")
    if (is.null(DNlimits)) {
      if ((mclass %in% c("glm") & mfamily %in% c("binomial", 
                                                 "quasibinomial")) | mclass == "lrm") {
        limits0 <- c(0, 1)
      }
      else {
        if (mclass %in% c("lm", "glm", "gam", "Gam")) {
          limits0 <- c(mean(model$model[, names(resp)]) - 
                         3 * sd(model$model[, names(resp)]), mean(model$model[, 
                                                                              names(resp)]) + 3 * sd(model$model[, names(resp)]))
        }
        if (mclass %in% c("ols", "lrm", "Glm")) {
          limits0 <- c(mean(model$y) - 3 * sd(model$y), 
                       mean(model$y) + 3 * sd(model$y))
        }
      }
      if (mclass %in% c("glm", "Glm") & mfamily %in% c("poisson", 
                                                       "quasipoisson", "Gamma")) {
        limits0[1] <- 0
      }
    }
    else {
      limits0 <- DNlimits
    }
    neededVar <- c(names(resp), names(preds))
    data <- data[, neededVar]
    input.data <- data[0, ]
    model <- update(model, data = data)
    wdir <- getwd()
    app.dir <- paste(wdir, "DynNomapp", sep = "/")
    message(paste("creating new directory: ", app.dir, sep = ""))
    dir.create(app.dir)
    setwd(app.dir)
    message(paste("Export dataset: ", app.dir, "/dataset.RData", 
                  sep = ""))
    save(data, model, preds, resp, mlinkF, getpred.DN, getclass.DN, 
         DNtitle, DNxlab, DNylab, DNlimits, limits0, terms, input.data, 
         file = "data.RData")
    message(paste("Export functions: ", app.dir, "/functions.R", 
                  sep = ""))
    dump(c("getpred.DN", "getclass.DN"), file = "functions.R")
    if (!is.null(DNlimits)) {
      limits.bl <- paste("limits <- reactive({ DNlimits })")
    }
    else {
      limits.bl <- paste("limits <- reactive({ if (input$limits) { limits <- c(input$lxlim, input$uxlim) } else {\n                         limits <- limits0 } })")
    }
    noSE.bl <- paste("by '", mclass, "'", sep = "")
    p1title.bl <- paste(clevel * 100, "% ", "Confidence Interval for Response", 
                        sep = "")
    p1msg.bl <- paste("Confidence interval is not available as there is no standard errors available by '", 
                      mclass, "' ", sep = "")
    if (m.summary == "formatted") {
      if (mclass == "lm") {
        sumtitle.bl <- paste("Linear Regression:", model$call[2], 
                             sep = " ")
      }
      if (mclass %in% c("glm")) {
        sumtitle.bl <- paste(mfamily, " regression (", model$family$link, 
                             "): ", model$formula[2], " ", model$formula[1], 
                             " ", model$formula[3], sep = "")
      }
      if (mclass %in% c("ols", "lrm", "Glm")) {
        sumtitle.bl <- paste("Linear Regression:", model$call[2], 
                             sep = " ")
      }
    }
    else {
      sumtitle.bl = NULL
    }
    if (m.summary == "formatted") {
      if (mclass %in% c("lm", "glm", "ols", "lrm", "Glm")) {
        sum.bi <- paste("stargazer(model, type = 'text', omit.stat = c('LL', 'ser', 'f'), ci = TRUE, ci.level = clevel, single.row = TRUE, title = '", 
                        sumtitle.bl, "')", sep = "")
      }
      if (mclass == "gam") {
        sum.bi <- paste("Msum <- list(summary(model)$formula, summary(model)$p.table, summary(model)$s.table)\n                invisible(lapply(1:3, function(i){ cat(sep='', names(Msum)[i], '\n') ; print(Msum[[i]])}))")
      }
      if (mclass == "Gam") {
        sum.bi <- paste("Msum <- list(model$formula, summary(model)$parametric.anova, summary(model)$anova)\n                invisible(lapply(1:3, function(i){ cat(sep='', names(Msum)[i], '\n')) ; print(Msum[[i]])}))")
      }
    }
    if (m.summary == "raw") {
      if (mclass %in% c("ols", "Glm", "lrm")) {
        sum.bi <- paste("print(model)")
      }
      else {
        sum.bi <- paste("summary(model)")
      }
    }
    if (mclass %in% c("ols", "Glm", "lrm")) {
      datadist.bl <- paste("t.dist <- datadist(data)\noptions(datadist = 't.dist')", 
                           sep = "")
    }
    else {
      datadist.bl <- ""
    }
    if (mclass %in% c("lm", "glm")) {
      library.bl <- ""
    }
    else {
      if (mclass %in% c("ols", "Glm", "lrm")) {
        library.bl <- paste("library(rms)")
      }
      if (mclass %in% c("Gam")) {
        library.bl <- paste("library(gam)")
      }
      if (mclass %in% c("gam")) {
        library.bl <- paste("library(mgcv)")
      }
    }
    GLOBAL = paste("library(ggplot2)\nlibrary(shiny)\nlibrary(plotly)\nlibrary(stargazer)\nlibrary(compare)\nlibrary(prediction)\n", 
                   library.bl, "\n\n#######################################################\n#### Before publishing your dynamic nomogram:\n####\n#### - You may need to edit the following lines if\n#### data or model objects are not defined correctly\n####
                   - You could modify ui.R or server.R for\n#### making any required changes to your app\n#######################################################\n\nload('data.RData')\nsource('functions.R')\n", 
                   datadist.bl, "\nm.summary <- '", m.summary, "'\ncovariate <- '", 
                   covariate, "'\nclevel <- ", clevel, "\n\n### Please cite the package if used in publication. Use:\n# Amirhossein Jalali, Davood Roshan, Alberto Alvarez-Iglesias and John Newell (2019). DynNom: Visualising statistical models using dynamic
                   nomograms.\n# R package version 5.0. https://CRAN.R-project.org/package=DynNom\n", 
                   sep = "")
    UI = paste("ui = bootstrapPage(fluidPage(\n    titlePanel('", 
               DNtitle, "'),\n    sidebarLayout(sidebarPanel(uiOutput('manySliders'),\n                               uiOutput('setlimits'),\n                               actionButton('add', 'Predict'),\n                               br(), br(),\n      
               helpText('Press Quit to exit the application'),\n                               actionButton('quit', 'Quit')\n    ),\n    mainPanel(tabsetPanel(id = 'tabs',\n 
               tabPanel('Graphical Summary', plotlyOutput('plot')),\n                          tabPanel('Numerical Summary', verbatimTextOutput('data.pred')),\n                          tabPanel('Model Summary', verbatimTextOutput('summary'))\n    )\n    )\n    )))", 
               sep = "")
    SERVER = paste("server = function(input, output){\nobserve({if (input$quit == 1)\n          stopApp()})\n\n", 
                   limits.bl, "\n\noutput$manySliders <- renderUI({\n  slide.bars <- list()\n               for (j in 1:length(preds)){\n               if (terms[j+1] == \"factor\"){\n               slide.bars[[j]] <- list(selectInput(paste(\"pred\", j, sep = \"\"),
                   names(preds)[j], preds[[j]]$v.levels, multiple = FALSE))\n               }\n               if (terms[j+1] == \"numeric\"){\n               if (covariate == \"slider\") {\n               
                   slide.bars[[j]] <- list(sliderInput(paste(\"pred\", j, sep = \"\"), names(preds)[j],\n               min = preds[[j]]$v.min, max = preds[[j]]$v.max, value = preds[[j]]$v.mean))\n               }\n               if (covariate == \"numeric\") {\n 
                   slide.bars[[j]] <- list(numericInput(paste(\"pred\", j, sep = \"\"), names(preds)[j], value = zapsmall(preds[[j]]$v.mean, digits = 4)))\n               }}}\n   
                   do.call(tagList, slide.bars)\n})\n\noutput$setlimits <- renderUI({\n        if (is.null(DNlimits)){\n               setlim <- list(checkboxInput(\"limits\", \"Set x-axis ranges\"),\n               conditionalPanel(condition = \"input.limits == true\",\n               
                   numericInput(\"uxlim\", \"x-axis upper\", zapsmall(limits0[2], digits = 2)),\n               numericInput(\"lxlim\", \"x-axis lower\", zapsmall(limits0[1], digits = 2))))\n        } else{ setlim <- NULL }\n        setlim\n})\n\na <- 0\nnew.d <- reactive({\n
                   input$add\n               input.v <- vector(\"list\", length(preds))\n               for (i in 1:length(preds)) {\n               input.v[[i]] <- isolate({\n
                   input[[paste(\"pred\", i, sep = \"\")]]\n               })\n               names(input.v)[i] <- names(preds)[i]\n               }\n               out <- data.frame(lapply(input.v, cbind))\n               if (a == 0) {\n               
                   input.data <<- rbind(input.data, out)\n               }\n               if (a > 0) {\n               if (!isTRUE(compare(old.d, out))) {\n               input.data <<- rbind(input.data, out)\n               }}\n               a <<- a + 1\n               
                   out\n})\n\np1 <- NULL\nold.d <- NULL\ndata2 <- reactive({\n               if (input$add == 0)\n               return(NULL)\n               if (input$add > 0) {\n              
                   if (!isTRUE(compare(old.d, new.d()))) {\n               isolate({\n               mpred <- getpred.DN(model, new.d(), set.rms=T)$pred\n               se.pred <- getpred.DN(model, new.d(), set.rms=T)$SEpred\n               if (is.na(se.pred)) {\n               
                   lwb <- \"No standard errors\"\n               upb <- \"", 
                   noSE.bl, "\"\n               pred <- mlinkF(mpred)\n               d.p <- data.frame(Prediction = zapsmall(pred, digits = 3),\n               Lower.bound = lwb, Upper.bound = upb)\n               } else {\n               lwb <- sort(mlinkF(mpred + 
                   cbind(1, -1) * (qnorm(1 - (1 - clevel)/2) * se.pred)))[1]\n               upb <- sort(mlinkF(mpred + cbind(1, -1) * (qnorm(1 - (1 - clevel)/2) * se.pred)))[2]\n              
                   pred <- mlinkF(mpred)\n               d.p <- data.frame(Prediction = zapsmall(pred, digits = 3),\n               Lower.bound = zapsmall(lwb, digits = 3),\n               Upper.bound = zapsmall(upb, digits = 3))\n               }\n               old.d <<- new.d()\n               
                   data.p <- cbind(d.p, counter = 1, count=0)\n               p1 <<- rbind(p1, data.p)\n               p1$counter <- seq(1, dim(p1)[1])\n               p1$count <- 0:(dim(p1)[1]-1) %% 11 + 1\n               p1\n               })\n               } else {\n
                   p1$count <- seq(1, dim(p1)[1])\n               }}\n               rownames(p1) <- c()\n               p1\n})\n\noutput$plot <- renderPlotly({\n  if (input$add == 0)\n
                   return(NULL)\n               if (is.null(new.d()))\n               return(NULL)\n               coll=c(\"#0E0000\", \"#0066CC\", \"#E41A1C\", \"#54A552\", \"#FF8000\", \"#BA55D3\",\n               \"#006400\", \"#994C00\", \"#F781BF\", \"#00BFFF\", \"#A9A9A9\")\n               
                   lim <- limits()\n               yli <- c(0 - 0.5, 10 + 0.5)\n               dat2 <- data2()\n               if (dim(data2())[1] > 11){\n               input.data = input.data[-c(1:(dim(input.data)[1]-11)),]\n               dat2 <- data2()[-c(1:(dim(data2())[1]-11)),]\n
                   yli <- c(dim(data2())[1] - 11.5, dim(data2())[1] - 0.5)\n               }\n               in.d <- input.data\n               
                   xx <- matrix(paste(names(in.d), \": \", t(in.d), sep = \"\"), ncol = dim(in.d)[1])\n               Covariates <- apply(xx, 2, paste, collapse = \"<br />\")\n               p <- ggplot(data = dat2, aes(x = Prediction, y = counter - 1, text = Covariates,\n
                   label = Prediction, label2 = Lower.bound, label3=Upper.bound)) +\n
                   geom_point(size = 2, colour = coll[dat2$count], shape = 15) +\n               ylim(yli[1], yli[2]) + coord_cartesian(xlim = lim) +\n               labs(title = \"", 
                   p1title.bl, "\",\n               x = \"", DNxlab, "\", y = \"", 
                   DNylab, "\") + theme_bw() +\n               theme(axis.text.y = element_blank(), text = element_text(face = \"bold\", size = 10))\n               if (is.numeric(dat2$Upper.bound)){\n               p <- p + geom_errorbarh(xmax = dat2$Upper.bound, xmin = dat2$Lower.bound,\n
                   size = 1.45, height = 0.4, colour = coll[dat2$count])\n               } else{\n               message(\"", 
                   p1msg.bl, "\")\n               }\n               gp <- ggplotly(p, tooltip = c(\"text\", \"label\", \"label2\", \"label3\"))\n               gp$elementId <- NULL\n               gp\n})\n\noutput$data.pred <- renderPrint({\n  if (input$add > 0) {\n               
                   if (nrow(data2()) > 0) {\n               if (dim(input.data)[2] == 1) {\n               in.d <- data.frame(input.data)\n               names(in.d) <- names(terms)[2]\n  
                   data.p <- cbind(in.d, data2()[1:3])\n               }\n               if (dim(input.data)[2] > 1) {\n               data.p <- cbind(input.data, data2()[1:3])\n               }}\n               stargazer(data.p, summary = FALSE, type = \"text\")\n}\n})\n\noutput$summary 
                   <- renderPrint({\n",                    sum.bi, "\n})\n}", sep = "")
    output = list(ui = UI, server = SERVER, global = GLOBAL)
    text <- paste("This guide will describe how to deploy a shiny application using scripts generated by DNbuilder:\n\n1. Run the shiny app by setting your working directory to the DynNomapp folder, and then run: shiny::runApp() If you are using the RStudio IDE, you can also run it by 
    clicking the Run App button in the editor toolbar after open one of the R scripts.\n\n2. You could modify codes to apply all the necessary changes. Run again to 
    confirm that your application works perfectly.\n\n3. Deploy the application by either clicking on the Publish button in the top right corner of the running app, or use the generated files and deploy it on your server if you host any.\n\nYou can find a full guide of how to deploy an 
    application on shinyapp.io server 
                  here:\nhttp://docs.rstudio.com/shinyapps.io/getting-started.html#deploying-applications\n\nPlease cite the package if using in publication.", 
                  sep = "")
    message(paste("writing file: ", app.dir, "/README.txt", 
                  sep = ""))
    writeLines(text, "README.txt")
    message(paste("writing file: ", app.dir, "/ui.R", sep = ""))
    writeLines(output$ui, "ui.R")
    message(paste("writing file: ", app.dir, "/server.R", sep = ""))
    writeLines(output$server, "server.R")
    message(paste("writing file: ", app.dir, "/global.R", sep = ""))
    writeLines(output$global, "global.R")
    setwd(wdir)
  }
}

##选择训练集数据文件
#logistic回归数据格式
#数据文件格式要求：注意数据文件要保存为txt版本
#数据格式要求：1. 数据中不要有空格；
#              2. 数据中第一行为变量名；
#              3. 数据中第一列为预测结局二分类指标（结局指标必须设置为0/1形式）；
#              4. 数据中从第二列开始为其他变量；
#              6. 空值以NA填充；二分类变量以原始注释表示；多分类变量不进行哑变量处理，直接以原始注释表示；
#模板为生成的仿真数据
#cox 回归数据格式
#数据文件格式要求：注意数据文件要保存为txt版本
#数据格式要求：1. 数据中不要有空格；
#              2. 数据中第一行为变量名；
#              3. 数据中前两列为预测的时间数据和截尾指标（截尾指标必须设置为0/1形式）；
#              4. 数据中从第二列开始为其他变量；
#              6. 空值以NA填充；二分类变量以原始注释表示；多分类变量 不 进行哑变量处理，直接以原始注释表示；
#模板为生成的仿真数据
#选择文件路径
library(rms)
library(Formula)
library(DynNom)
library(shiny)
library(plotly)
library(compare)
library(stargazer)
rm(list=ls())
#training_dataset_path <- choose.files(caption = "R",
#                                      multi = TRUE, filters = Filters,
#                                     index = nrow(Filters))


#training_dataset<- read.csv(training_dataset_path, header = TRUE,sep="\t", stringsAsFactors=FALSE)
#导入文件#建议手动
training_dataset<-R
#二分类变量赋值
training_dataset$DMIr <- factor(training_dataset$DMIr,levels = c(0,1),labels = c("No","Yes"))
#查看数据
print(paste0("该训练集有",dim(training_dataset)[1]," 个样本；",dim(training_dataset)[2]," 个变量"))
##1.0 构建多因素回归模型
ddist <- datadist(training_dataset)
options(datadist='ddist')
#查看数据集中的变量名
colnames(training_dataset)
#选择要纳入多因素回归的变量，并修改此处函数(用+号连接自变量)
f_lrm <- lrm(Label~SIRI+D.+Radscore, data=training_dataset, x=TRUE, y=TRUE,maxit=5000)
#查看多因素Logistic分析结果，最下方可见其对应的Coef、p值
print(f_lrm)

##2.0 尝试构建本地页面列线图
#原装代码如下
#Logistic
DynNom(f_lrm,training_dataset,
       clevel = 0.95) #生成动态列线图
#Cox
DynNom(f_cph,training_dataset,
       clevel = 0.95) #生成动态列线图
#如对lrm对象处理时出现如下报错：
#'names'属性的长度[6]必需和矢量的长度[5]一样
#可尝试改用我们修订后的代替函数DynNom_czx_lrm如下
DynNom_czx_lrm(model = f_lrm,data = training_dataset,
               clevel = 0.95) #生成动态列线图
#关掉页面后跳出来的debug页面不用管，直接"STOP"按钮关掉就好


#3.0 导出本地DynNomapp脚本文件
##0.3 设置导出的文件夹目录，注意完整地址中不可以有中文，否则后面无法正常上传
#注意文件夹路径的斜杠是双斜杠\\
output_dir <- choose.dir(default = "D:\\", caption = "D盘")
setwd(output_dir)
DNbuilder_czx_lrm(model = f_lrm,data = training_dataset,
                  clevel = 0.95)





####—训练+验证—####
library(devtools)
library(openxlsx)
library(survival)
library(survminer)
library(lattice)
library(Formula)
library(ggplot2)
library(Hmisc)
library(SparseM)
library(foreign)
library(rms)
library(caret)
library(car)
library(pROC)
library(timeROC)
library(ggDCA)
library(dcurves)
library(dplyr)
library(Rcpp)

####单因素Cox回归分析####
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
training_dataset <- data
#training_dataset$time=training_dataset$time*30

colnames(training_dataset)
ddist <- datadist(training_dataset)
options(datadist='ddist')
f_cph <- cph(Surv(time,status) ~ METTL3,
             x=T, y=T, surv=T, time.inc=3*12,
             data=training_dataset)
#查看单因素Cox分析结果，最下方可见其对应的Coef、p值
print(f_cph)




####——KM——####
data<-read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
fit3 <- survfit(Surv(time,status) ~FKBP9, data = data)
print(fit3)
plot(fit3)

data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
#究极进化
fit<- survfit(Surv(time, status) ~ FKBP9, data = data)
ggsurvplot(fit, data = data) #得到基础的生存曲线啦
res_cox<-coxph(Surv(time,status) ~ FKBP9, data=data)
p5 <- ggsurvplot(fit, data = data,
                 legend.title="Treatments", 
                 palette=c("#007fff","#FF0000"), #0,1"#007fff","#FF0000"
                 ylab="PFS (%)",xlab = "Months", #更改横纵坐标
                 conf.int = TRUE, #给生存曲线添加上置信区间
                 legend.labs = c("CCRT","CCRT&IO"), #在图中添加图例
                 legend.size=10,
                 #legend = c(0.90,0.90), # 调整图例位置
                 risk.table = TRUE,  #添加上风险表
                 break.x.by = 20)  #调整x轴刻度间距
p5$plot = p5$plot + ggplot2::annotate("text",x = 12, y = 0.17,label = paste("HR:",round(summary(res_cox)$conf.int[1],3))) + 
  ggplot2::annotate("text",x = 12, y = 0.10,label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],3),"-",round(summary(res_cox)$conf.int[4],3),")",sep = ""))+
  ggplot2::annotate("text",x = 12, y = 0.03,label = paste("P:",round(summary(res_cox)$coef[5],6)))
p5


#计算1年、3年、5年HR及CI####
library(ggsurvfit)
library(ggsurvfit)
data<-read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE)
cox_fit <- coxph(Surv(time,status) ~ FKBP9, data=data)
# 提取HR值和置信区间
hr <- round(summary(cox_fit)$conf.int[1],3)
hr_lower <- round(summary(cox_fit)$conf.int[3], 3)
hr_upper <- round(summary(cox_fit)$conf.int[4], 3)

# 提取N年生存率
fit <- survfit(Surv(time, status) ~ FKBP9, data = data)
N=12 #时间
surv_CTRL <- round(summary(fit, times = c(N))$surv[1], 3)*100
surv_CTRL_lower <- round(summary(fit, times = c(N))$lower[1], 3)*100
surv_CTRL_upper <- round(summary(fit, times = c(N))$upper[1], 3)*100
surv_Treatment <- round(summary(fit, times = c(N))$surv[2], 3)*100
surv_Treatment_lower <- round(summary(fit, times = c(N))$lower[2], 3)*100
surv_Treatment_upper <- round(summary(fit, times = c(N))$upper[2], 3)*100

#Events and median PFS
fit
df <- data.frame(Group = c("Placebo plus\nchemotherapy (138)", "Tislelizumab\nchemotherapy (90)"),
                 Events = c("112 (81%)", "53 (59%)"),
                 `Median OS` = c("8.87 (95% CI 6.97–11.43)", "14 (95% CI 10.2–18.1)"))
# 绘制生存曲线
p <- survfit2(Surv(time, status) ~ FKBP9, data = data ) %>%  ggsurvfit(linewidth = 0.8) +
              #add_risktable(risktable_height = 0.33, #风险表格
              #risktable_stats = c("{n.risk}\n({cum.censor})"), 
              #stats_label = list(n.risk = "Number at risk", cum.censor = "number censored"),
              #size = 5,theme =   # increase font size of risk table title and y-axis label
              #list(theme_risktable_default(axis.text.y.size = 12,plot.title.size = 12),
              #theme(plot.title = element_text(face = "bold"))) ) + 
  #add_risktable_strata_symbol(symbol = "\U25CF", size = 20) +
  #add_censor_mark(shape = 1, size = 2, stroke = 1) +  
  add_pvalue(location  = "annotation", x = 40, y = 0.70 , hjust = 1,size = 4, caption = "Stratified log-rank {p.value}") + # 添加log-rank p-value 
  annotate("text",    # 添加HR和置信区间
           x = 40,y = 0.75, label = paste("Hazard ratio,", round(hr,2), 
           "(95% CI,", round(hr_lower,2), "-", round(hr_upper,2), ")", sep = ""),
           size = 4, hjust = 1  ) +  
  geom_segment(mapping = aes(x = 12, y = 0, xend = 12, yend = 0.526), #X代表哪一年
               linetype = "dashed", linewidth = 0.4) +
  annotate("text",    # 添加半年生存率和置信区间-CTRL           
           x = 20,y = 0.28,label = paste(surv_CTRL,"%","\n","(95% CI ",surv_CTRL_lower, "-", surv_CTRL_upper, ")", sep = ""), 
           size = 4, hjust = 1  ) +  annotate("text", # 添加半年生存率和置信区间-Treatmet 
           x = 10,y = 0.60, label = paste(surv_Treatment,"%","\n","(95% CI ", surv_Treatment_lower, "-", surv_Treatment_upper, ")", sep = ""), 
           size = 4,hjust = 0  ) +  
  labs(title = "", x = "Time from randomisation (months)", 
                   y = "Overall survival (%)"  
       ) + scale_x_continuous(breaks = seq(0, 39, 3), expand = c(0.03, 0)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = seq(0, 100, 10)) + 
  scale_color_manual(values = c('#009FC3', '#B30437')) + 
  scale_fill_manual(values = c('#009FC3', '#B30437')) + 
  theme_classic() +  theme(axis.text = element_text(size = 14,color = "black"),
                           axis.title.y = element_text(size = 14,color = "black", vjust = 0.5),
                           axis.title.x = element_text(size = 14,color = "black"),
                           panel.grid =  element_blank(),
                           plot.margin = margin(r = 0.1, l = 0.1, unit = "cm"),
                           legend.text = element_text(size = 12,color = "black"), 
                           legend.background = element_blank(), 
                           legend.position = c(0.17,0.1) )
p
#添加注释表格：分组，事件和中位生存值
library(ggpp)
library(ggplot2)
library(tibble)
p + geom_table(data = df, aes(x = 33, y = 1, label = list(df)),size = 4,family = "sans",
               table.theme = ttheme_gtsimple, #ttheme_gtminimal
               table.hjust = 0, color = "black", fill = "transparent") 






####批量单因素COX####
rm(list=ls())
set.seed(123)#
data<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE,row.names = 1,)
credit_df <- data
index<-sample(1:nrow(credit_df),nrow(credit_df)*0.67,replace=F)#0.7代表37分
training_dataset<-credit_df[index,]
test_dataset<-credit_df[-index,]
dim(training_dataset)#训练集training_dataset
dim(test_dataset)#验证集test_dataset

R<- data[1:75,]
R<- data[76:150,]
R<- data[151:225,]
R<- data[226:300,]

R<- data[376:450,]
R<- data[451:525,]
R<- data[536:600,]

R<- data[751:825,]
R<- data[826:900,]

R<- data[1126:1200,]




can <- colnames(R)[3:38]#批量计算那些列
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

write.table(result, file = "COX单因素分析结果CD68-CD8.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)


####多因素Cox回归分析####
training_dataset<- read.csv("C:/Users/Administrator/Desktop/R.csv",header = TRUE,row.names = 1)
ddist <- datadist(training_dataset)
options(datadist='ddist')
#下面改公式
colnames(training_dataset)
f_cph <- cph(Surv(time,status) ~cN,
             x=T, y=T, surv=T, time.inc=3*12,
             data=training_dataset)
#查看多因素Cox分析结果，最下方可见其对应的Coef、p值
print(f_cph)

med  <- Quantile(f_cph)
surv <- Survival(f_cph)
#列线图####
plot(nomogram(f_cph, fun=function(x) med(lp=x), funlabel="Median Survival Time"))
#绘制1年、3年、5年生存时间的列线图。
nom<-nomogram(f_cph, fun=list(function(x) surv(12, x),
                              function(x) surv(3*12, x),
                              function(x) surv(5*12, x)),
              funlabel=c("1-year OS", 
                         "3-year OS",
                         "5-year OS"))
plot(nomogram(f_cph, fun=list(function(x) surv(12, x),
                              function(x) surv(3*12, x),
                              function(x) surv(5*12, x)),
              funlabel=c("1-year OS", 
                         "3-year OS",
                         "5-year OS"))
)

#计算总分
ddist <- datadist(training_dataset)
option<-options(datadist='ddist')
options(option)
library(nomogramFormula)
#计算总分，并作为新变量添加进原始数据
result1<-formula_rd(nomogram=nom)
training_dataset$points<-points_cal(formula=result1$formula,rd=training_dataset)
test_dataset$points<-points_cal(formula=result1$formula,rd=test_dataset)
head(training_dataset)
head(test_dataset)

write.table(training_dataset, file = "生存训练组列线图参数.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(test_dataset, file = "生存验证组列线图参数.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)



####森林图####
library(survminer)
library(forplo)
library(survival)
#黑白森林图
colnames(training_dataset)
multiCox <- coxph(Surv(time,status) ~Age+M+TNM+Pathscore,
                  data=training_dataset)
multiCox
multicox_ggforest <- ggforest(multiCox,  #coxph得到的Cox回归结果
                              data = training_dataset,  #数据集
                              main = 'Hazard ratio of multi cox',  #标题
                              cpositions = c(0.05, 0.15, 0.35),  #前三列距离
                              fontsize = 0.8, #字体大小
                              refLabel = 'reference', #相对变量的数值标签，也可改为1
                              noDigits = 3 #保留HR值以及95%CI的小数位数
)
multicox_ggforest


#简单森林图
#logistics回归
#mod1<-glm(status~rx+sex+age+obstruct+nodes,data=colon,family="binomial")
#Cox回归
mydata<- read.csv("D:/Desktop/R.csv",header = TRUE)
ddist <- datadist(mydata)
options(datadist='ddist')
colnames(mydata)

mydata$IOD<-factor(mydata$IOD,levels=c(0,1),labels = c("Low","High"))
#mydata$AFP<-factor(mydata$LI,levels=c(0,1),labels=c("negative","positive"))
mydata$LI<-factor(mydata$LI,levels=c(0,1),labels=c("No","Yes"))
mydata$DM<-factor(mydata$DM,levels=c(0,1),labels=c("No","Yes"))
mydata$LVSI<-factor(mydata$LVSI,levels=c(0,1),labels=c("No","Yes"))
mydata$LNM<-factor(mydata$LNM,levels=c(0,1),labels=c("No","Yes"))
mydata$Differentiation<-factor(mydata$Differentiation,levels=c(1,2,3),labels = c("Low","middle","High"))
mydata$TNM<-factor(mydata$TNM,levels=c(1,2,3,4),labels=c("I","II","III","IV"))

mul_cox <- coxph(Surv(time,status) ~Age+SCC+ADC+Radscore,
                 data=mydata)
mul_cox
library(forestmodel)
forest_model(mul_cox)



#复杂彩色森林图
#1.载入包
library(tableone)  
library(survival)
library(forestplot)
library(stringr)
#3.读入数据
aa<- read.csv("D:/Desktop/R.csv",header = TRUE)
ddist <- datadist(aa)
options(datadist='ddist')
mul_cox <- coxph(Surv(time,status) ~Age+SCC+ADC+Radscore,
                data=aa)
mul_cox
#一-2 multi1：用于提取：变量+HR+95%CI+95%CI，4列
mul_cox1 <- summary(mul_cox)
colnames(mul_cox1$conf.int)
multi1<-as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))
#一-3、multi2：用于提取：HR(95%CI)和p，2列
multi2<-ShowRegTable(mul_cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#一-4.将两次提取结果合并成表；取名result
result <-cbind(multi1,multi2);result
#一-5.行名转为表格第一列，并给予命名"Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result

fig1<- forestplot(result[,c(1,5,6)], #告诉函数，合成的表格result的第1，5，6列还是显示数字
                  mean=result[,2],   #告诉函数，表格第2列为HR，它要变成森林图的小方块
                  lower=result[,3],  #告诉函数表格第3列为5%CI，
                  upper=result[,4],  #表格第5列为95%CI，它俩要化作线段，穿过方块
                  zero=1,            #告诉函数，零线或参考线为HR=1即x轴的垂直线
                  boxsize=0.2,       #设置小黑块的大小
                  graph.pos=2)       #森林图应插在图形第2列
plot(fig1)

#3.加水平线、垂直辅助线、x轴标签、大标题、森林图占比
fig2<- forestplot(result[,c(1,5,6)], #合成的表格result的第1，5，6列还是显示数字
             mean=result[,2],   #告诉函数，表格第2列为HR，它要变成森林图的小方块
             lower=result[,3],  #告诉函数表格第3列为5%CI，
             upper=result[,4],  #表格第5列为95%CI，它俩要化作线段，穿过方块
             zero=1,            #告诉函数，零线或参考线为HR=1即x轴的垂直线
             boxsize=0.3,       #设置小黑块的大小
             graph.pos=2,       #"left"or "right"
             # ---------------#森林图的占比  
             graphwidth = unit(.3,"npc"),
             #----------------#X轴的垂直线
             #xlab="一只勤奋的科研喵",
             xticks=c(0,1,2,3,4,5,6) ,
             txt_gp=fpTxtGp(
               label=gpar(cex=1),
               ticks=gpar(cex=1.1), 
               xlab=gpar(cex=1.3), 
               title=gpar(cex=2)),
             #----------------#线条粗细（x轴、置信区间）
             lwd.zero=1,
             lwd.ci=1.5,
             lwd.xaxis=1.5, 
             lty.ci=1.5,
             ci.vertices =T,
             ci.vertices.height=0.2, 
             clip=c(0.1,8),
             #-----------------颜色
             col=fpColors(box = '#0d8931', #方块颜色，绿
                          lines = '#e5cc06', # 95CI线，黄
                          zero = "#00d5ff",# HR=1,即垂直线，浅蓝
                          text = "#021eaa",#文字，蓝
                          axes = "#f31017",#x轴，红
                          hrz_lines = "#de09dc"),#水平线，紫
             #其他方块变菱形或圆圈、页边距调整等
             fn.ci_norm="fpDrawCircleCI", #fpDrawNormalCI是方块
             #fpDrawDiamondCI是菱形、fpDrawCircleCI是圆形、fpDrawPointCI是圆圈
             mar=unit(rep(1.25, times = 4), "cm"), #图形页边距
             title="多因素Cox回归森林图" ,
             #---------------#添加水平线(可以不要)
             #hrzl_lines=list("1" = gpar(lty=1,lwd=1),#lwd=线粗
             #                "6" = gpar(lty=1))
)
fig2

#或者用这里的result在grampad prism上画森林图
write.table(result, file = "森林图原始数据2.xlsx",sep = "\t",row.names = T,col.names = NA,quote = F)





####ROC####
ROC <- timeROC(T = training_dataset$time,   
               delta = training_dataset$status,   
               marker = training_dataset$LVSI, # 这里就是risk score
               cause = 1,           
               weighting = "marginal",   
               times = c(12,36,60),       
               iid = TRUE)
ROC
plot(ROC, 
     time=12, col="red", lwd=2, title = "")  
plot(ROC,
     time=36, col="blue", add=TRUE, lwd=2)
plot(ROC,
     time=60, col="orange", add=TRUE, lwd=2)
legend("bottomright",c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],4)), 
                       paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],4)),
                       paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],4))),
       col=c("red", "blue", "orange"),
       lty=1, lwd=2,bty = "n") 

#验证集
ROC <- timeROC(T = test_dataset$time,   
               delta = test_dataset$status,   
               marker = test_dataset$points, # 这里就是risk score
               cause = 1,           
               weighting = "marginal",   
               times = c(12,36,60),       
               iid = TRUE)
ROC
plot(ROC, 
     time=12, col="red", lwd=2, title = "")  
plot(ROC,
     time=36, col="blue", add=TRUE, lwd=2)
plot(ROC,
     time=60, col="orange", add=TRUE, lwd=2)
legend("bottomright",c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],4)), 
                       paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],4)),
                       paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],4))),
       col=c("red", "blue", "orange"),
       lty=1, lwd=2,bty = "n") 



#c_index
f<-coxph(Surv(time,status==1)~Age+SCC+LNM+EQD2+Radscore, 
         data=training_dataset)
sum.surv<-summary(f)
c_index<-sum.surv$concordance
c_index

####校准图####
##校准曲线####
coxm<- cph(Surv(time,status) ~ Age+ADC+EBRT+Radscore3,x=T, y=T, surv=T,
              data=training_dataset)

cal1<-calibrate(coxm,cmethod='KM',method='boot',u=12,m=12,B=200)
cal2<-calibrate(coxm,cmethod='KM',method='boot',u=22,m=18,B=200)
cal3<-calibrate(coxm,cmethod='KM',method='boot',u=36,m=12,B=200)

plot(cal1,lwd=0,lty=0,errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
     xlim = c(0,1),ylim =c(0,1),
     xlab="Nomogram-Predicted Probability of DFS",
     ylab="Actual DFS (%)",
     col=c(rgb(0,118,192,maxColorValue = 255)))
plot(cal2,lwd=0,lty=0,errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
     xlim = c(0,1),ylim =c(0,1),
     xlab="Nomogram-Predicted Probability of DFS",
     ylab="Actual DFS (%)",
     col=c(rgb(0,118,192,maxColorValue = 255)))
plot(cal3,lwd=0,lty=0,##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlim = c(0,1),ylim =c(0,1),##x轴和y轴范围
     mgp = c(2, 1, 0), #控制坐标轴的位置
     xlab="Nomogram-Predicted Probability of DFS",
     ylab="Actual DFS (%)",
     col=c(rgb(0,118,192,maxColorValue = 255)))


lines(cal1[,c("mean.predicted","KM")],type="b",lwd=2,
      col=c(rgb(205,0,0,maxColorValue = 255)),pch=16)
lines(cal2[,c("mean.predicted","KM")],type="b",lwd=2,
      col=c(rgb(0,255,0,maxColorValue = 255)),pch=16)
lines(cal3[,c("mean.predicted","KM")],type="b",lwd=2,
      col=c(rgb(0,0,238,maxColorValue=255)),pch=16)

abline(0,1,lty=3,lwd=2,col=c(rgb(0,118,192,maxColorValue = 255)))




####决策曲线####
library(rmda)
library(tidyverse)
library(survival)
dca(Surv(time,status) ~pr_4_5, #4因子模型
    data=train,
    time=5,#代表几年
    thresholds = 1:100/100) %>%  #在1%到100%的概率阈值范围内计算决策曲线
  plot(smooth=T)
#旧，错误的
complex1<- decision_curve(status~SCC,data = training_dataset,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex2<- decision_curve(status~EBRT,data = training_dataset,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex3<- decision_curve(status~EQD2,data = training_dataset,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex4<- decision_curve(status~Radscore,data = training_dataset,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
complex5<- decision_curve(status~SCC+LNM+Age+Radscore1,data = training_dataset,
                          family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,study.design = 'case-control',
                          population.prevalence= 0.3)
List2<- list(complex1,complex2,complex3,complex4,complex5)
par(mai=c(1,1,0.5,1))
plot_decision_curve(List2,curve.names= c('SCC','ADC','EQD2','Radscore','Model'),
                    cost.benefit.axis =FALSE,
                    confidence.intervals =FALSE,standardize = FALSE) 
#如果不需要图例，可以加上legend.position = "none"。如下：
plot_decision_curve(List2,curve.names= c('1','2','3','4','5'),
                    cost.benefit.axis =FALSE,
                    confidence.intervals =FALSE,standardize = FALSE,legend.position = "none")

#临床影像曲线
plot_clinical_impact(complex5,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits= 8,col = c('red','blue'),
                     confidence.intervals= T)




####P次K折验####
#（第一次使用时才需要运行此部分代码。安装成功后不需要重复安装）
install.packages('car')
install.packages('survival')
install.packages('pROC')
install.packages('rms')
install.packages('tcltk')

rm(list = ls())
library(car)
library(survival)
library(pROC)
library(rms)
library(tcltk)
#下方为预设函数，请直接运行加???
#下方为内部验证相关预设函???,请直接运行加???
{
  nomogram_cross_validation<-function(f,dataset,N_fold=10,Iterations=200,seed=NULL){
    
    #设置随机数种???
    if(!is.null(seed)){
      set.seed(seed)
    }
    
    #提取公式中的因变量和自变???
    f<-as.formula(f)
    y_list<-as.character(f[[2]])
    x_list<-as.character(unlist(strsplit(as.character(f[3]),split=" \\+ ")))
    
    #检查列名对不对
    if(!all(c(y_list,x_list) %in% colnames(dataset))){
      output_result<-list()
      print("方程中有变量不在数据列名中！请检查后再运行！")
      return(output_result)
    }
    
    #筛选出式子中包含的列，并筛选数据完整的???
    clean_dataset<-dataset[,(c(y_list,x_list))]
    clean_dataset<-clean_dataset[(complete.cases(clean_dataset)),]
    #生成空的输出对象
    output_result<-list()
    output_result$f<-f
    output_result$y_list<-y_list
    output_result$x_list<-x_list
    output_result$N_fold<-N_fold
    output_result$Iterations<-Iterations
    output_result$clean_dataset<-clean_dataset
    
    output_result$AUCs<-data.frame(matrix(NA,output_result$Iterations,3),stringsAsFactors = FALSE)
    colnames(output_result$AUCs)<-c("Iterations","AUC","C_index")
    output_result$AUCs$Iterations<-1:output_result$Iterations
    
    #计算需要随机分组的次数
    round<-ceiling(Iterations/N_fold)
    
    #生成初始化的迭代计数???
    inter_counter=0
    
    #进度???
    pb <- tkProgressBar("进度","已完??? %", 0, 100)
    
    #开始计???
    if(output_result$N_fold>nrow(output_result$clean_dataset)){
      print("交叉验证折数大于数据中变量完整的行数???")
    }else{
      for(index_round in 1:round){
        #随机分组
        rand_seq<-rep(1:output_result$N_fold,ceiling(nrow(output_result$clean_dataset)/output_result$N_fold))[1:nrow(output_result$clean_dataset)]
        rand_seq<-sample(rand_seq, nrow(output_result$clean_dataset), replace = FALSE)
        #计算该round内各折的结果
        for(index_fold in 1:output_result$N_fold){
          #检验迭代次数是否足???
          if(inter_counter==output_result$Iterations){
            break
          }else{
            
            #迭代次数不够继续计算
            #划分训练集和验证???
            temp_training_dataset<-output_result$clean_dataset[(rand_seq!=index_fold),]
            temp_validation_dataset<-output_result$clean_dataset[(rand_seq==index_fold),]
            ddist <<- datadist(temp_training_dataset)
            options(datadist='ddist')
            #建模
            if(length(table(temp_validation_dataset[,output_result$y_list]))==1){
              #验证集的结局事情全为0或全???1时无法计算AUC，默认为NA
            }else{
              f_reg <- tryCatch(lrm(output_result$f, data=temp_training_dataset, x=TRUE, y=TRUE,maxit=1000),error=function(e){return(list(fail=TRUE))} )
              if(f_reg$fail | 'try-error' %in% class(f_reg)){           # 判断当前循环的try语句中的表达式是否运行正???
                # 若该回归分析无法进行则该轮循环不赋值，默认为NA
              }else{
                #验证集中计算AUC
                pred_f_validation<-predict(f_reg,temp_validation_dataset)
                modelroc <- roc(temp_validation_dataset[,output_result$y_list],pred_f_validation,quiet=TRUE)
                output_result$AUCs[(inter_counter+1),"AUC"]<-modelroc$auc[1]
                #验证集中计算C-index
                temp_pred_list<-data.frame(grop=temp_validation_dataset[,output_result$y_list],pred=pred_f_validation)
                #所有样本两两配???
                valid_num=0
                concordance_num=0
                for(index_temp_pred_list_1 in 1:nrow(temp_pred_list)){
                  for(index_temp_pred_list_0 in 1:nrow(temp_pred_list)){
                    if((temp_pred_list$grop[index_temp_pred_list_1]==1)&(temp_pred_list$grop[index_temp_pred_list_0]==0)){
                      #当匹配的对子中第一个为group=1（实验组???;第二个为group=0（对照组）时为有效对???
                      valid_num<-valid_num+1
                      if(temp_pred_list$pred[index_temp_pred_list_1]>temp_pred_list$pred[index_temp_pred_list_0]){
                        concordance_num<-concordance_num+1
                      }
                    }
                  }
                }
                #计算C-index
                output_result$AUCs[(inter_counter+1),"C_index"]<-concordance_num/valid_num
              }
            }
            #进度???
            info <- sprintf("已完??? %d%%",round(inter_counter*100/output_result$Iterations))
            setTkProgressBar(pb, inter_counter*100/output_result$Iterations, sprintf("进度 (%s)", info), info)
            inter_counter<-inter_counter+1
          }
        }
      }
    }
    output_result$mean_AUC<-mean(output_result$AUCs$AUC,na.rm = TRUE)
    output_result$mean_C_index<-mean(output_result$AUCs$C_index,na.rm = TRUE)
    output_result$print_report<-paste0("本次验证方法: ",output_result$N_fold," 折交叉验证\n"
                                       ,"迭代次数: ",output_result$Iterations," 次\n"
                                       ,"其中迭代失败次数: ",sum(is.na(output_result$AUCs$AUC))," 次\n"
                                       ,"计算所得的平均AUC: ",output_result$mean_AUC,"\n"
                                       ,"计算所得的平均C_index: ",output_result$mean_C_index,"\n")
    return(output_result)
    
  }
  nomogram_Jackknife_validation<-function(f,dataset){
    
    
    #提取公式中的因变量和自变???
    f<-as.formula(f)
    y_list<-as.character(f[[2]])
    x_list<-as.character(unlist(strsplit(as.character(f[3]),split=" \\+ ")))
    #检查列名对不对
    if(!all(c(y_list,x_list) %in% colnames(dataset))){
      output_result<-list()
      print("方程中有变量不在数据列名中！请检查后再运行！")
      return(output_result)
    }
    #筛选出式子中包含的列，并筛选数据完整的???
    clean_dataset<-dataset[,(c(y_list,x_list))]
    clean_dataset<-clean_dataset[(complete.cases(clean_dataset)),]
    
    N_fold=nrow(clean_dataset)
    Iterations=nrow(clean_dataset)
    
    #生成空的输出对象
    output_result<-list()
    output_result$f<-f
    output_result$y_list<-y_list
    output_result$x_list<-x_list
    output_result$N_fold<-N_fold
    output_result$Iterations<-Iterations
    output_result$clean_dataset<-clean_dataset
    
    output_result$predicted_value<-data.frame(matrix(NA,output_result$Iterations,3),stringsAsFactors = FALSE)
    colnames(output_result$predicted_value)<-c("Iterations","actual_group","predicted_value")
    output_result$predicted_value$Iterations<-1:output_result$Iterations
    
    #计算需要随机分组的次数
    round<-ceiling(Iterations/N_fold)
    
    #进度???
    pb <- tkProgressBar("进度","已完??? %", 0, 100)
    
    #生成初始化的迭代计数???
    inter_counter=0
    
    #开始计???
    if(output_result$N_fold>nrow(output_result$clean_dataset)){
      print("交叉验证折数大于数据中变量完整的行数???")
    }else{
      for(index_round in 1:round){
        #无需随机分组
        rand_seq<-1:nrow(output_result$clean_dataset)
        #计算该round内各折的结果
        for(index_fold in 1:output_result$N_fold){
          #检验迭代次数是否足???
          if(inter_counter==output_result$Iterations){
            break
          }else{
            #迭代次数不够继续计算
            #划分训练集和验证???
            temp_training_dataset<-output_result$clean_dataset[(rand_seq!=index_fold),]
            temp_validation_dataset<-output_result$clean_dataset[(rand_seq==index_fold),]
            ddist <<- datadist(temp_training_dataset)
            options(datadist='ddist')
            #建模
            f_reg <- tryCatch(lrm(output_result$f, data=temp_training_dataset, x=TRUE, y=TRUE,maxit=1000),error=function(e){return(list(fail=TRUE))} )
            if(f_reg$fail | 'try-error' %in% class(f_reg)){           # 判断当前循环的try语句中的表达式是否运行正???
              # 若该回归分析无法进行则该轮循环不赋值，默认为NA
            }else{
              #验证集中计算预测???
              pred_f_validation<-predict(f_reg,temp_validation_dataset)
              output_result$predicted_value[(inter_counter+1),"predicted_value"]<-pred_f_validation
              output_result$predicted_value[(inter_counter+1),"actual_group"]<-temp_validation_dataset[,output_result$y_list]
            }
            #进度???
            info <- sprintf("已完??? %d%%",round(inter_counter*100/output_result$Iterations))
            setTkProgressBar(pb, inter_counter*100/output_result$Iterations, sprintf("进度 (%s)", info), info)
            inter_counter<-inter_counter+1
          }
        }
      }
    }
    #计算AUC
    modelroc <- roc(output_result$predicted_value$actual_group,output_result$predicted_value$predicted_value,quiet=TRUE)
    output_result$AUC_value<-modelroc$auc[1]
    #计算C-index
    temp_pred_list<-data.frame(grop=output_result$predicted_value$actual_group,pred=output_result$predicted_value$predicted_value)
    #所有样本两两配???
    valid_num=0
    concordance_num=0
    for(index_temp_pred_list_1 in 1:nrow(temp_pred_list)){
      for(index_temp_pred_list_0 in 1:nrow(temp_pred_list)){
        if((temp_pred_list$grop[index_temp_pred_list_1]==1)&(temp_pred_list$grop[index_temp_pred_list_0]==0)){
          #当匹配的对子中第一个为group=1（实验组???;第二个为group=0（对照组）时为有效对???
          valid_num<-valid_num+1
          if(temp_pred_list$pred[index_temp_pred_list_1]>temp_pred_list$pred[index_temp_pred_list_0]){
            concordance_num<-concordance_num+1
          }
        }
      }
    }
    output_result$C_index_value<-concordance_num/valid_num
    output_result$print_report<-paste0("本次验证方法: Jackknife验证(留一法交叉验证，leave-one-out cross validation)\n"
                                       ,"迭代次数: ",output_result$Iterations," 次\n"
                                       ,"计算所得的AUC: ",output_result$AUC_value,"\n"
                                       ,"计算所得的C_index: ",output_result$C_index_value,"\n")
    return(output_result)
  }
  nomogram_Bootstrap_validation<-function(f,dataset,Iterations=200,seed=NULL){
    
    #设置随机数种???
    if(!is.null(seed)){
      set.seed(seed)
    }
    
    #提取公式中的因变量和自变???
    f<-as.formula(f)
    y_list<-as.character(f[[2]])
    x_list<-as.character(unlist(strsplit(as.character(f[3]),split=" \\+ ")))
    
    #检查列名对不对
    if(!all(c(y_list,x_list) %in% colnames(dataset))){
      output_result<-list()
      print("方程中有变量不在数据列名中！请检查后再运行！")
      return(output_result)
    }
    
    #筛选出式子中包含的列，并筛选数据完整的???
    clean_dataset<-dataset[,(c(y_list,x_list))]
    clean_dataset<-clean_dataset[(complete.cases(clean_dataset)),]
    #生成空的输出对象
    output_result<-list()
    output_result$f<-f
    output_result$y_list<-y_list
    output_result$x_list<-x_list
    
    output_result$Iterations<-Iterations
    output_result$clean_dataset<-clean_dataset
    
    output_result$AUCs<-data.frame(matrix(NA,output_result$Iterations,3),stringsAsFactors = FALSE)
    colnames(output_result$AUCs)<-c("Iterations","AUC","C_index")
    output_result$AUCs$Iterations<-1:output_result$Iterations
    
    #进度???
    pb <- tkProgressBar("进度","已完??? %", 0, 100)
    
    #生成初始化的迭代计数???
    inter_counter=0
    
    #开始计???
    
    for(index_iterations in 1:output_result$Iterations){
      #随机分组
      rand_seq<-1:nrow(output_result$clean_dataset)
      rand_seq<-sample(rand_seq, nrow(output_result$clean_dataset), replace = TRUE)
      
      if(inter_counter==output_result$Iterations){
        break
      }else{
        #迭代次数不够继续计算
        #划分训练集和验证???
        temp_training_dataset<-output_result$clean_dataset
        temp_validation_dataset<-output_result$clean_dataset[c(rand_seq),]
        ddist <<- datadist(temp_training_dataset)
        options(datadist='ddist')
        #建模
        if(length(table(temp_validation_dataset[,output_result$y_list]))==1){
          #验证集的结局事情全为0或全???1时无法计算AUC，默认为NA
        }else{
          f_reg <- tryCatch(lrm(output_result$f, data=temp_training_dataset, x=TRUE, y=TRUE,maxit=1000),error=function(e){return(list(fail=TRUE))} )
          if(f_reg$fail | 'try-error' %in% class(f_reg)){           # 判断当前循环的try语句中的表达式是否运行正???
            # 若该回归分析无法进行则该轮循环不赋值，默认为NA
          }else{
            #验证集中计算AUC
            pred_f_validation<-predict(f_reg,temp_validation_dataset)
            modelroc <- roc(temp_validation_dataset[,output_result$y_list],pred_f_validation,quiet=TRUE)
            output_result$AUCs[(inter_counter+1),"AUC"]<-modelroc$auc[1]
            #验证集中计算C-index
            temp_pred_list<-data.frame(grop=temp_validation_dataset[,output_result$y_list],pred=pred_f_validation)
            #所有样本两两配???
            valid_num=0
            concordance_num=0
            for(index_temp_pred_list_1 in 1:nrow(temp_pred_list)){
              for(index_temp_pred_list_0 in 1:nrow(temp_pred_list)){
                if((temp_pred_list$grop[index_temp_pred_list_1]==1)&(temp_pred_list$grop[index_temp_pred_list_0]==0)){
                  #当匹配的对子中第一个为group=1（实验组???;第二个为group=0（对照组）时为有效对???
                  valid_num<-valid_num+1
                  if(temp_pred_list$pred[index_temp_pred_list_1]>temp_pred_list$pred[index_temp_pred_list_0]){
                    concordance_num<-concordance_num+1
                  }
                }
              }
            }
            #计算C-index
            output_result$AUCs[(inter_counter+1),"C_index"]<-concordance_num/valid_num
          }
        }
        #进度???
        info <- sprintf("已完??? %d%%",round(inter_counter*100/output_result$Iterations))
        setTkProgressBar(pb, inter_counter*100/output_result$Iterations, sprintf("进度 (%s)", info), info)
        inter_counter<-inter_counter+1
      }
      
    }
    
    output_result$mean_AUC<-mean(output_result$AUCs$AUC,na.rm = TRUE)
    output_result$mean_C_index<-mean(output_result$AUCs$C_index,na.rm = TRUE)
    output_result$print_report<-paste0("本次验证方法: Bootstrap验证\n"
                                       ,"迭代次数: ",output_result$Iterations," 次\n"
                                       ,"其中迭代失败次数: ",sum(is.na(output_result$AUCs$AUC))," 次\n"
                                       ,"计算所得的平均AUC: ",output_result$mean_AUC,"\n"
                                       ,"计算所得的平均C_index: ",output_result$mean_C_index,"\n")
    return(output_result)
    
  }
}



####内部验证####
#logistic
a<- read.csv("D:/Desktop/R.csv",header = TRUE)
f<-"Label ~ NLR + D + Radscore"
#f <- cph(Surv(time, status) ~ agec + pr + pathscat + ln_yesno, x=T, y=T, surv=T, data=bc, time.inc=36)
b<-nomogram_cross_validation(f=f,dataset=a,N_fold=3,Iteration=200,seed=1234)


#cross validation, N_fold设置折数，Iterations设置迭代次数（一般大???200???
a<-nomogram_cross_validation(f=f,dataset=a,N_fold=3,Iterations=200,seed=1234)
cat(a$print_report)
a<-nomogram_cross_validation(f=f,dataset=a,N_fold=5,Iterations=200,seed=1234)
cat(a$print_report)
a<-nomogram_cross_validation(f=f,dataset=a,N_fold=10,Iterations=200,seed=1234)
cat(a$print_report)
#Jackknife validation
b<-nomogram_Jackknife_validation(f=f,dataset=a)
cat(b$print_report)
#Bootstrap validation，Iterations设置迭代次数（一般大???200???
c<-nomogram_Bootstrap_validation(f=f,dataset=R,Iterations=200,seed=NULL)
cat(c$print_report)



####多组ROC####
set1<-read.xlsx('C:/Users/Administrator/Desktop/Venn.xlsx',sheet= "Sheet1",sep=',')
set2<-read.xlsx('C:/Users/Administrator/Desktop/Venn.xlsx',sheet= "Sheet2",sep=',')
set3<-read.xlsx('C:/Users/Administrator/Desktop/Venn.xlsx',sheet= "Sheet3",sep=',')

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











