# install.packages("openxlsx")
# install.packages("devtools")
# install.packages("GSEABase")

library("openxlsx")
library("devtools")
library("GSEABase")
devtools::install_github("kgori/sigfit", args = "--preclean", build_vignettes = TRUE)

ci95<-function(x){
  error <- qt(0.975,df=length(x)-1)*sd(x)/sqrt(length(x))
  m<-round(mean(x),2)
  d<-round(mean(x)-error,2)
  u<-round(mean(x)+error,2)
  paste("mean=",m, ", 95%CI:",d,"-",u,sep="")
}

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/LungBrainMetastasis")
data<-read.xlsx("Result.xlsx",sheet=6,rowNames=T)
ci95(data[,9])
ci95(data[seq(1,26,by=2),9])
ci95(data[seq(2,26,by=2),9])
head(data)

data<-read.xlsx("Result.xlsx",sheet=7,rowNames=T)
ci95(data[1:14,8])

data<-read.xlsx("Result.xlsx",sheet=4,rowNames=F)
db<-read.xlsx("Result.xlsx",sheet=9,rowNames=F)

temp<-db[match(data[,1],as.character(db[,1])),]
out=na.omit(data.frame(data,temp))
out[1:5,1:5]
mean(apply(out[,2:27],2,function(x) sum(x>0)))
write.table(out,file="LungBrain-Pancancer.txt",sep="\t",quote=F)


data<-read.xlsx("Result.xlsx",sheet=4,rowNames=F)
head(data,n=30)
par(las=1,cex.axis=0.6)
barplot(data$Sum2[1:30]/26,col="red",names.arg=data[1:30,1],horiz = T,xlim=c(0,0.5),xlab="Freqency of Mutation")

# grep -v '#' *.vcf | perl -lane '{print "$1>$2" if $_=~/(\w)>(\w)/}' > sub.txt
data<-read.table("sub.txt")
head(data)

pdf("Figure2.pdf")
par(las=1,cex.axis=1)
barplot(table(data[,1]),col="red",horiz = T,xlab="Counts of different Mutation type")
dev.off()
write.table(table(data[,1]),file="Transition-transversion.txt",sep="\t",quote=F)


data<-read.xlsx("Result.xlsx",sheet=10,rowNames=F)
head(data,n=30)
pdf("pathway.enrichment.pdf")
par(las=2,cex.axis=0.5)
barplot(data$Fold.Enrichment,col="red",names.arg=data[,2],ylab="Fold of Enrichment")
dev.off()

data<-read.xlsx("Result.xlsx",sheet=11,rowNames=F)
head(data,n=30)
pdf("Kewwords.enrichment.pdf")
par(las=2,cex.axis=0.75)
barplot(data$Fold.Enrichment,col="red",names.arg=data[,2],ylab="Fold of Enrichment",ylim=c(0,3))
dev.off()


library("survival")
library("survminer")
data<-read.xlsx("Result.xlsx",sheet=1,rowNames=T)
mutation<-read.table("./bed/LB.MutationProfile.txt",head=T,sep="\t",row.names = 1,check.names = F)

Lmut<-mutation[,seq(1,26,by=2)]
Bmut<-mutation[,seq(2,26,by=2)]

head(Lmut)

P<-c()
j<-c()
for(i in 1:nrow(Lmut)){
  if(sum(Lmut[i,])>0){
    print(i)
    fit<-coxph(Surv(OS,status)~unlist(Lmut[i,]),data)
    fitv<-summary(fit)
    P<-rbind(P,fitv$coefficients)
    j<-c(j,i)
}
}
rownames(P)<-rownames(Lmut)[j]
write.table(P,file="Lung.Mutation.Survival.txt",sep="\t",quote=F,col.names = NA,row.names = T)

newp<-subset(P,P[,5]<0.05)
head(newp)
dim(newp)

for(i in 1:nrow(newp)){
print(i)
gene=rownames(newp)[i]
mut<-unlist(Lmut[match(gene,rownames(Lmut)),])
mut[mut>1]<-1
data$mut<-mut
fit <- survfit(Surv(OS,status)~mut, data = data)
survp<-ggsurvplot(fit, data = data,conf.int = F,pval = TRUE,
           fun = "pct",risk.table = TRUE,size = 1,linetype = "strata",
           palette = c("#E7B800","#2E9FDF"),
           legend = "bottom",legend.title = gene,
           legend.labs = c("wild","Mutation"))
ggsave(file = paste(gene,"Lung.pdf",sep="."), survp$plot)
}



P<-c()
j<-c()
for(i in 1:nrow(Bmut)){
  if(sum(Bmut[i,])>0){
    print(i)
    fit<-coxph(Surv(OS,status)~unlist(Bmut[i,]),data)
    fitv<-summary(fit)
    P<-rbind(P,fitv$coefficients)
    j<-c(j,i)
  }
}
rownames(P)<-rownames(Bmut)[j]
write.table(P,file="Brain.Mutation.Survival.txt",sep="\t",quote=F,col.names = NA,row.names = T)
head(P)

newp<-subset(P,P[,5]<0.05)
head(newp)
dim(newp)

for(i in 1:nrow(newp)){
  gene=rownames(newp)[i]
  print(gene)
  mut<-unlist(Bmut[match(gene,rownames(Bmut)),])
  mut[mut>1]<-1
  data$mut<-mut
  fit <- survfit(Surv(OS,status)~mut, data = data)
  survp<-ggsurvplot(fit, data = data,conf.int = F,pval = TRUE,
                    fun = "pct",risk.table = TRUE,size = 1,linetype = "strata",
                    palette = c("#E7B800","#2E9FDF"),
                    legend = "bottom",legend.title = gene,
                    legend.labs = c("wild","Mutation"))
  ggsave(file = paste(gene,"Brain.pdf",sep="."), survp$plot)
}


library("ComplexHeatmap")

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/LungBrainMetastasis/comutation")
mat = read.table(paste0(system.file("extdata", package = "ComplexHeatmap"), 
                        "/tcga_lung_adenocarcinoma_provisional_ras_raf_mek_jnk_signalling.txt"), 
                 header = TRUE,stringsAsFactors=FALSE, sep = "\t")
mat[1:3, 1:3]
mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat=  mat[, -ncol(mat)]
mat = t(as.matrix(mat))
mat[1:3, 1:3]

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
  }
)

col = c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue")

oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "OncoPrint for TCGA Lung Adenocarcinoma, genes in Ras Raf MEK JNK signalling",
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "HOMDEL", "MUT"), 
                                      labels = c("Amplification", "Deep deletion", "Mutation")))

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/LungBrainMetastasis/bed")
data<-read.table("LBM.MutationProfile.26.test.txt",head=T,sep="\t",row.names=1)
target<-c("TP53","KRAS","FAT4","STK11","EGFR","KMT2C","CHEK2P2","MIR4436A","BAGE2","BAGE3","BAGE4","BAGE5","LOC102723769","MIR6077","FER1L4","FRG1BP","FRG1DP","LOC100996724","AHNAK2","FRG1BP","HMCN2","TPTE")
input<-data[,na.omit(match(target,colnames(data)))]
input

oncoPrint(input, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "OncoPrint for Lung-Brain Metastasis",
          heatmap_legend_param = list(title = "Alternations", at = c("missense_variant", "intron_variant", "frameshift_variant"), 
                                      labels = c("missense_variant", "intron_variant", "frameshift_variant")))


data<-read.table("LBM.MutationProfile.test.txt",head=T,sep="\t",row.names = 1)
data[1:3, 1:3]
target<-unique(c("TP53","KRAS","FAT4","STK11","EGFR","KMT2C","CHEK2P2","ERF","MIR4436A","BAGE2","BAGE3","BAGE4","BAGE5",
          "LOC102723769","MIR6077","FER1L4","FRG1BP","FRG1DP","LOC100996724","AHNAK2","FRG1BP","HMCN2","TPTE","BRSK1",
          "CR1L","CROCCP2","DUS3L","FRG1","MUC17","NOTCH2NL","SDHAP2","SSPO","ALPP","B4GALNT4","BAIAP3","BRSK2","MUC5B"))
input<-t(data[,c(na.omit(match(target,colnames(data))))])
input[1:3, 1:3]

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  missense_variant = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
  },
  regulatory_region_variant = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  splice_related_variants = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#008000", col = NA))
  },
  frameshift_variant = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "green", col = NA))
  },
  Others = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "purple", col = NA))
  }
  
)

col = c("missense_variant" = "red", "regulatory_region_variant" = "blue", "splice_related_variants" = "#008000","frameshift_variant"="green","Others"="purple")

oncoPrint(input, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col=col,
          remove_empty_columns = TRUE,
          heatmap_legend_param = list(title = "Alternations", at = c("missense_variant", "regulatory_region_variant", "splice_related_variants","frameshift_variant","Others"), 
                                      labels = c("missense_variant", "regulatory_region_variant", "splice_related_variants","frameshift_variant","Others")))



devtools::install_github("kgori/sigfit", args = "--preclean", build_vignettes = TRUE)




install.packages("survminer")
library("survival")
library("survminer")
library("openxlsx")
library(survminer)
require("survival")
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/LungBrainMetastasis")
data<-read.xlsx("Result.xlsx",sheet=1,rowNames=T)
P<-c()
for(i in 1:13){
xdata<-data[,c(2,3,4,6,7,8,9,12,13,15)]
xdata<-xdata[-i,]
model <- coxph( Surv(OS,status) ~ .,data =xdata)
ggforest(model,data=xdata)
fit<-summary(model)
P<-c(P,fit$coefficients[8,5])
}
i=11
xdata<-data[,c(2,3,4,6,7,8,9,12,13,15)]
xdata<-xdata[-11,]
model <- coxph( Surv(OS,status) ~ .,data =xdata)
fit<-summary(model)
ggforest(model,data=xdata)
library("meta")
library("ggplot2")
input<-data.frame(fit$conf.int)
write.table(input[,c(1,3,4)],file="ShareRatio_CoxRegression_95CI.txt",col.names = NA,row.names = T,sep="\t",quote=F)
input[3,4]<-10^25
input$Group<-rownames(input)
input$Condition<-1:nrow(input)
input$RiskRatio<-input[,1]
input$LowerLimit<-input[,3]
input$UpperLimit<-input[,4]
p = ggplot(data=input,aes(x = Group,y = RiskRatio, ymin = LowerLimit, ymax = UpperLimit ))+
  geom_pointrange(aes(col=Group))+
  geom_hline(aes(fill=Group),yintercept =1, linetype=2)+
  xlab('Group')+ ylab("Risk Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=LowerLimit, ymax=UpperLimit,col=Group),width=0.5,cex=1)+ 
  facet_wrap(~Condition,strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip()
p + scale_y_log10()+ guides(col = guide_legend(reverse = TRUE))



devtools::install_github("kgori/sigfit", args = "--preclean", build_vignettes = TRUE)
install.packages("factoextra")
library("sigfit")
library("devtools")
library("factoextra")
data<-read.table("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/LungBrainMetastasis/vaf/map2matrix.txt",head=T)
head(data)
name=gsub("X001","",colnames(data))
for(i in 1:12){
temp<-data[,c(2*i-1,2*i)]
temp<-temp[-which(rowSums(temp)==0),]
km.res <- kmeans(temp,8, nstart = 25)
write.table(km.res$cluster,file=paste(name[2*i-1],".cluster8.txt",sep=""),sep="\t",col.names = F,row.names = T,quote=F)
p<-fviz_cluster(km.res, data = temp,geom = c("point"),pointsize = 2)+theme_bw()
ggsave(p,file=paste(name[2*i-1],".8.pdf",sep=""))
}

###
data<-read.table("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/LungBrainMetastasis/vaf/map2matrix.txt",head=T)
head(data)
catNum<-c()
for(i in 1:(ncol(data)/2)){
  temp<-data[,c(2*i-1,2*i)]
  temp<-temp[-which(rowSums(temp)==0),]
  L<-subset(temp,temp[,1]>0 & temp[,2]==0)
  LB<-subset(temp,temp[,1]>0 & temp[,2]>0)
  B<-subset(temp,temp[,1]==0 & temp[,2]>0)
  catNum<-cbind(catNum,c(nrow(L),nrow(LB),nrow(B)))
}
rownames(catNum)<-c("Lung","LB","Brain")
dim(catNum)
length(name)
name=gsub("X001","",colnames(data))
name=gsub("_Tumor1","",name)
colnames(catNum)<-name[seq(1,length(name),by=2)]
catNum<-catNum[,-c(4,13)]
write.table(catNum,file="Mutation.Number.txt",col.names = NA,row.names = T,quote=F)
library("reshape2")
input<-melt(catNum,id=c("A","B","C","E","F","G","H","I","J","K","M","O"))
head(input)
library("ggplot2")
s<-ggplot(input,aes(Var2,fill=value))
s<-s + geom_bar(aes(weight = value,fill=Var1,position="stack"))
s<- s+theme_bw()
ggsave(file="Mutation.Number.png")

####
data<-read.table("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/LungBrainMetastasis/vaf/map2matrix.txt",head=T)
head(data)
symbol<-unlist(lapply(strsplit(rownames(data),":"),function(x) x[3]))
RealDrive<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/LungBrainMetastasis/127NatureDrivenMutation.txt")
data<-data[which(symbol %in% RealDrive[,2]),]

catNum<-c()
for(i in 1:(ncol(data)/2)){
  temp<-data[,c(2*i-1,2*i)]
  temp<-temp[-which(rowSums(temp)==0),]
  L<-subset(temp,temp[,1]>0 & temp[,2]==0)
  LB<-subset(temp,temp[,1]>0 & temp[,2]>0)
  B<-subset(temp,temp[,1]==0 & temp[,2]>0)
  catNum<-cbind(catNum,c(nrow(L),nrow(LB),nrow(B)))
}
rownames(catNum)<-c("Lung","LB","Brain")
dim(catNum)
length(name)
name=gsub("X001","",colnames(data))
name=gsub("_Tumor1","",name)
colnames(catNum)<-name[seq(1,length(name),by=2)]
catNum<-catNum[,-c(4,13)]
write.table(catNum,file="Mutation.Number.1.txt",col.names = NA,row.names = T,quote=F)
library("reshape2")

input<-melt(catNum,id=c("A","B","C","E","F","G","H","I","J","K","M","O"))
head(input)
library("ggplot2")
s<-ggplot(input,aes(Var2,fill=value))
s<-s + geom_bar(aes(weight = value,fill=Var1,position="stack"))
s<- s+theme_bw()
s
ggsave(file="Mutation.Number.pdf")

catNum<-catNum[,-c(match(c("A","F","G"),colnames(catNum)))]
write.table(catNum,file="Mutation.Number.2.txt",col.names = NA,row.names = T,quote=F)
input<-melt(catNum,id=c("B","C","E","H","I","J","K","M","O"))
head(input)
library("ggplot2")
s<-ggplot(input,aes(Var2,fill=value))
s<-s + geom_bar(aes(weight = value,fill=Var1,position="stack"))
s<- s+theme_bw()
s
ggsave(file="Mutation.Number.2.pdf")



url<-"https://raw.githubusercontent.com/Shicheng-Guo/LungBrainMetastasis/master/MutationProfileSignature/Lung_COSMIC_sign_contributions.txt"
MS<-read.table(url,sep="\t",head=T)
head(MS)
MS<-data.frame(MS[,3:ncol(MS)])
head(MS)
colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
pdf("Lung.COSMIC.heatmap.pdf")
heatmap.2(data.matrix(MS),trace="none",density.info="none",keysize=0.9,col=colors)
dev.off()

url<-"https://raw.githubusercontent.com/Shicheng-Guo/LungBrainMetastasis/master/MutationProfileSignature/Brain_COSMIC_sign_contributions.txt"
MS<-read.table(url,sep="\t",head=T)
head(MS)
MS<-data.frame(MS[,3:ncol(MS)])
head(MS)
colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
pdf("Brain.COSMIC.heatmap.pdf")
heatmap.2(data.matrix(MS),trace="none",density.info="none",keysize=0.9,col=colors)
dev.off()


library("factoextra")
data<-read.table("vafmatrix.txt",head=T,check.names = F)
head(data)
name=gsub("001","",colnames(data))
for(i in 1:14){
  temp<-data[,c(2*i-1,2*i)]
  temp<-temp[-which(rowSums(temp)==0),]
  km.res <- kmeans(temp,8, nstart = 25)
  write.table(km.res$cluster,file=paste(name[2*i-1],".cluster8.txt",sep=""),sep="\t",col.names = F,row.names = T,quote=F)
  p<-fviz_cluster(km.res, data = temp,geom = c("point"),pointsize = 2)+theme_bw()
  ggsave(p,file=paste(name[2*i-1],".8.pdf",sep=""))
  print(name[2*i-1])
}


install.packages("NbClust")
library("factoextra")
library("NbClust")
pkgs <- c("factoextra",  "NbClust")
onc<-fviz_nbclust(temp, kmeans,method = "gap_stat")
fviz_nbclust(temp, FUNcluster, method = c("silhouette", "wss", "gap_stat"))
fviz_nbclust(df, kmeans, method = "wss") + geom_vline(xintercept = 4, linetype = 2)+labs(subtitle = "Elbow method")
# Silhouette method
fviz_nbclust(df, kmeans, method = "silhouette")+labs(subtitle = "Silhouette method")
# Gap statistic
set.seed(123)
fviz_nbclust(df, kmeans, nstart = 25,  method = "gap_stat", nboot = 100,verbose = FALSE)+labs(subtitle = "Gap statistic method")