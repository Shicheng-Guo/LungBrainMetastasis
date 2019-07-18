<<<<<<< HEAD
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/LungBrainMetastasis")
data<-read.xlsx("Result.xlsx",sheet=6,rowNames=T)

=======
>>>>>>> 133380ca3c1afef90aafffd78b5621a330ecf703
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
<<<<<<< HEAD
  }
=======
}
>>>>>>> 133380ca3c1afef90aafffd78b5621a330ecf703
}
rownames(P)<-rownames(Lmut)[j]
write.table(P,file="Lung.Mutation.Survival.txt",sep="\t",quote=F,col.names = NA,row.names = T)

newp<-subset(P,P[,5]<0.05)
head(newp)
dim(newp)

for(i in 1:nrow(newp)){
<<<<<<< HEAD
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
=======
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
>>>>>>> 133380ca3c1afef90aafffd78b5621a330ecf703
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
<<<<<<< HEAD
}
=======
}
>>>>>>> 133380ca3c1afef90aafffd78b5621a330ecf703
