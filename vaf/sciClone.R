
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/LungBrainMetastasis/vaf")
library("factoextra")
data<-read.table("map2matrix.txt",head=T,check.names = F)
head(data)
name=gsub("001","",colnames(data))
for(i in 1:14){
  temp<-data[,c(2*i-1,2*i)]
  temp<-temp[-which(rowSums(temp)==0),]
  km.res <- kmeans(temp,8, nstart = 25)
  write.table(km.res$cluster,file=paste(name[2*i-1],".cluster8.txt",sep=""),sep="\t",col.names = F,row.names = T,quote=F)
  p<-fviz_cluster(km.res, stand=F,data = temp,choose.vars = colnames(temp),geom = c("point"),pointsize = 2)+theme_bw()
  print(p)
  ggsave(p,file=paste(name[2*i-1],".8.pdf",sep=""))
  print(name[2*i-1])
}
