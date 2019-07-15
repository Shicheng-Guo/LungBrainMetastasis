# devtools::install_github("kgori/sigfit", args = "--preclean", build_vignettes = TRUE)
# install.packages("factoextra")
library("sigfit")
library("devtools")
library("factoextra")

data<-read.table("C:\Users\shg047\Documents\GitHub\LungBrainMetastasis\vaf\map2matrix.txt",head=T)
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
