data<-read.table("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/LungBrainMetastasis/vaf/vafmatrix.txt",head=T)
install.packages("factoextra")
library("factoextra")
name=gsub("X001","",colnames(data))
for(i in 1:12){
temp<-na.omit(data[,c(2*i-1,2*i)])
km.res <- kmeans(temp, 5, nstart = 25)
p<-fviz_cluster(km.res, data = temp,geom = c("point"),pointsize = 2)+theme_bw()
ggsave(p,file=paste(name[2*i-1],".pdf",sep=""))
}
