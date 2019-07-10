
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
