
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
