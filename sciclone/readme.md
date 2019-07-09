```
FF<-c("A","B","C","E","F","G","H","I","J","K","M","O")
for(ff in FF){
  f1<-paste("001",ff,"_Tumor1.vafs",sep="")
  f2<-paste("001",ff,"_Tumor2.vafs",sep="")
  f1
  f2
  v1 = read.table(f1,header=F)
  v2 = read.table(f2,header=F)
  v1 = v1[,c(1,2,10)]
  v2 = v2[,c(1,2,10)]

  names = c("Lung","Brain")
  print(ff)
  sc = sciClone(vafs=list(v1,v2),sampleNames=names[1:2],maximumClusters=20)
  sc.plot2d(sc,paste(ff,"vaf.2d.pdf",sep="."))
}
```
