BiocManager::install("devtools")
BiocManager::install("IRanges")
BiocManager::install("digest")
BiocManager::install("IRanges")
install.packages("https://cran.rstudio.com/bin/windows/contrib/3.6/processx_3.3.1.zip", repos = NULL)
install.packages("https://cran.rstudio.com/bin/windows/contrib/3.6/callr_3.2.0.zip", repos = NULL)

library('processx')
library("callr")
library("devtools")
library("IRanges")
install_github("genome/bmm")
install_github("genome/sciClone")

library(sciClone)

# wget https://raw.githubusercontent.com/genome/sciclone-meta/master/tests/data/copy_number_tum1 
# wget https://raw.githubusercontent.com/genome/sciclone-meta/master/tests/data/copy_number_tum2
# wget https://raw.githubusercontent.com/genome/sciclone-meta/master/tests/data/copy_number_tum3
# wget https://raw.githubusercontent.com/genome/sciclone-meta/master/tests/data/exclude.chr1
# wget https://raw.githubusercontent.com/genome/sciclone-meta/master/tests/data/regionsToExclude
# wget https://raw.githubusercontent.com/genome/sciclone-meta/master/tests/data/vafs.dat
# wget https://raw.githubusercontent.com/genome/sciclone-meta/master/tests/shortTest.R


library(sciClone)

setwd("C:/Users/shg047/Documents/GitHub/LungBrainMetastasis/annovar")
annovar<-list.files(pattern="*.T1.hg19_multianno.txt")
annovar
laml<-annovarToMaf(annovar=annovar, Center = NULL, refBuild = "hg19",
                   tsbCol = NULL, table = "refGene", basename = NULL, sep = "\t",
                   MAFobj = F, sampleAnno = NULL)

laml

setwd("C:/Users/shg047/Documents/GitHub/LungBrainMetastasis/sciclone/")
v = read.table("vafs.dat",header=T);
v1 = v[1:100,c(1,2,8,9,10)]
v2 = v[1:100,c(1,2,11,12,13)]
v3 = v[1:100,c(1,2,14,15,16)]

head(v1)
head(v2)
head(v3)

regions = read.table("exclude.chr1")
cn1 = read.table("copy_number_tum1")
cn1 = cn1[,c(1,2,3,5)]
cn2 = read.table("copy_number_tum2")
cn2 = cn2[,c(1,2,3,5)]
cn3 = read.table("copy_number_tum3")
cn3 = cn3[,c(1,2,3,5)]

head(cn1)
head(cn2)
head(cn3)

names = c("Sample1","Sample2","Sample3")
reg1 = read.table("regionsToExclude")

#make an output directory, deleting old results first if they exist
suppressWarnings(dir.create("results"))
unlink("results/*", recursive=TRUE)

#run one sample
sc = sciClone(vafs=v1,copyNumberCalls=cn1,sampleNames=names[1],regionsToExclude=reg1)
writeClusterTable(sc, "results/clusters1")
sc.plot1d(sc,"results/clusters1.1d.pdf")

#run two samples
sc = sciClone(vafs=list(v1,v2),copyNumberCalls=list(cn1,cn2),sampleNames=names[1:2])
writeClusterTable(sc, "results/clusters2")
sc.plot1d(sc,"results/clusters2.1d.pdf")
sc.plot2d(sc,"results/clusters2.2d.pdf")

#run three samples
sc = sciClone(vafs=list(v1,v2,v3),copyNumberCalls=list(cn1,cn2,cn3),sampleNames=names,regionsToExclude=list(reg1,reg1))
if(!(is.null(sc))){
  print("ERROR - this should have failed, because there are no cn-neutral points in all three samples")
}


cat("\n")
cat("=========================================================\n")
cat("Test 3.1 - three samples - should succeed")
cat("\n")
#run two samples
sc = sciClone(vafs=list(v1,v2,v3),
              copyNumberCalls=list(cn1,cn2,cn2),
              sampleNames=names,
              regionsToExclude=list(reg1,reg1))
writeClusterTable(sc, "results.new/clusters3")
sc.plot1d(sc,"results.new/clusters3.1d.pdf")
sc.plot2d(sc,"results.new/clusters3.2d.pdf")
sc.plot3d(sc, sc@sampleNames, size=700, outputFile="results.new/clusters3.3d.gif")


#### LBMT

v = read.table("vafs.dat",header=T);
v1 = v[1:100,c(1,2,8,9,10)]
names = c("Sample1","Sample2")
sc = sciClone(vafs=list(v1,v2),sampleNames=names[1:2])
sc.plot2d(sc,"results/clusters2.2d.pdf")


setwd("C:/Users/shg047/Documents/GitHub/LungBrainMetastasis/vaf")
FF<-c("A","B","C","E","F","G","H","I","J","K","M","O")

ff = FF[1]
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


