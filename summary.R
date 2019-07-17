library("mclust")
library("maftools")
library("callr")
library("maftools")
library("BiocParallel")
library('pheatmap')
library("devtools")

setwd("C:\\Users\\shg047\\Documents\\GitHub\\LungBrainMetastasis\\annovar")

annovar<-list.files(pattern="*_T1.hg19_multianno.txt")
laml<-annovarToMaf(annovar=annovar, Center = NULL, refBuild = "hg19",
                   tsbCol = NULL, table = "refGene", basename = NULL, sep = "\t",
                   MAFobj = T, sampleAnno = NULL)
write.mafSummary(maf = laml, basename = 'Lung.maf')

annovar<-list.files(pattern="*_T2.hg19_multianno.txt")
laml<-annovarToMaf(annovar=annovar, Center = NULL, refBuild = "hg19",
                   tsbCol = NULL, table = "refGene", basename = NULL, sep = "\t",
                   MAFobj = T, sampleAnno = NULL)
write.mafSummary(maf = laml, basename = 'Brain.maf')


annovar<-list.files(pattern="*.hg19_multianno.txt")
laml<-annovarToMaf(annovar=annovar, Center = NULL, refBuild = "hg19",
                   tsbCol = NULL, table = "refGene", basename = NULL, sep = "\t",
                   MAFobj = T, sampleAnno = NULL)
write.mafSummary(maf = laml, basename = 'Lung_Brain.maf')

