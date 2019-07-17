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
write.table(getSampleSummary(laml),file="Lung_Sample_MutationSummary.txt",sep="\t",col.names = NA,row.names = T,quote=F)
write.table(getGeneSummary(laml),file="Lung_Gene_MutationSummary.txt",sep="\t",col.names = NA,row.names = T,quote=F)

annovar<-list.files(pattern="*_T2.hg19_multianno.txt")
laml<-annovarToMaf(annovar=annovar, Center = NULL, refBuild = "hg19",
                   tsbCol = NULL, table = "refGene", basename = NULL, sep = "\t",
                   MAFobj = T, sampleAnno = NULL)
write.table(getSampleSummary(laml),file="Brain_Sample_MutationSummary.txt",sep="\t",col.names = NA,row.names = T,quote=F)
write.table(getGeneSummary(laml),file="Brain_Gene_MutationSummary.txt",sep="\t",col.names = NA,row.names = T,quote=F)

annovar<-list.files(pattern="*_T*.hg19_multianno.txt")
laml<-annovarToMaf(annovar=annovar, Center = NULL, refBuild = "hg19",
                   tsbCol = NULL, table = "refGene", basename = NULL, sep = "\t",
                   MAFobj = T, sampleAnno = NULL)
write.table(getSampleSummary(laml),file="Lung_Brain_Sample_MutationSummary.txt",sep="\t",col.names = NA,row.names = T,quote=F)
write.table(getGeneSummary(laml),file="Lung_Brain_Gene_MutationSummary.txt",sep="\t",col.names = NA,row.names = T,quote=F)
