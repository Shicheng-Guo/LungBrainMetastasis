# install.packages("BiocManager")
# BiocManager::install("maftools")
# BiocManager::install("digest")
library("maftools")

setwd("~/hpc/project/LungBrainMetastasis/vcf/annovar")
setwd("C:\\Users\\shg047\\Documents\\GitHub\\LungBrainMetastasis\\annovar")
annovar<-list.files(pattern="*.T1.hg19_multianno.txt")
laml<-annovarToMaf(annovar=annovar, Center = NULL, refBuild = "hg19",
                   tsbCol = NULL, table = "refGene", basename = NULL, sep = "\t",
                   MAFobj = T, sampleAnno = NULL)
getSampleSummary(laml)
getGeneSummary(laml)
pdf("Lung.pdf")
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = T)
dev.off()
pdf("Lung-S1.pdf")
oncoplot(maf = laml, top = 50,draw_titv = TRUE)
dev.off()

annovar<-list.files(pattern="*.T2.hg19_multianno.txt")
laml<-annovarToMaf(annovar=annovar, Center = NULL, refBuild = "hg19",
                   tsbCol = NULL, table = "refGene", basename = NULL, sep = "\t",
                   MAFobj = T, sampleAnno = NULL)
getSampleSummary(laml)
getGeneSummary(laml)
pdf("Lung_Brain.pdf")
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
pdf("Lung_Brain-S1.pdf")
oncoplot(maf = laml, top = 50, draw_titv = TRUE,fontSize = 0.4)
dev.off()

? oncoplot


annovar<-list.files(pattern="*.T*.hg19_multianno.txt")
laml<-annovarToMaf(annovar=annovar, Center = NULL, refBuild = "hg19",
                   tsbCol = NULL, table = "refGene", basename = NULL, sep = "\t",
                   MAFobj = T, sampleAnno = NULL)
getSampleSummary(laml)
getGeneSummary(laml)
pdf("Brain.pdf")
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
pdf("Brain-S1.pdf")
oncoplot(maf = laml, top = 50, draw_titv = TRUE,fontSize = 0.4)
dev.off()


