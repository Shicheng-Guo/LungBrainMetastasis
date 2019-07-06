# install.packages("BiocManager")
# BiocManager::install("maftools")
# BiocManager::install("digest")
library("maftools")

setwd("~/hpc/project/LungBrainMetastasis/vcf/annovar")
setwd("C:\\Users\\shg047\\Documents\\GitHub\\LungBrainMetastasis\\annovar")

### Lung
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
pdf("Lung-S2.pdf")
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))
dev.off()

### Brain
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
pdf("Brain-S2.pdf")
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))
dev.off()

### Lung + Brain
annovar<-list.files(pattern="*.T*.hg19_multianno.txt")
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
pdf("Lung_Brain-S2.pdf")
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))
dev.off()

annovar<-list.files(pattern="*.T*.hg19_multianno.txt")
laml<-annovarToMaf(annovar=annovar, Center = NULL, refBuild = "hg19",
                   tsbCol = NULL, table = "refGene", basename = NULL, sep = "\t",
                   MAFobj = T, sampleAnno = NULL)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
oncoplot(maf = laml, top = 50, draw_titv = TRUE,fontSize = 0.4)
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
rainfallPlot(maf = laml, detectChangePoints = TRUE, pointSize = 0.6)
laml.mutload = tcgaCompare(maf = laml, cohortName = 'LBMT')
plotVaf(maf = laml, vafCol = 'i_TumorVAF_WU')
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))
laml.sig = oncodrive(maf = laml, AACol = 'aaChange', minMut = 5, pvalMethod = 'zscore')
laml.sig
write.table(laml.sig,file="oncodrive.txt",sep="\t",col.names = NA,row.names = T,quote=F)
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE)
lollipopPlot(maf = laml, gene = 'FRG1', AACol = 'aaChange', showMutationRate = TRUE)
lollipopPlot(maf = laml, gene = 'AHNAK2', AACol = 'aaChange', showMutationRate = TRUE)
lollipopPlot(maf = laml, gene = 'KMT2C', AACol = 'aaChange', showMutationRate = TRUE)

laml.pfam = pfamDomains(maf = laml, AACol = 'aaChange', top = 10)
laml.pfam$domainSummary[,1:3, with = FALSE]


annovar<-list.files(pattern="*.T1.hg19_multianno.txt")
laml1<-annovarToMaf(annovar=annovar, Center = NULL, refBuild = "hg19",
                   tsbCol = NULL, table = "refGene", basename = NULL, sep = "\t",
                   MAFobj = T, sampleAnno = NULL)
annovar<-list.files(pattern="*.T2.hg19_multianno.txt")
laml2<-annovarToMaf(annovar=annovar, Center = NULL, refBuild = "hg19",
                   tsbCol = NULL, table = "refGene", basename = NULL, sep = "\t",
                   MAFobj = T, sampleAnno = NULL)

pt.vs.rt <- mafCompare(m1 = laml1, m2 = laml2, m1Name = 'Primary', m2Name = 'Relapse', minMut = 2)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.5, color = c('royalblue', 'maroon'), geneFontSize = 0.4)

genes = c("KMT2C", "BAGE2", "ANKRD36C", "AHNAK2", "ADAMTSL4","MST1L","PRAMEF4","PDE4DIP","FLT3LG","DMBT1")
coOncoplot(m1 = laml1, m2 = laml2, m1Name = 'Lung', m2Name = 'Brain', genes = genes, removeNonMutated = TRUE)

