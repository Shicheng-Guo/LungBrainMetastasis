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






laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
lollipopPlot(maf = laml, gene = 'DNMT3A', AACol = 'Protein_Change', showMutationRate = TRUE)
lollipopPlot(maf = laml, gene = 'KIT', AACol = 'Protein_Change', labelPos = 816, refSeqID = 'NM_000222')
rainfallPlot(maf = laml, detectChangePoints = TRUE, pointSize = 0.6)
laml.mutload = tcgaCompare(maf = laml, cohortName = 'LBMT')
plotVaf(maf = laml, vafCol = 'i_TumorVAF_WU')
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))
laml.sig = oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE)

getFields(laml)
laml$aaChange

