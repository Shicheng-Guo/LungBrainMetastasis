# install.packages("BiocManager")
# BiocManager::install("maftools")
# BiocManager::install("digest")
library("maftools")
library("callr")
library("maftools")
install.packages("callr")
install.packages("processx")
install.packages("processx_3.4.0.tar.gz")
devtools::install_github(repo = "PoisonAlien/TCGAmutations")
BiocManager::install("BSgenome")
library("BSgenome", quietly = TRUE)


setwd("~/hpc/project/LungBrainMetastasis/vcf/annovar")
setwd("C:\\Users\\shg047\\Documents\\GitHub\\LungBrainMetastasis\\annovar\\multianno")

### Lung
annovar<-list.files(pattern="*_T1.hg19_multianno.txt")
laml<-annovarToMaf(annovar=annovar, Center = NULL, refBuild = "hg19",
                   tsbCol = NULL, table = "refGene", basename = NULL, sep = "\t",
                   MAFobj = T, sampleAnno = NULL)
Lung<-laml
getSampleSummary(laml)
getGeneSummary(laml)
pdf("Lung.pdf")
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = T)
dev.off()
pdf("Lung-S1.pdf")
oncoplot(maf = laml, top = 50,draw_titv = TRUE,fontSize = 0.8)
dev.off()
pdf("Lung-S2.pdf")
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))
dev.off()


### Brain
annovar<-list.files(pattern="*_T2.hg19_multianno.txt")
laml<-annovarToMaf(annovar=annovar, Center = NULL, refBuild = "hg19",
                   tsbCol = NULL, table = "refGene", basename = NULL, sep = "\t",
                   MAFobj = T, sampleAnno = NULL)
Brain<-laml
getSampleSummary(laml)
getGeneSummary(laml)
pdf("Brain.pdf")
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
pdf("Brain-S1.pdf")
oncoplot(maf = laml, top = 50, draw_titv = TRUE,fontSize = 0.4)
dev.off()
pdf("Brain-S2.pdf")
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))
dev.off()

### Lung + Brain
annovar<-list.files(pattern="*.hg19_multianno.txt")
laml<-annovarToMaf(annovar=annovar, Center = NULL, refBuild = "hg19",
                   tsbCol = NULL, table = "refGene", basename = NULL, sep = "\t",
                   MAFobj = T, sampleAnno = NULL)
getSampleSummary(laml)
getGeneSummary(laml)
pdf("Lung_Brain.pdf")
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
pdf("Lung_Brain-S1.pdf")
oncoplot(maf = laml, top = 50, sampleOrder=SO,removeNonMutated=F,draw_titv = TRUE,fontSize = 0.4,showTumorSampleBarcodes=T)
dev.off()

pdf("Lung_Brain-S2.pdf")
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))
dev.off()


SO<-unlist(lapply(annovar,function(x) unlist(strsplit(x,"[.]"))[1]))
oncoplot(maf = laml, top = 50, sampleOrder=SO,removeNonMutated=F,draw_titv = TRUE,fontSize = 0.4,showTumorSampleBarcodes=T)
oncoplot(maf = laml, top = 50, draw_titv = TRUE,fontSize = 0.4)

? oncoplot
genes = c("KMT2C", "BAGE2", "ANKRD36C", "AHNAK2", "ADAMTSL4","MST1L","PRAMEF4","PDE4DIP","FLT3LG","DMBT1")
coOncoplot(m1 = Lung, m2 = Brain, m1Name = 'Lung', m2Name = 'Brain', genes = genes, removeNonMutated = TRUE)


library('NMF')
BSgenome::available.genomes()
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")


library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)

library('pheatmap')

laml.tnm= trinucleotideMatrix(maf = Lung, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
laml.sign = extractSignatures(mat = laml.tnm, nTry = 6, plotBestFitRes = FALSE)
plotSignatures(laml.sign, title_size = 0.8, )
pheatmap::pheatmap(mat = laml.sign$coSineSimMat, cluster_rows = FALSE, main = "cosine similarity against validated signatures")
laml.se = signatureEnrichment(maf = laml, sig_res = laml.sign)




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

devtools::install_github(repo = "PoisonAlien/TCGAmutations")

library("TCGAmutations")

tcga_load(study = "LUAD") 
tcga_load(study = "LUSC") 
tcga_load(study = "GBM") 
tcga_load(study = "LGG") 

tcga_luad_mc3
tcga_lusc_mc3
tcga_lgg_mc3
tcga_gbm_mc3

pdf("TCGA_luad.pdf")
plotmafSummary(maf = tcga_luad_mc3, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

pdf("TCGA_lusc.pdf")
plotmafSummary(maf = tcga_lusc_mc3, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

pdf("TCGA_lgg.pdf")
plotmafSummary(maf = tcga_lgg_mc3, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

pdf("TCGA_gbm.pdf")
plotmafSummary(maf = tcga_gbm_mc3, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

