################ DEFINICION EXPERIMENTO #################################

source("./code/require.r", chdir=T)
source("./code/auxiliares.r", chdir=T)
source("./code/params.r", chdir=T)
source("./code/graficos.r", chdir=T)
source("./code/ChAMP/quality control plots.r", chdir=T)
source('C:/PROYECTOS/BIO12_AL_007/code/Analisis/Locus-by-locus analysis/multiRegMultistat.R')

################ DEFINICION EXPERIMENTO #################################

source("./code/experiment.r", chdir=T)

################ CARGA DATOS ############################################

sampleDescriptions = read.table(SampleDescFileName,header=TRUE,sep=",")
# meditacion
#sampleDescriptions = read.table(SampleDescFileName,header=TRUE,sep="\t")

colnames(sampleDescriptions) <- c("grupo")
sampleDescriptions$SampleID <- rownames(sampleDescriptions)

colnames(sampleDescriptions) <- c("SampleID","estadio","sexo","edad")
rownames(sampleDescriptions) <- sampleDescriptions$SampleID
sampleDescriptions$SampleLabel <- sampleDescriptions$SampleID

meth450_data_original <-methylumiR(filename=MethylDataFileName, sampleDescriptions=sampleDescriptions )
explore450(meth450_data_original, group="estadio")
sapply(colnames(meth450_data_original), function(x){print predictSex(meth450_data_original,x)})
load(paste("./data",projectName,"meth450_data_original.Rdata",sep="/"))

################ REMOVE BAD SAMPLES ############################################

meth450_data <- meth450_data_original[,!(pData(meth450_data_original)[,1] %in% remove_samples)]
explore450(meth450_data, group="sexo")

################ FILTER PROBES ############################################

################
# controles Illumina
################

meth450_data <- meth450_data[!is.na(as.array(meth450_data@featureData@data$MAPINFO)),] 

################
# ch probes
################

table(substr(rownames(meth450_data_original),1,2))

# cg      ch     rs 
# 482421  3091   65 

meth450_data <- meth450_data[substr(rownames(meth450_data),1,2)!='ch',] 

################
# bad quality probes
################
# 1328 sites were removed as beadcount <3 in 5 % of samples 406 sites having 1 % of samples with a detection p-value greater than 0.05 were removed
meth450_data <- pfilter(meth450_data)

# sum(meth450_data@featureData@data$PROBE_SNPS_10!="")
#[1] 36535
# sum(meth450_data@featureData@data$PROBE_SNPS!="")
#[1] 59234

# sum(!is.na(SNPs.138$Probe_rs))
#[1] 87317

################
# PROBE_SNPS_10
################

meth450_data <- meth450_data[as.array(meth450_data@featureData@data$PROBE_SNPS_10)=="",]

################
# chrX
################

meth450_data <- meth450_data[as.array(!(meth450_data@featureData@data$CHR)=="X"),]

################
# chrY
################

meth450_data <- meth450_data[as.array(!(meth450_data@featureData@data$CHR)=="Y"),]

################
# non_specific_probes
################

non_specific_probes = read.table("./data/lista_29233.txt")
sum(rownames(meth450_data) %in% non_specific_probes[,1])
meth450_data <- meth450_data[!(rownames(meth450_data) %in% non_specific_probes[,1]),]

################
# noisy probes
################

extended_annotation = read.table("./data/1471-2164-15-51-s2.csv",sep=",",header=T)
probes_to_filter <- extended_annotation$probe[(extended_annotation$Flag.discard.keep. == "discard")]

sum(rownames(meth450_data) %in% probes_to_filter)

meth450_data <- meth450_data[!(rownames(meth450_data) %in% probes_to_filter),]

################
# probes left
################

nrow(meth450_data)
# 265373 left

explore450(meth450_data, group="estadio")
save(meth450_data, file=paste("./data",projectName,"meth450_data.dat",sep="/"))
load(paste("./data",projectName,"meth450_data.dat",sep="/"))

################
# Normalize
################

meth450_data <- BMIQ(meth450_data)
meth450_data <- SWAN(meth450_data)
meth450_data <- dasen(meth450_data)

explore450(meth450_data, group="estadio")

save(meth450_data, file=paste("./data",projectName,"meth450_norm.dat",sep="/"))
load(paste("./data",projectName,"meth450_norm.dat",sep="/"))

# meth450_norm <- meth450_norm[,!(pData(meth450_norm)[,2] %in% remove_samples)]

################################
# Differential methylation
################################

samps <- pData(meth450_data)
b <- betas(meth450_data)
M <- exprs(meth450_data)
M <- getM(meth450_data) # for signature ‘"MethylSet"’

## samps$estadio <- sampleDescriptions[rownames(samps),"estadio"] ## correccion

# alzheimer
clases = c("ALZHEIMER","CONTROL")[as.factor(samps$estadio=="EA_CONTROL")]
clases <-as.factor(clases)

##
## correlation with other variables
##

cor_M_edad <- apply(M,1,cor, y =samps$edad)
hist(cor_M_edad)

amieloide <- read.csv('data/alzheimer/AMILOIDE_081014.csv')
head(M[,amieloide$SampleName])
cor_M_amieloide <- apply(M[,amieloide$SampleName],1,cor, y =amieloide$Mean.total.area)
hist(cor_M_amieloide)

###################################################
### Wilcox 
###################################################

data.test.M <- sitetest_local(beta=M, grouplev=clases, gcase="ALZHEIMER", gcontrol="CONTROL",concov=concov,testmethod = testmethod, Padj=Padj, paired = paired)

plotMDS(M[rownames(topW),], labels=samps$sampleID, col=as.integer(clases))

###################################################
### code chunk number 7: mdsplot
###################################################
par(mfrow=c(1,1))
plotMDS(M, labels=samps$sampleID, col=as.integer(clases))
legend("topleft",legend=c("ALZHEIMER","CONTROL"),pch=16,cex=1.2,col=1:2)


###################################################
### code chunk number 8: design
###################################################
library(limma)
id <- factor(samps$sampleID)
design <- model.matrix(~clases)
design


###################################################
### code chunk number 9: diffmeth
###################################################
fit <- lmFit(M,design)
fit <- eBayes(fit)


###################################################
### code chunk number 10: results
###################################################
summary(decideTests(fit))
top<-topTable(fit,coef=2)
head(top)


plotMDS(b[rownames(top),], labels=samps$sampleID, col=as.integer(clases))
b_heatmap.3(b, rownames(top), clases)


###
### volcano
###
xgene <- function(x) {strsplit(as.character(x),";")[[1]][1]}


gene_list <- topTable(fit,coef=2,number=10000)
gene_list$threshold = as.factor(abs(gene_list$logFC) > 1 & gene_list$adj.P.Val < 0.01)
gene_list$geneName <- sapply(meth450_data_original[rownames(gene_list)]@featureData@data$UCSC_REFGENE_NAME,xgene)

volcano = ggplot(data=gene_list, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  xlim(c(-2.5, 2.5)) + ylim(c(0, 3)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  geom_text(data=subset(gene_list, as.logical(gene_list$threshold)), aes(x=logFC, y=-log10(adj.P.Val),
                  label=geneName, size=1.2), colour="black")

volcano


###
### Beta difference
###

colnames(data.test.M) <- c("pvalue","qvalue","BetaDifference", "Median_ALZHEIMER", "Median_CONTROL")

ggplot(as.data.frame(data.test.M), aes(x=Beta-Difference)) + geom_histogram()
ggplot(as.data.frame(data.test.M[lista,]), aes(x=BetaDifference,y = ..density..)) + geom_histogram(colour = "darkgreen", fill = "white", binwidth = 0.1) + geom_density()
hist(data.test.M[,"BetaDifference"],breaks=200)

###
### Encode report
###

ranges <-  select_for_encode(meth450_data_original, data.test.M, lista)
report <- histones_report(ranges)


gene_list$threshold = as.factor(abs(gene_list$logFC) > 2 & gene_list$P.Value < 0.05/no_of_genes)


# meditacionq
clases <- as.factor(as.character(samps$grupo))

data.test <- sitetest_local(beta=b, grouplev=clases, gcase="ALZHEIMER", gcontrol="CONTROL",concov=concov,testmethod = testmethod, Padj=Padj, paired = paired)
data.test.M <- sitetest_local(beta=M, grouplev=clases, gcase="M", gcontrol="C",concov=concov,testmethod = testmethod, Padj=Padj, paired = paired)
save(data.test, file=paste("./data",projectName,"data.test.dat",sep="/"))
save(data.test.M, file=paste("./data",projectName,"data.test.M.dat",sep="/"))
load(paste("./data",projectName,"data.test.M.dat",sep="/"))
load(paste("./data",projectName,"data.test.dat",sep="/"))

hist(data.test[data.test[,"Adjust Pval"]<0.1,"Adjust Pval"])
hist(data.test[,"Beta-Difference"], breaks=100)
hist(data.test[abs(data.test[,"Beta-Difference"])>0.1,"Beta-Difference"])

lista <- names(which(abs(data.test[,"Beta-Difference"])>0.1 & data.test[,"Adjust Pval"]<0.1))
nsel <- length(lista)

hist(data.test.M[data.test.M[,"Adjust Pval"]<1,"Adjust Pval"])
lista.05 <- names(which(abs(data.test.M[,"Beta-Difference"])>0.1 & data.test.M[,"Adjust Pval"]<0.05))
lista <- names(which(abs(data.test.M[,"Beta-Difference"])>0.1 & data.test.M[,"Adjust Pval"]<0.1))

lista.jager = read.table("./results/alzheimer/lista_jager.txt")
lista431 = read.table("./results/alzheimer/lista_431.txt")
lista.anotada <- meth450_norm@featureData[lista,]
write.csv(pData(lista.anotada),file="./results/alzheimer/lista.anotada.csv")
i <- intersect(lista, lista431[,1])
lista.anotada.2 <- meth450_norm@featureData[i,]
write.csv(pData(lista.anotada.2),file="./results/alzheimer/lista.anotada.2.csv")

resultado <- cbind(data.test.M[lista.1,], data.test[lista.1,],b[lista.1,order(clases)], pData(lista.anotada))
write.csv(resultado,file="./results/resultado.csv")

lista.anotada <- meth450_norm@featureData[lista431,]
resultado <- cbind(data.test.M[lista431,], data.test[lista431,],b[lista431,order(clases)], pData(lista.anotada))
write.csv(resultado,file="./results/resultado431.csv")

interesultado.jager <- cbind(data.test.M[lista.jager,], b[lista.jager,order(clases)], pData(meth450_norm@featureData[lista.jager,]))
write.csv(resultado.jager,file="./results/resultado.jager.csv")

qcPlots(b[lista,],rownames(samps),clases)
b_heatmap.3(b,lista,clases)
MDSplot(t(b[lista,]), col.var = clases)
MDSplot(t(b), col.var = clases)
### using minfi DMP y DMR
data.rs <- RatioSet(Beta = b,annotation=c(array= "IlluminaHumanMethylation450k",  annotation = "ilmn12.hg19")); # create RatioSet
data.grs <- mapToGenome(data.rs); # create GenomicRatioSet
samples <- sampleNames(data.rs); # get sample names (ie column headers) from input file of beta values

pheno <- factor(phenoData(meth450_norm)$estadio == "EA_CONTROL", labels= c("DISEASE", "CONTROL"))
dmp <- dmpFinder(b, pheno = pheno, type = "categorical" )
head(dmp)

#j <- match(samples,make.names(sampleDescriptions$SampleID)); # get row indices from sample sheet that correspond to sample names (column headers) in beta value input file
#pheno <- sampleDescriptions$estadio[j] # extract corresponding phenotype (e.g. "tumor" or "normal") from Sample_Group column of Sample Sheet
designMatrix <- model.matrix(~ pheno); # specify design matrix for regression
dmr <- bumphunter(data.grs, design = designMatrix, cutoff = 0.05, B=100); # run bumphunter with B permutations
head(dmr$table)

tab<-dmr$table
Index=(tab[1,7]-3):(tab[1,8]+3)
matplot(pos[Index],meth[Index,,drop=TRUE],col=cols,pch=1,xlab="genomic location",ylab="Methylation",ylim=c(0,1))
plot(pos[Index],dmr$fitted[Index,1],xlab="genomic location",ylab="Methylation difference",ylim=c(-1,1))
abline(h=c(-0.1,0,.1),lty=2)

require(impute)
require(made4)

qcPlots(b[lista314[,1],],samples,clases)

lista431 = read.table("./results/lista_431.txt")
lista431 <- lista431[,1]
lista314 = read.table("./results/lista_314.txt")
lista314 <- lista314[,1]

m314 <- b[(lista314[,1]),]
m314 <- impute.knn(m314)$data
heatplot(m314)
b_heatmap(b,lista314[,1],clases)


library("ggplot2")
library("plyr")
library("reshape2")
library("scales")

ggplot(m314, aes(variable, Name)) + 
  geom_tile(aes(fill = rescale), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue") + 
  scale_x_discrete("", expand = c(0, 0)) + 
  scale_y_discrete("", expand = c(0, 0)) + 
  theme_grey(base_size = 9) + 
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_text(angle = 330, hjust = 0))


#### Filtramos filas #########################################################################

data_methylumi.pf <- data_methylumi[!is.na(as.array(meth450_data@featureData@data$MAPINFO)),]
data_methylumi.pf <- data_methylumi.pf[as.array(data_methylumi.pf@featureData@data$PROBE_SNPS_10)=="",]
data_methylumi.pf <- data_methylumi.pf[as.array(!(data_methylumi.pf@featureData@data$CHR)=="X"),]
data_methylumi.pf <- data_methylumi.pf[as.array(!(data_methylumi.pf@featureData@data$CHR)=="Y"),]

data_methylumi.pf<-pfilter(data_methylumi)

### Normalizamos
data_methylumi.pf.norm <- BMIQ(data_methylumi.pf)

data_methylumi.pf.norm <- dasen(data_methylumi.pf)

### Extraemos betas
betas <- betas(data_methylumi.pf.norm);
grouplev <- phenodataframe[,2]

data.test <- sitetest_local(beta=b, grouplev=phenodataframe[,2], gcase=gcase, gcontrol=gcontrol,concov=concov,testmethod = testmethod, Padj=Padj, paired = paired)
hist(data.test[data.test[,"Adjust Pval"]<0.1,"Adjust Pval"])
hist(data.test[,"Beta-Difference"], breaks=100)
hist(data.test[abs(data.test[,"Beta-Difference"])>0.1,"Beta-Difference"])

# Load package
library("FDb.InfiniumMethylation.hg19")

# To use the function please supply or use default settings:
# data:
# The DMR finder in the default setting assumes a R object with two columns: 
# one with Illumina identifiers and one with zeros and ones indicating DMP status 
# (Zero: non-DMP and One:DMP)
# To make the DMR finder also available for other types of data one can set 'illumina' to false and use a more generic input:
# Three columns: chromosome (for example "chr1"), genomic position and DMP status (Zero: non-DMP and One:DMP)
# chromosome = default is all chromosomes but a subset is also possible
# mismatches= maximum number of of allowed non-DMPs within DMRs
# icd = inter CpG distance
# The tDMRs can be found in the tDMR object (GRanges).
source("./code/DMRfinder.R", chdir=T)
data.DMR <- data.frame(id=rownames(b),dmp=as.numeric(rownames(b) %in% lista))

tDMRs = DMRfinder(data.DMR,chromosome=c(1:22), mismatches=3, icd=1000,illumina=TRUE)
save(tDMRs, file=paste("./data",projectName,"tDMRs.dat",sep="/"))
load(paste("./data",projectName,"tDMRs.dat",sep="/"))

library(Gviz)
data(cpgIslands)
atr <- AnnotationTrack(cpgIslands, name="CpG")
gtr <- GenomeAxisTrack()
itr <- IdeogramTrack(genome="mm9", chromosome="chr7")
data(geneModels)
grtr <- GeneRegionTrack(geneModels, name="Gene Model", showId=TRUE)
plotTracks(list(itr, gtr, atr, grtr))

library(GenomicRanges)
data(cpgIslands)
class(cpgIslands)
chr <- as.character(unique(seqnames(cpgIslands)))
gen <- genome(cpgIslands)
atrack <- AnnotationTrack(cpgIslands, name="CpG")
itrack <- IdeogramTrack(genome="hg19", chromosome="chr7")
plotTracks(list(itrack, gtrack, atrack))

gen = "hg19"
chr = "17"
start = 41277600
end = 41280000

biomartTrack <- BiomartGeneRegionTrack(genome="hg19", chromosome="chr17", start=41277600, end=41280000, name="ENSEMBL")



martfunc <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
ensfunc <- 
  getBM(c("ensembl_gene_id","ensembl_transcript_id","exon_chrom_start",
          "exon_chrom_end","strand","gene_biotype","external_gene_name"),
          filters = c("chromosome_name", "start", "end"),
          values = list(chr, start, end), mart=martfunc)
data_trackfunc <- AnnotationTrack(chr=chr,strand=ensfunc[,5],start=ensfunc[,3],end=ensfunc[,4],
                                   feature=ensfunc[,6],group=ensfunc[,1],id=ensfunc[,7],
                                   name = "genes ENSEMBL")

plotTracks(list(biomartTrack, data_trackfunc))
