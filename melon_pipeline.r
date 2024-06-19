################ DEFINICION EXPERIMENTO #################################

setwd("C:/PROYECTOS/BIO12_AL_007")

source("./code/require.r", chdir=T)
source("./code/auxiliares.r", chdir=T)
source("./code/params.r", chdir=T)
source("./code/graficos.r", chdir=T)
source("./code/ChAMP/quality control plots.r", chdir=T)
source('C:/PROYECTOS/BIO12_AL_007/code/Analisis/Locus-by-locus analysis/multiRegMultistat.R')

################ DEFINICION EXPERIMENTO #################################

source("./code/experiment4.r", chdir=T)

################ CARGA DATOS ############################################

load(paste("./data",projectName,"meth450_data_original.Rdata",sep="/"))
load(paste("./data",projectName,"meth450_data.dat",sep="/"))
load(paste("./data",projectName,"meth450_norm.dat",sep="/"))
load(paste("./data",projectName,"cor_b_edad.dat",sep="/"))
load(paste("./data",projectName,"results","res.dat",sep="/"))
################ CARGA DATOS ############################################

#alzheimer
sampleDescriptions = read.table(SampleDescFileName,header=TRUE,sep=",")
colnames(sampleDescriptions) <- c("SampleID","estadio","sexo","edad")
rownames(sampleDescriptions) <- sampleDescriptions$SampleID

# meditacion
sampleDescriptions = read.table(SampleDescFileName,header=TRUE,sep=",")
sampleDescriptions$SampleID <- rownames(sampleDescriptions)
# colnames(sampleDescriptions) <- c("status","SampleID")

# liquida
sampleDescriptions = read.table(SampleDescFileName,header=TRUE,sep=",")
sampleDescriptions$SampleLabel <- sampleDescriptions$SampleID

meth450_data_original <-methylumiR(filename=MethylDataFileName, sampleDescriptions=sampleDescriptions )
explore450(meth450_data_original, group="Source")

# sapply(colnames(meth450_data_original), function(x){print predictSex(meth450_data_original,x)})
save(meth450_data_original, file=paste("./data",projectName,"meth450_data_original.Rdata",sep="/"))
load(paste("./data",projectName,"meth450_data_original.Rdata",sep="/"))

################ REMOVE BAD SAMPLES ############################################

meth450_data <- meth450_data_original
explore450(meth450_data, group="Lote")

#meditacion pData(meth450_data_original)[,2]
meth450_data <- meth450_data_original[,!(pData(meth450_data_original)[,2] %in% remove_samples)]
explore450(meth450_data, group="Source")

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


# sum(meth450_data@featureData@data$PROBE_SNPS_10!="")
#[1] 36535
# sum(meth450_data@featureData@data$PROBE_SNPS!="")
#[1] 59234

# sum(!is.na(SNPs.138$Probe_rs))
#[1] 87317

meth450_data_2 <- meth450_data

################
# PROBE_SNPS_10
################

meth450_data_3 <- meth450_data_2[as.array(meth450_data_2@featureData@data$PROBE_SNPS_10)=="",]

################
# chrX
################

meth450_data_4 <- meth450_data_3[as.array(!(meth450_data_3@featureData@data$CHR)=="X"),]

################
# chrY
################

meth450_data_5 <- meth450_data_4[as.array(!(meth450_data_4@featureData@data$CHR)=="Y"),]

################
# non_specific_probes
################

non_specific_probes = read.table("./data/lista_29233.txt")
sum(rownames(meth450_data_5) %in% non_specific_probes[,1])
meth450_data_6 <- meth450_data_5[!(rownames(meth450_data_5) %in% non_specific_probes[,1]),]

# 1328 sites were removed as beadcount <3 in 5 % of samples 406 sites having 1 % of samples with a detection p-value greater than 0.05 were removed
meth450_data_7 <- pfilter(meth450_data_6)

################
# noisy probes
################

extended_annotation = read.table("./data/1471-2164-15-51-s2.csv",sep=",",header=T)
probes_to_filter <- extended_annotation$probe[(extended_annotation$Flag.discard.keep. == "discard")]

sum(rownames(meth450_data_7) %in% probes_to_filter)

meth450_data_8 <- meth450_data_7[!(rownames(meth450_data_7) %in% probes_to_filter),]

################
# probes left
################

nrow(meth450_data_8)
# 265373 left

explore450(meth450_data, group="Source")

save(meth450_data, file=paste("./data",projectName,"meth450_data.dat",sep="/"))
load(paste("./data",projectName,"meth450_data.dat",sep="/"))

meth450_data <- meth450_data_8

################
# Normalize
################

normalization <- "dasen"

meth450_norm <- normalize450(meth450_data, method=normalization)

explore450(meth450_norm, group="Source")

save(meth450_norm, file=paste("./data",projectName,"meth450_norm.dat",sep="/"))
load(paste("./data",projectName,"meth450_norm.dat",sep="/"))

# meth450_norm <- meth450_norm[,!(pData(meth450_norm)[,2] %in% remove_samples)]

M <- exprs(meth450_data)
M <- removeBatchEffect(M,batch=samps$Lote, design=design)
b <- 2^M / (2^M + 1)

save(meth450_norm, file=paste("./data",projectName,"meth450_norm.dat",sep="/"))
load(paste("./data",projectName,"meth450_norm.dat",sep="/"))

################################
# Differential methylation
################################

samps <- pData(meth450_norm)
b <- betas(meth450_norm)
M <- exprs(meth450_norm)
M <- getM(meth450_norm) # for signature '"MethylSet"'

modBatch = model.matrix(~as.factor(status) + as.factor(Lote),data=samps)
mod0Batch = model.matrix(~as.factor(Lote),data=samps)
pValuesBatch = f.pvalue(M,modBatch,mod0Batch)
qValuesBatch = p.adjust(pValuesBatch,method="BH")

M2 <- removeBatchEffect(M,batch=samps$Lote, design=design)
Mbatch <- M
M <- M2

################################
# DMRcate
################################

myMs.noSNPs <- rmSNPandCH(M, dist=2, mafcut=0.05)
print(paste(nrow(M) - nrow(myMs.noSNPs), "sondas eliminadas",sep=" "))

#alzheimer
clases = c("CONTROL", "ALZHEIMER")[as.factor(samps$estadio !="EA_CONTROL")]
clases <-as.factor(clases)

#meditacion
clases <- as.factor(samps$status)

design <- model.matrix(~clases)

myannotation <- cpg.annotate(myMs.noSNPs, analysis.type="differential",
                             design=design, coef=2)

dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
nrow(dmrcoutput$results)

View(head(dmrcoutput$results[(dmrcoutput$results$maxbetafc)>0.05,],20))
View(head(dmrcoutput$results[order(-abs(dmrcoutput$results$maxbetafc)),],20))

dmrs <- subset(dmrcoutput$results,dmrcoutput$results$no.probes>1 & abs(dmrcoutput$results$maxbetafc) > 0.05)
nrow(dmrs)
View(head(dmrs[order(-abs(dmrs$maxbetafc)),],20))
DMR.plot(dmrcoutput=dmrcoutput, dmr=which(rownames(dmrcoutput$results)==1363), betas=b, 
         phen.col=c("blue","red")[clases], 
         pch=16, toscale=TRUE, plotmedians=TRUE)

# MYL5
chr= "chr4"
startcpg = 675137 
stopcpg = 675827 

# KBTBD11
chr="chr8"
startcpg = 1954777 
stopcpg = 1955196

#NR4A2
chr="chr2"
startcpg = 157182707
stopcpg = 157187500

DMR.plot.region(dmrcoutput=dmrcoutput, chr, startcpg, stopcpg, betas=b, 
         phen.col=c("blue","red")[clases], 
         pch=16, toscale=TRUE, plotmedians=TRUE)

save(dmrcoutput, file=paste("./data",projectName,"dmrcoutput.Rdata",sep="/"))
gene.list <- unlist(sapply(unique(dmrs$gene_assoc),commas))
names(gene.list)<-NULL
write.table(gene.list,file=paste("./data",projectName,"genes.txt",sep="/"), row.names=FALSE,col.names=FALSE,quote=FALSE)

load(paste("./data",projectName,"dmrcoutput.Rdata",sep="/"))

################################
# Annotation
################################

commas <- function(x) {strsplit(as.character(x),",")[[1]]}
xgene <- function(x) {strsplit(as.character(x),";")[[1]][1]}
featureData <- meth450_norm@featureData@data
genome_build <- featureData[1,"GENOME_BUILD"]
annot <- featureData[,c("CHR",  "MAPINFO", "UCSC_REFGENE_NAME","UCSC_REFGENE_GROUP", "RELATION_TO_UCSC_CPG_ISLAND")]
annot$geneName <- sapply(as.character(annot$UCSC_REFGENE_NAME),xgene)

genes_gr <- GRanges(seqnames = paste("chr",as.character(res$chromosome_name),sep=""), ranges = IRanges(start = res$start_position - 5000, end = res$end_position + 5000), names = res$hgnc_symbol)

mi <- meth450_data@featureData@data
mi <- mi[,c("TargetID", "CHR","MAPINFO", "UCSC_REFGENE_NAME")]
probes_gr <- GRanges(seqnames = paste("chr",as.character(mi$CHR),sep=""), ranges = IRanges(start = mi$MAPINFO, end = mi$MAPINFO + 1), names = mi$TargetID)


## samps$estadio <- sampleDescriptions[rownames(samps),"estadio"] ## correccion

# alzheimer
clases = c("CONTROL","ALZHEIMER")[as.factor(samps$estadio!="EA_CONTROL")]
clases <-as.factor(clases)

# meditacion
clases <- as.factor(as.character(samps$status))

# biopsia liquida
clases <- as.factor(as.character(samps$Source))

##
## correlation with other variables
##

cor_M_edad <- apply(M,1,cor, y =samps$edad)
hist(cor_M_edad)

amieloide <- read.csv('data/alzheimer/AMILOIDE_081014.csv')
head(M[,amieloide$SampleName])


cor_M_amieloide <- apply(M[,amieloide$SampleName],1,cor, y =amieloide$Mean.total.area)
hist(cor_M_amieloide)


cor_M_amieloide <- apply(M[,amieloide$SampleName],1, function(x,y) {
  a <- cor.test(x, y,
                method="pearson")[c("estimate", "p.value", "conf.int")]
  out <- data.frame(cor = a$estimate, p.value = a$p.value, lowerCI = a$conf.int[1], upperCI = a$conf.int[2])
}, y =amieloide$Mean.total.area)
cm <- do.call("rbind", cor_M_amieloide)

cor_b_edad <- apply(b,1, function(x,y) {
  a <- cor.test(x, y,
                method="pearson")[c("estimate", "p.value", "conf.int")]
  out <- data.frame(cor = a$estimate, p.value = a$p.value, lowerCI = a$conf.int[1], upperCI = a$conf.int[2])
}, y=samps$edad)
cor_b_edad <- do.call("rbind", cor_b_edad)

colnames(cor_b_edad) <- paste('cor.b.edad',colnames(cor_b_edad),sep='.')
save(cor_b_edad, file=paste("./data",projectName,"cor_b_edad.dat",sep="/"))
load(paste("./data",projectName,"cor_b_edad.dat",sep="/"))

profilePlot(b,names(cor_M_amieloide),clases)

###################################################
### Wilcox 
###################################################

data.test <- sitetest_local(beta=b, grouplev=clases, gcase="ALZHEIMER", gcontrol="CONTROL", concov=concov,testmethod = testmethod, Padj=Padj, paired = paired)
top2 <- topTest(data.test)
plotMDS(b[rownames(top2),], labels=samps$sampleID, col=as.integer(clases))

###################################################
### code chunk number 7: mdsplot
###################################################

par(mfrow=c(1,1))
plotMDS(M, labels=samps$sampleID, col=as.integer(clases))
legend("topleft",legend=c("Blood","Brain"),pch=16,cex=1.2,col=1:2)

###################################################
### code chunk number 8: design
###################################################

id <- samps$sampleID
id <- factor(gsub(" ", "", as.character(id)))

###################################################
### code chunk number 9: diffmeth
###################################################
design <- model.matrix(~clases)
fit <- lmFit(M,design)
fit <- eBayes(fit)

###################################################
### code chunk number 10: results
###################################################
summary(decideTests(fit))
res<-topTable(fit,coef=2, number=nrow(b), adjust.method="BH", p.value=1)
head(res)

#meditacion
res <- add.betaDiff(res, b, clases, l1='M', l2='C')

#alzheimer
res <- add.betaDiff(res, b, clases, l1='ALZHEIMER', l2='CONTROL')

#biopsia liquida
res <- add.betaDiff(res, b, clases, l1='Blood', l2='Brain')

resn <- rownames(res)
annot.o <- annot[resn,]

top <- topTest(res, diff=0.1, pval = 0.1)
top <- topTest(res, diff=0.7, pval= 0.01,  pval.var= "P.Value", diff.var="Beta.Difference")
selected <- rownames(top)
selected <- selected[1:255]
top<- top[1:255,]

# meditacion
load('data/meditacion/results/table_selected.dat')
selected <- rownames(table_selected)

res$selected <- rownames(res) %in% selected
res$geneName <- annot.o$geneName

table_selected <- cbind(top,annot[selected,])
write.csv(table_selected, file=paste("./data",projectName,"results","table_selected_26012017.csv",sep="/"))

lista2 <- read.table(file="data/meditacion/results/selectedGenes.txt")

save(res, file=paste("./data",projectName,"results","res.dat",sep="/"))
load(paste("./data",projectName,"results","res.dat",sep="/"))

dmrs <- read.csv('data/meditacion/results/dmrs_30042015.csv',sep="\t", stringsAsFactors=FALSE)
View(dmrs)


library(stringr)
dmrs_coords <- as.data.frame(str_match(dmrs$hg19coord, "^(.*):(.*)-(.*)$")[,-1], stringsAsFactors=FALSE)
colnames(dmrs_coords) <- c("chr","start","end")
dmrs_gr <- GRanges(seqnames = dmrs_coords$chr, ranges = IRanges(start = as.numeric(dmrs_coords$start), end = as.numeric(dmrs_coords$end)), names = res$dmrs$hg19coord)
overlaps<- findOverlaps(probes_gr, dmrs_gr, type="within")
cpgs <- (probes_gr[queryHits(overlaps),])$names

volcano450(res,  adjust.p = 0.05, diff.thresh=0.05, pval.var="P.Value", diff.var="Beta.Difference", names=TRUE, lista = cpgs)


#res[subjectHits(overlaps),]


###################################################
### code chunk number 12: plot selected
###################################################

plotMDS(b[selected,], labels=samps$sampleID, col=as.integer(clases))
b_heatmap.3(b, selected, clases)
volcano450(res,  adjust.p = 0.0005, diff.thresh=0.05, pval.var="P.Value", diff.var="Beta.Difference")
beta_FC(res, pval.var= "P.Value", adjust.p = 0.0005)

###################################################
### code chunk number 13: BACA functional analysis
###################################################
# http://cran.r-project.org/web/packages/BACA/vignettes/BACA.html

pos<-unique(top$geneName[which(top$Beta.Difference>0)])
neg<-unique(top$geneName[which(top$Beta.Difference<0)])

library("biomaRt")
ensembl=useMart("ensembl")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
getBM(attributes=c('affy_hg_u133_plus_2', 'entrezgene'), filters = 'affy_hg_u133_plus_2', values = affyids, mart = ensembl)
filters = listFilters(ensembl)
View(filters)

entrez.pos <- getBM(attributes=c('hgnc_symbol', 'entrezgene'), filters = 'hgnc_symbol', values = pos, mart = ensembl)
entrez.neg <- getBM(attributes=c('hgnc_symbol', 'entrezgene'), filters = 'hgnc_symbol', values = neg, mart = ensembl)

diff.genes <- list()
diff.genes$pos <- entrez.pos$entrezgene
diff.genes$neg <- entrez.neg$entrezgene


library(RDAVIDWebService)
library(ggplot2)
library(BACA)
david.obj <- DAVIDWebService$new(email="vittorio.fortino@ttl.fi")

result.annot <- DAVIDsearch(diff.genes, david.user = "vittorio.fortino@ttl.fi", idType="ENTREZ_GENE_ID", annotation="GOTERM_BP_3")

bbplot.annot <- BBplot(result.annot, max.pval = 0.2, min.ngenes = 2, 
                      name.com = c("Alzheimer's disease"), 
                      labels = c("down", "up"), colors = c("#009E73", "red"), 
                      title = "BBplot - KEGG")

plot(bbplot.annot)
#######################################################
### code chunk number 12: build and save table result
#######################################################

res <- cbind(res,b[rownames(res),])
res <- cbind(res,annot[rownames(res),])


top <- res #[which(res$color %in% c("D+","D-")),]
selected <- rownames(top)
table_selected <- cbind(top,annot[selected,])

table_selected <- cbind(table_selected,cor_b_edad[selected,c('cor.b.edad.cor','cor.b.edad.p.value')])

write.csv(table_selected, file=paste("./data",projectName,"results","table_total_17032015.csv",sep="/"))
save(table_selected, file=paste("./data",projectName,"results","table_selected.dat",sep="/"))
load(paste("./data",projectName,"results","table_selected.dat",sep="/"))

table_selected <- top
table_selected <- read.csv(file=paste("./data",projectName,"results","table_selected_17032015.csv",sep="/"), row.names=T)
selected <- rownames(table_selected)  
table_selected <- cbind(table_selected,annot[selected,])
table_selected <- cbind(table_selected,b[selected,])
                        
###################################################
### code chunk number 12: diffvar
###################################################
fitvar <- varFit(M, design = design, coef = c(1,2))
summary(decideTests(fitvar))
topDV <- topVar(fitvar, coef=2)
selected <- rownames(topDV)



###
### volcano
###

gene_list <- topTable(fit,coef=2,number=10000)
gene_list <- res
gene_list <- as.data.frame(top2)
gene_list <- as.data.frame(data.test)

gene_list$geneName <- sapply(meth450_data_original[rownames(gene_list)]@featureData@data$UCSC_REFGENE_NAME,xgene)
nombres <- sapply(meth450_data_original[rownames(gene_list)]@featureData@data$UCSC_REFGENE_NAME,xgene)
gene_list$geneName <- nombres

cn <- colnames(gene_list)
cn
cn[1]
cn[1] <- "p.value"
cn[2] <- "adj.P.Val"
cn[3] <- "beta.diff"
colnames(gene_list) <- cn


#gene_list$threshold = as.factor(abs(gene_list$logFC) > 0.1 & gene_list$adj.P.Val < 0.1)
gene_list$threshold = as.factor(abs(gene_list[,"beta.diff"]) > 0.1 & gene_list["adj.P.Val"] < 0.1)

head(gene_list)


res$selected = as.factor(abs(top[,"Beta.Difference"]) > 0.1 & top[,"adj.P.Val"] < 0.1)
res$geneName <- sapply(as.character(annot$UCSC_REFGENE_NAME)[rownames(res)],xgene)

res$selected <- rownames(res) %in% selected
res$geneName <- annot$geneName[rownames(res)]


res$color[(res[,"Beta.Difference"]) > 0.1 & res[,"adj.P.Val"] < 0.2] = "D+"
res$color[(res[,"Beta.Difference"]) < - 0.1 & res[,"adj.P.Val"] < 0.2] = "D-"
res$color[abs(res[,"Beta.Difference"]) < 0.1 | res[,"adj.P.Val"] > 0.2] = "NoDM"


volcano = ggplot(data=res, aes(x=Beta.Difference, y=-log10(adj.P.Val), colour=color)) +
  scale_colour_manual(values=c("#56B4E9", "#FF9999", "#888888")) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  xlim(c(-0.2, 0.2)) + ylim(c(0, 3)) +
  xlab("beta difference") + ylab("-log10 p-value") +
  geom_text(data=subset(res, as.logical(res$color == "D+")), aes(x=Beta.Difference, y=-log10(adj.P.Val),
                                                                 label=geneName, size=1.2), colour="#FF9999") +
  geom_text(data=subset(res, as.logical(res$color == "D-")), aes(x=Beta.Difference, y=-log10(adj.P.Val),
                                                                 label=geneName, size=1.2), colour="#56B4E9")
volcano


lista <- rownames(gene_list[which(gene_list$threshold),])

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
gtrack <- GeneRegionTrack(geneModels, name="Gene Model", showId=TRUE)
plotTracks(list(itrack, gtrack, atrack))

gen = "hg19"
chr = "17"
start = 41277600
end = 41280000

library(biomaRt)

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
