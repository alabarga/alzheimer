################ DEFINICION EXPERIMENTO #################################

setwd("C:/PROYECTOS/BIO12_AL_007")

source("./code/require.r", chdir=T)
source("./code/auxiliares.r", chdir=T)
source("./code/params.r", chdir=T)
source("./code/graficos.r", chdir=T)
source("./code/ChAMP/quality control plots.r", chdir=T)
source('C:/PROYECTOS/BIO12_AL_007/code/Analisis/Locus-by-locus analysis/multiRegMultistat.R')

################ DEFINICION EXPERIMENTO #################################

source("./code/experiment2.r", chdir=T)

################ CARGA DATOS ############################################

load(paste("./data",projectName,"meth450_data_original.Rdata",sep="/"))
load(paste("./data",projectName,"meth450_data.dat",sep="/"))
load(paste("./data",projectName,"meth450_norm.dat",sep="/"))
load(paste("./data",projectName,"cor_b_edad.dat",sep="/"))
load(paste("./data",projectName,"results","res.dat",sep="/"))
################ CARGA DATOS ############################################

sampleDescriptions = read.table(SampleDescFileName,header=TRUE,sep=",")
sampleDescriptions$SampleLabel <- sampleDescriptions$SampleID

meth450_data_original <-methylumiR(filename=MethylDataFileName, sampleDescriptions=sampleDescriptions )
explore450(meth450_data_original, group="Source")

save(meth450_data_original, file=paste("./data",projectName,"meth450_data_original.Rdata",sep="/"))
load(paste("./data",projectName,"meth450_data_original.Rdata",sep="/"))

################ REMOVE BAD SAMPLES ############################################

meth450_data <- meth450_data_original[,!(pData(meth450_data_original)[,1] %in% remove_samples)]
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
meth450_data_6 <- meth450_data_5[!(rownames(meth450_data_5) %in% non_specific_probes[,1]),]

# 1328 sites were removed as beadcount <3 in 5 % of samples 406 sites having 1 % of samples with a detection p-value greater than 0.05 were removed
meth450_data_7 <- pfilter(meth450_data_6)

################
# noisy probes
################

extended_annotation = read.table("./data/1471-2164-15-51-s2.csv",sep=",",header=T)
probes_to_filter <- extended_annotation$probe[(extended_annotation$Flag.discard.keep. == "discard")]
meth450_data_8 <- meth450_data_7[!(rownames(meth450_data_7) %in% probes_to_filter),]

################
# probes left
################

nrow(meth450_data_8)
# 264902 left

explore450(meth450_data, group="Source")

meth450_data <- meth450_data_8
save(meth450_data, file=paste("./data",projectName,"meth450_data.dat",sep="/"))
load(paste("./data",projectName,"meth450_data.dat",sep="/"))

################################
# Remove intermediate results
################################

rm(meth450_data_8)
rm(meth450_data_7)
rm(meth450_data_6)
rm(meth450_data_5)
rm(meth450_data_4)
rm(meth450_data_3)
rm(meth450_data_2)

################
# Normalize
################

normalization <- "dasen"

meth450_norm <- normalize450(meth450_data, method=normalization)

explore450(meth450_norm, group="Source")

save(meth450_norm, file=paste("./data",projectName,"meth450_norm.dat",sep="/"))
load(paste("./data",projectName,"meth450_norm.dat",sep="/"))

################################
# Differential methylation
################################

samps <- pData(meth450_norm)
b <- betas(meth450_norm)
M <- exprs(meth450_norm)

clases <- as.factor(as.character(samps$Source))

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

#biopsia liquida
res <- add.betaDiff(res, b, clases, l1='Blood', l2='Brain')

################################
# Annotation
################################

commas <- function(x) {strsplit(as.character(x),",")[[1]]}
xgene <- function(x) {strsplit(as.character(x),";")[[1]][1]}
featureData <- meth450_norm@featureData@data
genome_build <- featureData[1,"GENOME_BUILD"]
annot <- featureData[,c("CHR",  "MAPINFO", "UCSC_REFGENE_NAME","UCSC_REFGENE_GROUP", "RELATION_TO_UCSC_CPG_ISLAND")]
annot$geneName <- sapply(as.character(annot$UCSC_REFGENE_NAME),xgene)

resn <- rownames(res)
annot.o <- annot[resn,]

top <- topTest(res, diff=0.7, pval= 0.01,  pval.var= "P.Value", diff.var="Beta.Difference")
selected <- rownames(top)

res$selected <- rownames(res) %in% selected
res$geneName <- annot.o$geneName

plotMDS(b[selected,], labels=samps$sampleID, col=as.integer(clases))

table_selected <- cbind(top,annot[selected,])
write.csv(table_selected, file=paste("./data",projectName,"results","table_selected_01022017.csv",sep="/"))

