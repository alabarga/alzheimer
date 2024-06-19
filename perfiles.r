setwd("C:/PROYECTOS/BIO12_AL_007")
source("./code/require.r", chdir=T)
source("./code/auxiliares.r", chdir=T)
source("./code/params.r", chdir=T)
source("./code/graficos.r", chdir=T)
source("./code/ChAMP/quality control plots.r", chdir=T)
source('C:/PROYECTOS/BIO12_AL_007/code/Analisis/Locus-by-locus analysis/multiRegMultistat.R')
################ DEFINICION EXPERIMENTO #################################
source("./code/experiment.r", chdir=T)
################ CARGA DATOS ############################################
load(paste("./data",projectName,"meth450_data_original.Rdata",sep="/"))
load(paste("./data",projectName,"meth450_data.dat",sep="/"))
load(paste("./data",projectName,"meth450_norm.dat",sep="/"))
samps <- pData(meth450_norm)
b <- betas(meth450_norm)
M <- exprs(meth450_norm)


commas <- function(x) {strsplit(as.character(x),",")[[1]]}
xgene <- function(x) {strsplit(as.character(x),";")[[1]][1]}
featureData <- meth450_norm@featureData@data
genome_build <- featureData[1,"GENOME_BUILD"]
annot <- featureData[,c("CHR",  "MAPINFO", "UCSC_REFGENE_NAME","UCSC_REFGENE_GROUP", "RELATION_TO_UCSC_CPG_ISLAND")]
annot$geneName <- sapply(as.character(annot$UCSC_REFGENE_NAME),xgene)

clases = c("CONTROL","ALZHEIMER")[as.factor(samps$estadio!="EA_CONTROL")]
clases <-as.factor(clases)
design <- model.matrix(~clases)

fit <- lmFit(M,design)
fit <- eBayes(fit)

myMs.noSNPs <- rmSNPandCH(M, dist=2, mafcut=0.05)
print(paste(nrow(M) - nrow(myMs.noSNPs), "sondas eliminadas",sep=" "))
myannotation <- cpg.annotate(myMs.noSNPs, analysis.type="differential",
design=design, coef=2)
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
nrow(dmrcoutput$results)
my_dmrs = dmrcoutput$results[abs(dmrcoutput$results$maxbetafc)>0.1,]

c <- colsplit(my_dmrs$hg19coord,":",c("chr","start"))
x0 <- (data.frame(chr=c$chr, colsplit(c$start,"-",c("start","stop"))))
gr <-with(x0, GRanges(chr, IRanges(start=start, end=stop)))

x1 <- data.frame(chr=paste("chr",annot$CHR,sep=""),start=annot$MAPINFO,stop=annot$MAPINFO+1)
ir <- with(x1, GRanges(chr, IRanges(start=start, end=stop)))
hits = findOverlaps(ir, gr)

ov <- rownames(annot)[hits@queryHits]

a <- data.frame(M=M[ov[1],],estadio=samps$estadio)
ggplot(a,aes(x=estadio, y=M, fill=estadio)) + geom_boxplot() + ggtitle(ov[1])
Mav <- aggregate(M[ov,],by=list(genes=hits@subjectHits),FUN=mean)

Mav$gene <- my_dmrs$gene_assoc
colnames(Mav)<- c('ID',paste(colnames(M),samps$estadio,sep="|"),"gene")
Mmelt <- melt(Mav,id.vars=c("ID","gene"))
Mmelt2 <- with(Mmelt,data.frame(ID,gene,colsplit(variable,"\\|",c("sample","estadio")),M=value))

xx <-dcast(Mmelt2,ID ~ estadio, fun.aggregate = mean, na.rm=TRUE)
xx$gene <- my_dmrs$gene_assoc

x <- melt(xx,id.vars=c("ID","gene"))
p <- ggplot(data = x, aes(x = variable, y = value,group=ID)) + geom_line() 

p + geom_text(data = x[x$variable == "EA3",], aes(label = gene), hjust = 0.7, vjust = 1)
p + geom_text(data = x[x$variable == "EA3",], aes(label = ID), hjust = 0.7, vjust = 1)


a <- data.frame(M=M[ov[i],],estadio=samps$estadio)
ggplot(xx[],aes(x=estadio, y=M, fill=estadio)) + geom_boxplot() + ggtitle(gene)

i <- 19
ggplot(data=Mmelt2[Mmelt2$ID==i,],aes(x=estadio, y=M, fill=estadio)) + geom_boxplot() +ggtitle(my_dmrs$gene_assoc[i])
