library(wateRmelon)
library(IMA)
library(reshape)
library(minfi)
library(RColorBrewer)
library(limma)
library(missMethyl)
library(DMRcate)
library(BiocInstaller)
library(Gviz)
library(ggplot2)
library(plyr)
library('FDb.InfiniumMethylation.hg19')

hm450 <- get450k()
probenames <- c("cg16392865", "cg00395291", "cg09310185", "cg21749424")
probes <- hm450[probenames]


#biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other
names(anno)
