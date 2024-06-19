library(minfi); # currently requires development version 1.9.11
library(IlluminaHumanMethylation450kanno.ilmn12.hg19); # Note: had to be installed manually
library(foreach); # used for parallel processing
library(doRNG); # is this actually used??
library(doParallel); # "backend" for foreach package
registerDoParallel(cores = 4);

basepath <- "C:/PROYECTOS/BIO12_AL_007/data"
GS.raw <- minfi:::read.GenomeStudio(file.path(basepath, "SampleMethFinalReport_nonorm.txt"))

betas <- as.matrix(GS.raw$beta)
#betas <- as.matrix(betas.frame); # convert data frame to matrix
data.rs <- RatioSet(Beta = betas,annotation=c(array= "IlluminaHumanMethylation450k",  annotation = "ilmn12.hg19")); # create RatioSet
data.grs <- mapToGenome(data.rs); # create GenomicRatioSet
samples <- sampleNames(data.rs); # get sample names (ie column headers) from input file of beta values

sampleDescriptions <- sampleDescriptions[!rownames(sampleDescriptions) %in% remove_samples,]
betas <- betas[,make.names(sampleDescriptions$SampleID)]
pheno <- factor(sampleDescriptions$estadio == "EA_CONTROL", labels= c("DISEASE", "CONTROL"))
dmp <- dmpFinder(betas, pheno = pheno, type = "categorical" )
head(dmp)
#j <- match(samples,make.names(sampleDescriptions$SampleID)); # get row indices from sample sheet that correspond to sample names (column headers) in beta value input file
#pheno <- sampleDescriptions$estadio[j] # extract corresponding phenotype (e.g. "tumor" or "normal") from Sample_Group column of Sample Sheet
designMatrix <- model.matrix(~ pheno); # specify design matrix for regression
dmr <- bumphunter(data.grs, design = designMatrix, cutoff = 0.05, B=100); # run bumphunter with B permutations
head(dmr$table)

