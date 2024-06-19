### R code from vignette source 'vignettes/DMRcate/inst/doc/DMRcate.Rnw'

###################################################
### code chunk number 1: bioconductor (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("DMRcate")


###################################################
### code chunk number 2: libr
###################################################
library(DMRcate)


###################################################
### code chunk number 3: loaddata
###################################################
data(dmrcatedata)
myMs <- logit2(myBetas)


###################################################
### code chunk number 4: filter
###################################################
nrow(illuminaSNPs)
nrow(myMs)
myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05)
nrow(myMs.noSNPs)


###################################################
### code chunk number 5: annotate
###################################################
patient <- factor(sub("-.*", "", colnames(myMs)))
type <- factor(sub(".*-", "", colnames(myMs)))
design <- model.matrix(~patient + type) 
myannotation <- cpg.annotate(myMs.noSNPs, analysis.type="differential",
    design=design, coef=39)


###################################################
### code chunk number 6: dmrcate
###################################################
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)


###################################################
### code chunk number 7: plotting
###################################################
head(dmrcoutput$results)
DMR.plot(dmrcoutput=dmrcoutput, dmr=2, betas=myBetas, 
         phen.col=c(rep("orange", 38), rep("blue", 38)), 
         pch=16, toscale=TRUE, plotmedians=TRUE)


###################################################
### code chunk number 8: sessionInfo
###################################################
sessionInfo()


