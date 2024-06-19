##########################################
############## tDMR analysis #############
##########################################
##########################################
### Script written by Roderick Slieker ###
### Molecular Epidemiology             ###
### Leiden University Medical Center   ###
### r.c.slieker@lumc.nl                ###
##########################################

# Load package
library("FDb.InfiniumMethylation.hg19")

# Change working directory to a writable directory!

source("DMRfinder.R")

# The slots of the function include:
# data:
    # The DMR finder in the default setting assumes a R object with two columns: 
    # one with Illumina identifiers and one with zeros and ones indicating DMP status 
    # (Zero: non-DMP and One:DMP)
    # To make the DMR finder also available for other types of data one can set 'illumina' to false and use a more generic input:
    # Three columns: chromosome (for example "chr1"), genomic position and DMP status (Zero: non-DMP and One:DMP)

# mismatches= maximum number of of allowed non-DMPs within DMRs
# icd = inter CpG distance
# The tDMRs can be found in the tDMR object (GRanges).

tDMRs = DMRfinder(data,mismatches=3, icd=1000,illumina=TRUE)