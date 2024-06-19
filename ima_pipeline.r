rm(list=ls())
options(stringAsfactors = FALSE)

############################################################################################################
############################################################################################################
########################Prepare all the parameters here#####################################################
############################################################################################################
############################################################################################################

#######################options in IMA.methy450R#############################################################
######################Load data#############################################################################
#libPaths = "/home/danwang/R/myR"                              ### Specify the location of your R library

MethyFileName = "../data/SampleMethFinalReport_nonorm.txt"
PhenoFileName = "../data/Pheno_without46ER.txt"

############################################################################################################

###########output file######################################################################################
siteleveltest = "../results/0_test/sitelevle.test.xls" ### Specify the path and name for the site-level testing result
list11excel = "../results/0_test/list11excel.xls"     ### Specify the path and name for region-level analysis results
list11Rdata = "../results/0_test/list11.Rdata"        ### Specify the path of Rdata file which stores the region-level analysis results
############################################################################################################

#################Preprocessing:IMA.methy450PP ##############################################################
samplefilterdetectP = 1e-5   ### The cutoff for sample-level detection Pvalue
samplefilterperc = 0.75      ### The percent of loci with detection Pvalue less than "samplefilterdetectP" in each sample
sitefilterdetectP = 0.05     ### The cutoff for site-level detection Pvalue
sitefilterperc = 0.5         ### The percent of samples with detection Pvalue less than "sitefilterdetectP" for each site
na.omit = TRUE               ### Remove the sites containing missing beta value
XYchrom = FALSE              ### Remove the sites on chromosome X
peakcorrection = FALSE       ### If TRUE, peak correction is performed
normalization = FALSE        ### If TRUE, quantile normalization performed
transfm = FALSE              ### If FALSE, no transform is performed; if "arcsinsqr", arcsin square root transformation is performed; if "logit", logit transformation is performed;
locidiff = FALSE             ### If FALSE, don't filter sites by the difference of group beta value. Otherwise, remove the sites with beta value difference smaller than the specified value
locidiffgroup = c("CONTROL","MUESTRA") ### Specify which two groups are considered to check the loci difference (if "locidiff" is not true)
snpfilter = FALSE            ### If FALSE, keep the loci whose methylation level are measured by probes containing SNP(s) at/near the targeted CpG site; otherwise, filter out the list of SNP containing loci by specifying the snp file name and location
                             ### A list of SNP-containing probes (based on dbSNP v132) could be accessed by the command: snpfilter = system.file("extdata/snpsites.txt",package ="IMA")
##############################################################################################################

############sitetest/regionwrapper############################################################################
testmethod = "wilcox"      ### Other options of differential testing methods: "wilcox"/"pooled"/"satterthwaite" for the comparison between two group
concov = "OFF"             ### If "ON", covariates is continuous variable
gcase = "MUESTRA"          ### Specify the case group index in the sample.txt file (if "concov" is "ON")
gcontrol = "CONTROL"       ### Specify the control group index in the sample.txt file (if "concov" is "ON")
Padj = "BH"                ### Options for multiple testing correction. The user can choose the methods provided by p.adjust function of R stat package 
indexmethod ="median"      ### Options for deriving an index of overall methylation value of each region. mean/median/tbrm: "tbrm" is Tukey's Biweight robust average 
paired = FALSE             ### If ture, the differential test methods would change to the corresponding paired-test methods
##############################################################################################################

####################################output the differential sites#############################################
rawpcut = 1             ### cut off for raw pvalue 
adjustpcut = 0.01       ### cut off for adjusted pvalue
betadiffcut = 0.1       ### cut off for beta value difference
##############################################################################################################

##############################################################################################################
###############################End of the parameter specification#############################################
##############################################################################################################


################################ Analysis Routes #############################################################
##############################################################################################################

.libPaths(libPaths) ## Specify your R library
library(IMA)        ## load the IMA package

data =IMA.methy450R(fileName = MethyFileName,columnGrepPattern=list(beta=".AVG_Beta",detectp=".Detection.Pval"),groupfile = PhenoFileName) ## load the data
dataf = IMA.methy450PP(data,na.omit = na.omit,normalization=normalization,peakcorrection = peakcorrection,transfm = transfm,samplefilterdetectP = samplefilterdetectP,samplefilterperc = samplefilterperc,sitefilterdetectP = sitefilterdetectP,locidiff = locidiff, locidiffgroup = locidiffgroup,XYchrom = XYchrom,snpfilter = snpfilter) ## QC filtering

sitetestALL = sitetest(dataf,gcase=gcase,gcontrol=gcontrol,concov=concov,testmethod = testmethod, Padj=Padj, paired = paired) ## site-level testing with the "BH" adjustment
filtrados = outputDMfunc(sitetestALL,rawpcut = rawpcut,adjustpcut =adjustpcut,betadiffcut = betadiffcut)

regionswrapper(dataf,indexmethod =indexmethod,gcase = gcase,gcontrol=gcontrol,testmethod = testmethod,Padj=Padj,concov = concov,list11excel=list11excel,list11Rdata=list11Rdata,rawpcut = rawpcut,adjustpcut = adjustpcut,betadiffcut = betadiffcut,paired = paired)   ## region-level testing for all 11 categories of annotated regions
############################### End of Analysis Routes #########################################################
################################################################################################################

hist(sitetest[,"Beta-Difference"],breaks=100)

