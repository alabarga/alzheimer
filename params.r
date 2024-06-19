############sitetest/regionwrapper############################################################################
testmethod = "wilcox"            ### Other options of differential testing methods: "limma"/"pooled"/"satterthwaite" for the comparison between two group
concov = "OFF"                   ### If "ON", covariates is continuous variable
Padj = "BH"                      ### Options for multiple testing correction. The user can choose the methods provided by p.adjust function of R stat package 
indexmethod ="median"            ### Options for deriving an index of overall methylation value of each region. mean/median/tbrm: "tbrm" is Tukey's Biweight robust average 
paired = FALSE                   ### If ture, the differential test methods would change to the corresponding paired-test methods
##############################################################################################################

####################################output the differential sites#############################################
rawpcut = 0.05             ### cut off for raw pvalue 
adjustpcut = 1          ### cut off for adjusted pvalue
betadiffcut = 0         ### cut off for beta value difference
##############################################################################################################

#################Preprocessing:IMA.methy450PP ##############################################################
samplefilterdetectP = 1e-5   ### The cutoff for sample-level detection Pvalue
samplefilterperc = 0.75      ### The percent of loci with detection Pvalue less than "samplefilterdetectP" in each sample
sitefilterdetectP = 0.99     ### The cutoff for site-level detection Pvalue
sitefilterperc = 0.75         ### The percent of samples with detection Pvalue less than "sitefilterdetectP" for each site
na.omit = FALSE               ### Remove the sites containing missing beta value
XYchrom = FALSE               ### Remove the sites on chromosome X
peakcorrection = FALSE       ### If TRUE, peak correction is performed
normalization = FALSE        ### If TRUE, quantile normalization performed
transfm = FALSE              ### If FALSE, no transform is performed; if "arcsinsqr", arcsin square root transformation is performed; if "logit", logit transformation is performed;
locidiff = FALSE             ### If FALSE, don't filter sites by the difference of group beta value. Otherwise, remove the sites with beta value difference smaller than the specified value
locidiffgroup = c("CONTROL","MUESTRA") ### Specify which two groups are considered to check the loci difference (if "locidiff" is not true)
snpfilter = FALSE            ### If FALSE, keep the loci whose methylation level are measured by probes containing SNP(s) at/near the targeted CpG site; otherwise, filter out the list of SNP containing loci by specifying the snp file name and location
##############################################################################################################
### A list of SNP-containing probes (based on dbSNP v132) could be accessed by the command: snpfilter = system.file("extdata/snpsites.txt",package ="IMA")
##############################################################################################################
