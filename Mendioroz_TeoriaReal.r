library(wateRmelon)
library(IMA)
library(reshape)

###########porcentaje NA's######################################################################################
#porcentajeNAenRegion = 10
#porcentajeNAenGrupo = 5
############################################################################################################

source("mendioroz.R")
source("auxiliares.R")

MethyFileName_Ctrl_Bkg = "./data/alzheimer/SampleMethFinalReport_ctrl-bkg.txt"
MethyFileName_NoNorm = "./data/alzheimer/SampleMethFinalReport_nonorm.txt"
PhenoFileName = "./data/alzheimer/Pheno_lumi.txt"
PhenoFileNameIMA = "./data/alzheimer/Pheno_without46ER.txt"
#PhenoFileName = "./data/Pheno_lumi_Fases.txt"
PhenoFileNameIMA = "./data/alzheimer/Pheno_without46ER_soloFases.txt"
PhenoFileNameAGE = "./data/alzheimer/Pheno_without46ER_age.txt"


################ CARGA DATOS ############################################
#load("./fullannotInd_all.rda") 
all_snps<-read.delim("snpsites_240000_integer.txt",header=TRUE)
pheno_ages = read.table(PhenoFileNameAGE,sep="\t", header=TRUE)
ages <- pheno_ages$age
phenodataframe = read.table(PhenoFileName,sep="\t",header=TRUE)

################ CARGA DATOS ############################################

data_methylumi<-methylumiR(filename= MethyFileName_NoNorm, sampleDescriptions=phenodataframe )

#### Filtramos filas #########################################################################

data_methylumi_filtered <- data_methylumi[!is.na(as.array(data_methylumi@featureData@data$MAPINFO)),]
data_methylumi_filtered <- data_methylumi_filtered[as.array(data_methylumi_filtered@featureData@data$PROBE_SNPS_10)=="",]
data_methylumi_filtered <- data_methylumi_filtered[as.array(!(data_methylumi_filtered@featureData@data$CHR)=="X"),]
data_methylumi_filtered <- data_methylumi_filtered[as.array(!(data_methylumi_filtered@featureData@data$CHR)=="Y"),]

#### Filtramos muestras erroneas o que no pertenecen al estudio
data_methylumi_filtered <- data_methylumi_filtered[,pData(data_methylumi_filtered)[,1]!= "46 ER"]

### Normalizamos
bmiq_data <- BMIQ(data_methylumi_filtered)

### Extraemos betas
bmiq_betas <- betas(bmiq_data);

### aux cleanup

data_methylumi_backup <- NULL
data_methylumi_backup <- data_methylumi
bmiq_betas_samesamples <- bmiq_betas
rm(bmiq_betas)
rm(bmiq_data)

### Starting Quality Control... (from lumi data)

# assayData(data_methylumi_filtered)$pvals
# assayData(data_methylumi_filtered)$betas
# sampleNames(phenoData(data_methylumi_filtered))

### Starting Quality Control... (from IMA data)
# data@detectP
# data@bmatrix
# data@groupinfo[, 1]

qcPlots <- function(betas,pvalues,samples){
  
  require(bioDist)
  eset = na.omit(betas)
  
  hc1 <- hclust(cor.dist(t(eset)), method = "average")
  
  plot(hc1, samples, xlab = "Sample", main = "Clustering samples by all the CpG loci ", 
       lwd = 2, font.axis = 2, font.lab = 2)
  boxplot(betas, ylab = "beta Value")
  avgPval = apply(pvalues, 2, function(x) {
    sum(x >= 1e-05) * 100/length(x)
  })
  barplot(avgPval, ylab = "% of detect pvalue >1e-5")

}

## Extract annotation

annotation <- pData(featureData(data_methylumi))
annot_id<-cbind(Id=row.names(annotation),annotation);
annot_id_columnsSelected <- annot_id[,c(1,13,14,19,20,23)]
rm(annot_id)

################################################# IMA Pipeline
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Filtramos igual que anteriormente
#CUIDADO!!!! Los Indices de regiones no se han podido filtrar. Con la subseleccion solo se puede analizar el site analysis por el momento!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

data =IMA.methy450R(fileName = MethyFileName_Ctrl_Bkg,columnGrepPattern=list(beta=".AVG_Beta",detectp=".Detection.Pval"),groupfile = PhenoFileNameIMA)
# file.rename("./QC.pdf","./QC_Ctrl_Bkg.pdf")
dataf = IMA.methy450PP(data,na.omit = na.omit,normalization=normalization,peakcorrection = peakcorrection,transfm = transfm,samplefilterdetectP = samplefilterdetectP,samplefilterperc = samplefilterperc,sitefilterdetectP = sitefilterdetectP,locidiff = locidiff, locidiffgroup = locidiffgroup,XYchrom = XYchrom,snpfilter = snpfilter) ## QC filtering
rm(data)

dataf2 <- NULL;
dataf2 <- dataf;
#SNP10
for(sloti in slotNames(dataf)[c(1,2,3)]){
  sloti
  slot(dataf2, sloti) <- slot(dataf,sloti)[is.na(annot_id_columnsSelected[,5]) & !is.na(annot_id_columnsSelected[,3]) & annot_id_columnsSelected[,2] != "X" & annot_id_columnsSelected[,2] != "Y"]

}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
beta <- dataf2@bmatrix;
betas_converted <- convertBMIQ2IMA(bmiq_betas_samesamples,beta)
dataf_backup <- NULL
dataf_backup <- dataf
#dataf <- dataf2;
dataf@bmatrix <- betas_converted

rm(betas_id)
rm(annot_id_columnsSelected)
rm(dataf_backup)
rm(betas_converted)
rm(bmiq_betas_samesamples)
rm(beta)
rm(dataf2)

################################################# 

betas_id<-cbind(Id=row.names(dataf@bmatrix),dataf@bmatrix)
betas_id_annot<-merge(betas_id,annot_id_columnsSelected,by.x="Id",by.y="Id")

################################################# 

#md <- cmdscale(dist(t(betas_id_annot[betas_id_annot$CHR == "X",seq(2,38)])))
#plot(md,pch=c('F','M')[phenodataframe_sex$sex],col=c("red","blue")[phenodataframe_sex$sex])

################# SITE ANALYSIS PVALUE 0.05#######################################

############sitetest/regionwrapper############################################################################
testmethod = "wilcox"       ### Other options of differential testing methods: "limma"/"pooled"/"satterthwaite" for the comparison between two group
concov = "OFF"             ### If "ON", covariates is continuous variable
gcase = "MUESTRA"               ### Specify the case group index in the sample.txt file (if "concov" is "ON")
gcontrol = "CONTROL"            ### Specify the control group index in the sample.txt file (if "concov" is "ON")
#gcase = "EA1"               ### Specify the case group index in the sample.txt file (if "concov" is "ON")
#gcontrol = "EA_CONTROL"            ### Specify the control group index in the sample.txt file (if "concov" is "ON")
Padj = "BH"                ### Options for multiple testing correction. The user can choose the methods provided by p.adjust function of R stat package 
indexmethod ="median"        ### Options for deriving an index of overall methylation value of each region. mean/median/tbrm: "tbrm" is Tukey's Biweight robust average 
paired = FALSE             ### If ture, the differential test methods would change to the corresponding paired-test methods
##############################################################################################################


####################################output the differential sites#############################################
rawpcut = 0.01             ### cut off for raw pvalue 
adjustpcut = 1          ### cut off for adjusted pvalue
betadiffcut = 0.1         ### cut off for beta value difference
##############################################################################################################

snpsitetestALL = sitetest_median_AL(dataf,gcase=gcase,gcontrol=gcontrol,concov=concov,testmethod = testmethod,Padj=Padj,rawpcut = rawpcut,adjustpcut =adjustpcut,betadiffcut = betadiffcut,paired = paired) ## site-level testing with the "BH" adjustment
#snpsitetestALL = sitetest(dataf,gcase=gcase,gcontrol=gcontrol,concov=concov,testmethod = testmethod,Padj=Padj,rawpcut = rawpcut,adjustpcut =adjustpcut,betadiffcut = betadiffcut,paired = paired) ## site-level testing with the "BH" adjustment
sitetest_out = outputDMfunc(snpsitetestALL,rawpcut=rawpcut,adjustpcut=adjustpcut,betadiffcut=betadiffcut)
sitetest_id<-cbind(Id=row.names(sitetest_out),sitetest_out)
tmp1<-dataf@groupinfo

print("Comenzando test adades")
pvalue_ages <- apply(dataf@bmatrix,1,function(x){summary(lm(unlist(x)~ages))$coefficients["ages","Pr(>|t|)"]})
rsquared_ages <- apply(dataf@bmatrix,1,function(x){summary(lm(unlist(x)~ages))$r.squared})

pvalue_ages_id <- cbind(Id=row.names(dataf@bmatrix),pvalue_ages,rsquared_ages);
print("Finalizando test edades")

pheno_fases <- read.delim(PhenoFileNameIMA,header=TRUE,sep="\t")

#betas_ordered<-cbind(dataf@bmatrix[,tmp1$group=="CONTROL"],dataf@bmatrix[,tmp1$group=="MUESTRA"])
betas_ordered<-cbind(dataf@bmatrix[,pheno_fases$group=="EA_CONTROL"],dataf@bmatrix[,pheno_fases$group=="EA1"],dataf@bmatrix[,pheno_fases$group=="EA2"],dataf@bmatrix[,pheno_fases$group=="EA3"])

betas_ordered_id<-cbind(Id=row.names(betas_ordered),betas_ordered)

sitetest_ages<-merge(sitetest_id,pvalue_ages_id,by.x="Id",by.y="Id")
rm(sitetest_id)
rm(pvalue_ages_id)
#sitetest_merge<-merge(sitetest_id,betas_ordered_id,by.x="Id",by.y="Id")
sitetest_merge<-merge(sitetest_ages,betas_ordered_id,by.x="Id",by.y="Id")
rm(sitetest_ages)
rm(betas_ordered_id)
annot_id<-cbind(Id=row.names(dataf@annot),dataf@annot);
sitetest_merge_annot<-merge(sitetest_merge,annot_id,by.x="Id",by.y="Id")
rm(sitetest_merge)
rm(annot_id)
rm(betas_id_annot)
rm(betas_ordered)
#write.table(sitetest_merge_annot,file="SITE_LEVEL_ANALYSIS_005.txt",sep="\t",row.names=FALSE) ## saving the reults (note that writeXLS won't work on the data exceeds 65535 rows or 256 columns)
sitetest_merge_annot_rounded <- sitetest_merge_annot;
#incluyendo pvalue en el round a 3###
#columnas_paraRound_site <- c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45)
#columnas_paraRound_site <- c(4,5,6,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47)
columnas_paraRound_site <- c(4,5,6,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45)

for(index_round_site in columnas_paraRound_site){  
  sitetest_merge_annot_rounded[,index_round_site] <- round(as.numeric(levels(sitetest_merge_annot[,index_round_site]))[sitetest_merge_annot[, index_round_site]],digits=3)  
}
columnas_paraRound_pvalue <- c(2,3,7,8)
for(index_round in columnas_paraRound_pvalue){  
  sitetest_merge_annot_rounded[,index_round] <- round(as.numeric(levels(sitetest_merge_annot_rounded[,index_round]))[sitetest_merge_annot_rounded[, index_round]],digits=10)  
}


#sitetest_merge_annot_rounded<- sitetest_merge_annot_rounded[,-c(44,46,47,48,49,50,51,52,53,54,57,58,59,60,63,64)]

sitetest_merge_annot_rounded<- sitetest_merge_annot_rounded[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,57,58,67,63,64,47,68,69,70,71,72,73,74,75,76,77,78,79)]

#print("Eliminando cromosomas X e Y");
#sitetest_merge_annot_rounded_XY <- sitetest_merge_annot_rounded[sitetest_merge_annot_rounded$CHR != "X",]
#sitetest_merge_annot_rounded_XY <- sitetest_merge_annot_rounded_XY[sitetest_merge_annot_rounded_XY$CHR != "Y",]
#sitetest_merge_annot_rounded_XY <- sitetest_merge_annot_rounded_XY[sitetest_merge_annot_rounded_XY$CHR != "X/4",]

rm(sitetest_merge_annot)

###Preparamos el listado por zonas para que no contengan los filtrados
# zonas_filtrado_raw<- sitetest_merge_annot_rounded[,c(1,47,48)]
# indice = 0;
# dataf_copia<-dataf@UTR3Ind
# for(pos in dataf@UTR3Ind$SID){
#   indice = indice + 1;
#   mantener = FALSE;
#   for (pos2 in pos){
#     
#     print(pos2)
#     idx <- grep(pos2,zonas_filtrado_raw$Id)
#     if ((length(idx) > 0) & !(mantener)){
#       mantener = TRUE;
#       next;           
#     }
#         
#   }
#   
#   if (!mantener){
#     
#     dataf_copia$SID <- dataf_copia$SID[-indice]
#     
#   }
  
  
}
###
#rm(sitetest_merge_annot_rounded)
#write.table(sitetest_merge_annot_rounded,file="SITE_LEVEL_ANALYSIS_005.txt",row.names=FALSE,sep="\t",dec=",",quote=FALSE, col.names=TRUE) 
#write.table(sitetest_merge_annot_rounded_XY,file="SITE_LEVEL_ANALYSIS_005_ages_REDO.txt",row.names=FALSE,sep="\t",dec=",",quote=FALSE, col.names=TRUE) 
write.table(sitetest_merge_annot_rounded,file="SITE_LEVEL_ANALYSIS_005_ages_noSNP10_noXY_noControls.txt",row.names=FALSE,sep="\t",dec=",",quote=FALSE, col.names=TRUE) 

#rm(sitetest_merge_annot_rounded_XY)
rm(sitetest)
rm(sitetest_merge_annot)
rm(snpsitetestALL)
#######################################################################


########## SITE ANALYSIS ALL##############################
snpsitetestALL = sitetest_median_AL(dataf,gcase=gcase,gcontrol=gcontrol,concov=concov,testmethod = testmethod,Padj=Padj,paired = paired) ## site-level testing with the "BH" adjustment

#snpsitetestALL = sitetest(dataf,gcase=gcase,gcontrol=gcontrol,concov=concov,testmethod = testmethod,Padj=Padj,rawpcut = rawpcut,adjustpcut =adjustpcut,betadiffcut = betadiffcut,paired = paired) ## site-level testing with the "BH" adjustment
#sitetest_sigma <- outfunc_bestSigma_AL(snpsitetestALL,rawpcut=NULL,adjustpcut=NULL,sigma=2)
#sitetest <- sitetest_sigma

sitetest_id<-cbind(Id=row.names(snpsitetestALL),snpsitetestALL)
tmp1<-dataf@groupinfo
#betas_ordered<-cbind(dataf@bmatrix[,tmp1$group=="CONTROL"],dataf@bmatrix[,tmp1$group=="MUESTRA"])
betas_ordered<-cbind(dataf@bmatrix[,pheno_fases$group=="EA_CONTROL"],dataf@bmatrix[,pheno_fases$group=="EA1"],dataf@bmatrix[,pheno_fases$group=="EA2"],dataf@bmatrix[,pheno_fases$group=="EA3"])

betas_ordered_id<-cbind(Id=row.names(betas_ordered),betas_ordered)
sitetest_merge<-merge(sitetest_id,betas_ordered_id,by.x="Id",by.y="Id")
annot_id<-cbind(Id=row.names(dataf@annot),dataf@annot[,c(22,12,13,24,26)]);
rm(betas_ordered)
rm(betas_ordered_id)
rm(all_snps)
rm(beta)
rm(betas_converted)
rm(betas_id_annot)
#rm(bmiq_betas_samesamples)
rm(snpsitetestALL)
rm(sitetest_id)
sitetest_merge_annot<-merge(sitetest_merge,annot_id,by.x="Id",by.y="Id")

rm(sitetest_merge)
rm(annot_id)

sitetest_merge_annot_rounded <- sitetest_merge_annot;

rm(sitetest_merge_annot)

rm(data_methylumi)
rm(data_methylumi_backup)
rm(data_methylumi_filtered)
rm(dataf)
rm(dataf2)
rm(dataf_backup)
columnas_paraRound <- c(4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43)
for(index_round in columnas_paraRound){  
  sitetest_merge_annot_rounded[,index_round] <- round(as.numeric(levels(sitetest_merge_annot_rounded[,index_round]))[sitetest_merge_annot_rounded[, index_round]],digits=3)  
}
columnas_paraRound_pvalue <- c(2,3)
for(index_round in columnas_paraRound_pvalue){  
  sitetest_merge_annot_rounded[,index_round] <- round(as.numeric(levels(sitetest_merge_annot_rounded[,index_round]))[sitetest_merge_annot_rounded[, index_round]],digits=10)  
}

# print("Eliminando cromosomas X e Y");
# sitetest_merge_annot_rounded_XY <- sitetest_merge_annot_rounded[sitetest_merge_annot_rounded$CHR != "X",]
# sitetest_merge_annot_rounded_XY <- sitetest_merge_annot_rounded_XY[sitetest_merge_annot_rounded_XY$CHR != "Y",]
# sitetest_merge_annot_rounded_XY <- sitetest_merge_annot_rounded_XY[sitetest_merge_annot_rounded_XY$CHR != "X/4",]

#write.table(cbind(Id=row.names(TSS1500sig),TSS1500sig),file="TSS1500_annot.txt",row.names=FALSE) ## saving the reults (note that writeXLS won't work on the data exceeds 65535 rows or 256 columns)
write.table("[Parametros]",file="SITE_LEVEL_ANALYSIS_400000_noSNO10_noXY_noControls.txt", row.names=FALSE,sep="\t",dec=",",quote=FALSE,col.names=FALSE)
write.table("All sites",file="SITE_LEVEL_ANALYSIS_400000_noSNO10_noXY_noControls.txt", row.names=FALSE,sep="\t",dec=",",quote=FALSE,col.names=FALSE,append=TRUE)
write.table("",file="SITE_LEVEL_ANALYSIS_400000_noSNO10_noXY_noControls.txt", row.names=FALSE,sep="\t",dec=",",quote=FALSE,col.names=FALSE, append=TRUE)

#write.table(sitetest_merge_annot_rounded_XY,file="SITE_LEVEL_ANALYSIS_400000_noSNO10_noXY_noControls.txt", row.names=FALSE,sep="\t",dec=",",quote=FALSE,col.names=TRUE, append=TRUE) ## saving the reults (note that writeXLS won't work on the data exceeds 65535 rows or 256 columns)
write.table(sitetest_merge_annot_rounded,file="SITE_LEVEL_ANALYSIS_400000_noSNO10_noXY_noControls.txt", row.names=FALSE,sep="\t",dec=",",quote=FALSE,col.names=TRUE, append=TRUE) ## saving the reults (note that writeXLS won't work on the data exceeds 65535 rows or 256 columns)

selected_site <- sitetest_merge_annot_rounded[sitetest_merge_annot_rounded$"Adjust Pval" <= 0.05,-c(2,3,4,5,6,44,45,46)]
mdata_indiv <- melt(selected_site, id="Id");
colnames(mdata_indiv) <- c("Id","Sample","Beta")

#####################################


####################################differential methylation regions PARAMETERS#############################################
indexmethodRegion = "median"
testmethodRegion = "wilcox"


####################################REGION ANALYSI#############################################
tmp1<-dataf@groupinfo
betas_ordered<-cbind(dataf@bmatrix[,tmp1$group=="CONTROL"],dataf@bmatrix[,tmp1$group=="MUESTRA"])
colnames(betas_ordered) <- paste(colnames(betas_ordered),"_Indiv",sep="")
betas_ordered_id<-cbind(Id=row.names(betas_ordered),betas_ordered)

slots<-slotNames(dataf);
regiones <- slots[-c(1,2,3,4)];

regiones <- regiones[11]
for (region_name in regiones){
print(region_name);
region_analisis_filename <- paste("PRUEBA_REGION ANALYSIS_",region_name,".txt",sep="");
region_analisis_all_filename <- paste("PRUEBA_REGION ANALYSIS_ALL_",region_name,".txt",sep="");
print("Cargando betas");
betaTSS1500 <- indexregionfunc_AL(indexlist=slot(dataf,region_name),beta=dataf@bmatrix,indexmethod=indexmethodRegion)
print("Cargando porcentajes");
porcentajeregion_TSS1500 <- indexregionfunc_porcentajeregion_after_AL(indexlist=slot(dataf,region_name),beta=dataf@bmatrix,indexmethod=indexmethodRegion)
mean_porcentajeregion_TSS1500 <- apply(porcentajeregion_TSS1500,1,function(x){return(mean(x[x != 0]))})
porcentajeregion_TSS1500_id <- cbind(Id=row.names(porcentajeregion_TSS1500),mean_porcentajeregion_TSS1500)
colnames(betaTSS1500) <- paste(colnames(betaTSS1500),"_Region",sep="")
betaTSS1500_id <- cbind(Id=row.names(betaTSS1500),betaTSS1500[,tmp1$group=="CONTROL"],betaTSS1500[,tmp1$group=="MUESTRA"])
print("Comenzando el test")
TSS1500testALL <- testfunc_wilcox_median_AL(eset = betaTSS1500,testmethod=testmethodRegion,Padj="BH",concov="OFF",groupinfo = dataf@groupinfo,gcase ="MUESTRA",gcontrol="CONTROL",paired = FALSE,porcentajeNAenGrupo)
print("Finalizado el test")
TSS1500test <- outputDMfunc(TSS1500testALL,rawpcut=rawpcut,adjustpcut=adjustpcut,betadiffcut=betadiffcut)
TSS1500test_id<-cbind(Id=row.names(TSS1500test),TSS1500test)
TSS1500_tamanyoregion_porcentajereunion <- merge(TSS1500test_id,porcentajeregion_TSS1500_id,by.x="Id",by.y="Id",all.x=TRUE)
TSS1500testALL_id<-cbind(Id=row.names(TSS1500testALL),TSS1500testALL)
write.table(TSS1500testALL_id,file=region_analisis_all_filename,row.names=FALSE,sep="\t",dec=",",quote=FALSE) 
porcentajeTSS1500selected <- 100*dim(TSS1500test)[1]/dim(TSS1500testALL)[1];
listtoannot <- rownames(TSS1500test)
fullannotInd <- fullannot
fullIndexannot <- get(region_name)
filteredannot <- dataf@annot
filteredIndexannot <- slot(dataf,region_name)
TSS1500sig <- annotfunc(listtoannot,fullannot,filteredannot,fullIndexannot,filteredIndexannot,category = "region")
TSS1500sig_id <- cbind(Id=row.names(TSS1500sig),TSS1500sig)
TSS1500sig_id_cgs <-TSS1500sig_id[,4]

print("Montando estructura de probes")
TSS1500_cgs<-data.frame();
index_tss1500<-1;
index<-1;
cgs<-matrix();
for(i in TSS1500sig_id_cgs){
  listado_cgs<-strsplit(i, "/");
  for(j in listado_cgs){    
    TSS1500_cgs =rbind(TSS1500_cgs,data.frame(Probe=j,Id=TSS1500sig_id[index_tss1500,1]))
  };
  index_tss1500 = index_tss1500 + 1
}
TSS1500_pruebas<-merge(TSS1500_cgs,TSS1500sig_id,by.x="Id",by.y="Id",all.x=TRUE)

print("Merge de todas las estructuras de datos")
TSS1500_annot_snp<-merge(TSS1500_pruebas[,c(1,16,17,26,22,23,2,3,5,6,27,28,29,30,31,32,33,34,35,36,37,38)],all_snps[,c(1,10,15,16,17,18,19)],by.x="Probe",by.y="PROBE",all.x=TRUE)
TSS1500_annot_snp_means<-merge(TSS1500_tamanyoregion_porcentajereunion,TSS1500_annot_snp,by.x="Id",by.y="Id",all.x=TRUE)
TSS1500test_merge_betaTSS1500_id <-merge(TSS1500_annot_snp_means,betaTSS1500_id,by.x="Id",by.y="Id")
TSS1500_merge_betas<-merge(TSS1500test_merge_betaTSS1500_id,betas_ordered_id,by.x="Probe",by.y="Id",all.x=TRUE)
TSS1500_merge_betas_rounded <- TSS1500_merge_betas;
#columnas_paraRound <- c(3,4,5,6,7,8,9,10,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110)
columnas_paraRound_3decimales <- c(5,6,7,8,9,10,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110)
columnas_paraRound_7decimales <- c(3,4)
for(index_round in columnas_paraRound_3decimales){  
  TSS1500_merge_betas_rounded[,index_round] <- round(as.numeric(levels(TSS1500_merge_betas[,index_round]))[TSS1500_merge_betas[, index_round]],digits=3)  
}
for(index_round in columnas_paraRound_7decimales){  
  TSS1500_merge_betas_rounded[,index_round] <- round(as.numeric(levels(TSS1500_merge_betas[,index_round]))[TSS1500_merge_betas[, index_round]],digits=7)  
}

# print("Eliminando cromosomas X e Y");
# TSS1500_merge_betas_rounded_XY <- TSS1500_merge_betas_rounded[TSS1500_merge_betas_rounded$CHR != "X",]
# TSS1500_merge_betas_rounded_XY <- TSS1500_merge_betas_rounded_XY[TSS1500_merge_betas_rounded_XY$CHR != "Y",]
# TSS1500_merge_betas_rounded_XY <- TSS1500_merge_betas_rounded_XY[TSS1500_merge_betas_rounded_XY$CHR != "X/4",]

print("Escribiendo resultados")
write.table("[Parametros]",file=region_analisis_filename, row.names=FALSE,sep="\t",dec=",",quote=FALSE,col.names=FALSE)
write.table(paste("Porcentaje Seleccionados",porcentajeTSS1500selected,sep="\t"),file=region_analisis_filename,row.names=FALSE,sep="\t",dec=",",quote=FALSE,append=TRUE, col.names=FALSE)
write.table(paste("Estadistico cada region",indexmethodRegion,sep="\t"),file=region_analisis_filename,row.names=FALSE,sep="\t",dec=",",quote=FALSE,append=TRUE, col.names=FALSE)
write.table(paste("Estadistico grupo diferencial","median",sep="\t"),file=region_analisis_filename,row.names=FALSE,sep="\t",dec=",",quote=FALSE,append=TRUE, col.names=FALSE)
write.table(paste("Metodo diferenciacion",testmethodRegion,sep="\t"),file=region_analisis_filename,row.names=FALSE,sep="\t",dec=",",quote=FALSE,append=TRUE, col.names=FALSE)
write.table(paste("NA's","Incluidos",sep="\t"),file=region_analisis_filename,row.names=FALSE,sep="\t",dec=",",quote=FALSE,append=TRUE, col.names=FALSE)
write.table(paste("% NA's admitidos region","Todos indicando porcentaje medio de NA's",sep="\t"),file=region_analisis_filename,row.names=FALSE,sep="\t",dec=",",quote=FALSE,append=TRUE, col.names=FALSE)
write.table(paste("% NA's admitidos grupo","Incluida columna para filtrado posterior",sep="\t"),file=region_analisis_filename,row.names=FALSE,sep="\t",dec=",",quote=FALSE,append=TRUE, col.names=FALSE)
write.table("",file=region_analisis_filename,row.names=FALSE,sep="\t",dec=",",quote=FALSE,append=TRUE,col.names=FALSE)
write.table(TSS1500_merge_betas_rounded_XY,file=region_analisis_filename,row.names=FALSE,sep="\t",dec=",",quote=FALSE,append=TRUE, col.names=TRUE) 


print("Listado y graficas seleccionados p-value adjusted")

seleccionados <- NULL;

seleccionados <- TSS1500_merge_betas_rounded_XY[TSS1500_merge_betas_rounded_XY[,4] <=0.05,]
if (dim(seleccionados)[1] != 0){
selection_indiv <- seleccionados[,c(2,seq(74,110))]
selection_group <- seleccionados[,c(2,seq(37,73))]
selection_region <- seleccionados[,c(2,4,5,6,7)]

mdata_indiv <- melt(selection_indiv, id="Id");
colnames(mdata_indiv) <- c("Id","Sample","Beta")
mdata_group <- melt(selection_group, id="Id")
colnames(mdata_group) <- c("Id","Sample","Beta")

### Los ids deben ser identicos para grupo, indiv y region
ids <- unique(sort(selection_indiv$Id))

graficas_analisis_filename <- paste("GRAFICAS_PVALUEADJUSTED_",region_name,".pdf",sep="");
pdf(graficas_analisis_filename)

for (id in ids){
  
  selection_indiv_i <- mdata_indiv[mdata_indiv$Id == id,];
  selection_group_i <- mdata_group[mdata_group$Id == id,];
  selection_region_i <- selection_region[selection_region$Id == id,];
  #boxplot(Beta ~ Sample, data=selection_i,col=c(rep("green",13),rep("red",24)),ylim=c(0,1))

  print(ggplot(selection_indiv_i, aes(factor(Sample),Beta))+ geom_boxplot(fill=c(rep("blue",12),rep("red",25)))+ geom_hline(yintercept=mean(selection_region_i$Median_MUESTRA),colour="red", size=3) + geom_hline(yintercept=mean(selection_region_i$Median_CONTROL),colour="blue", size =3) + ggtitle(paste(id,mean(selection_region_i[,3]),mean(selection_region_i[,2]))) + geom_point() + ylim(0,1))
  
}

dev.off()
}
print("Borrando variables intermedias")
rm(TSS1500_pruebas)
rm(TSS1500_tamanyoregion_porcentajereunion)
rm(TSS1500testALL_id)
rm(TSS1500test_id)
rm(region_analisis_filename)
rm(betaTSS1500_id)
rm(mean_porcentajeregion_TSS1500)
rm(mean_porcentajeregion_TSS1500)
rm(betaTSS1500)
rm(porcentajeregion_TSS1500)
rm(TSS1500testALL)
rm(TSS1500test)
rm(TSS1500sig)
rm(TSS1500sig_id)
rm(TSS1500sig_id_cgs)
rm(TSS1500_cgs)
rm(TSS1500_annot_snp)
rm(TSS1500_annot_snp_means)
rm(TSS1500_merge_betas)
rm(TSS1500test_merge_betaTSS1500_id)
 rm(TSS1500_merge_betas_rounded)
 rm(TSS1500_merge_betas_rounded_XY)
}
###########################################################################






