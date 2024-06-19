
doSVA <- function(meth450_data) {
  
  pheno = pData(meth450_data)
  edata = getM(meth450_data)
  
  mod = model.matrix(~as.factor(c("ALZHEIMER","CONTROL")[as.factor(estadio=="EA_CONTROL")]) + as.numeric(edad), data=pheno)
  
  mod0 = model.matrix(~as.numeric(edad),data=pheno)
  
  n.sv = num.sv(edata,mod,method="leek")
  
  svobj = sva(edata,mod,mod0)

  pValues = f.pvalue(edata,mod,mod0)
  qValues = p.adjust(pValues,method="BH")
  
}

# To return DMPS from different values of k:
runSva <- function(rgSet, pheno){
  
  library(sva)
  library(minfi)
  
  
  if (ncol(rgSet)!=length(pheno)){
    stop("Pheno does not have the right length")
  }
  
  n <- length(pheno)
  sampleNames <- colnames(rgSet)
  if (sum(names(pheno) %in% sampleNames)!=n){
    stop("pheno does not match sample names")
  }
  
  # To make sure pheno as the right order
  sampleNames <- colnames(rgSet)
  pheno <- pheno[match(sampleNames, names(pheno))]
  
  cat("[SVA450k] Extraction of the m-values \n")
  raw <- preprocessRaw(rgSet)
  mvalues <- getM(raw)
  
  
  # Removing missing values:
  mvalues <- mvalues[complete.cases(mvalues),]
  
  
  # Construction of the phenotype matrix
  cat("[SVA450k] Constructing the pheno matrix \n")
  pheno <- as.factor(as.numeric(as.factor(pheno))) # Make sure all levels have values
  phenoMatrix <- model.matrix(~as.factor(pheno))
  
  # Is that correct for all phenotypes ?
  mod <- phenoMatrix
  mod0 <- matrix(phenoMatrix[,1],ncol=1)
  
  
  matrix <- mvalues
  
  # To remove rows with infinite values:
  inf.values <- rowSums(is.infinite(matrix))
  inf.row <- which(inf.values!=0)
  matrix <- matrix[-inf.row,]
  if (length(inf.row)!=0){
    matrix <- matrix[-inf.row,]
  }
  
  
  
  cat("[SVA450K] Running SVA")
  sva.object <- sva(dat = matrix, mod=mod, mod0 = mod0, method="irw")
  
  
  mod.sv <- cbind(mod,sva.object$sv)
  mod0.sv <- cbind(mod0,sva.object$sv)
  cat("[SVA450k] Computing F-statistics")
  results.sva <- sva:::fstats(matrix,mod.sv,mod0.sv)
  results.sva <- results.sva[order(results.sva, decreasing=T),]
  return(list(sva.object=sva.object, results.sva=results.sva))		
  
}


cor_M_edad <- apply(M,1,cor, y =samps$edad)
hist(cor_M_edad)
amieloide <- read.csv('data/alzheimer/AMILOIDE_081014.csv')
head(M[,amieloide$SampleName])
cor_M_amieloide <- apply(M[,amieloide$SampleName],1,cor, y =amieloide$Mean.total.area)

hist(cor_M_amieloide)
hist(cor_M_edad)
edad <- samps$edad
pvalue_ages <- apply(M,1,function(x){summary(lm(unlist(x)~edad))$coefficients["edad","Pr(>|t|)"]})
hist(pvalue_ages)
qvalue_ages <- p.adjust(pvalue_ages, method='fdr')
b_edad <- b[names(qvalue_ages[qvalue_ages<0.05]),]
names(edad) <- rownames(samps)

edad_sorted <- sort(edad)

sdata <- b_edad[order(qvalue_ages[qvalue_ages<0.00001]), names(edad_sorted)]
sdata <- melt(sdata)


profile_plot <- ggplot(sdata, aes(x=X2, y=value, group=X1)) + geom_line()
geom_point(aes, size = 2.5) +
geom_line(aes, size = 1)

install_github("ropensci/plotly")
require(plotly)
py <- plotly(username="alabarga", key="a935mlupe5")
py$ggplotly(profile_plot)

install_github('ramnathv/rCharts')
require(rCharts)
rPlot(value ~ X2, color='X1', data=sdata, type='line')



