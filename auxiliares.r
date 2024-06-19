require(RCurl)

disgenet_assoc<-function(entrezId='7157'){
url <- "http://www.disgenet.org/oql"
oqlTemplate <- "DEFINE
c0='/data/gene_disease_score_onexus',
c1='/data/diseases',
c2='/data/genes',
c3='/data/sources'
ON
'http://bitbucket.org/janis_pi/disgenet_onexus.git'
SELECT
c1 (cui, name, cui, name, diseaseClassName, STY, cui, name),
c2 (geneId, name, geneId, name, uniprotId, description, pathName, pantherName, geneId, name),
c0 (score, diseaseId, score, geneId, score, diseaseId, geneId, score, diseaseId, geneId, score, pmids)
FROM
c0
WHERE
(
c2 = '%ENTREZID%'
AND
c3 = 'ALL'
)
ORDER BY
c0.score DESC" 


oql <- gsub(oqlTemplate,'%ENTREZID%', entrezId)
dataTsv <- getURLContent(url, readfunction =charToRaw(oql), upload = TRUE, customrequest = "POST")
disease_assoc <- read.csv(textConnection(dataTsv), header = TRUE, sep="\t")

disease_assoc
}

####################################generate annotation#############################################
# 
# data_nonorm =IMA.methy450R(fileName = MethyFileName_NoNorm,columnGrepPattern=list(beta=".AVG_Beta",detectp=".Detection.Pval"),groupfile = PhenoFileName)
# file.rename("./QC.pdf","./QC_nonorm.pdf")
#   
# dataf = IMA.methy450PP(data,na.omit = na.omit,normalization=normalization,peakcorrection = peakcorrection,transfm = transfm,samplefilterdetectP = samplefilterdetectP,samplefilterperc = samplefilterperc,sitefilterdetectP = sitefilterdetectP,locidiff = locidiff, locidiffgroup = locidiffgroup,XYchrom = XYchrom,snpfilter = snpfilter) ## QC filtering
# fullannot = dataf2@annot 
# temp = c("TSS1500Ind","TSS200Ind","UTR5Ind", "EXON1Ind","GENEBODYInd","UTR3Ind","ISLANDInd","NSHOREInd","SSHOREInd","NSHELFInd", "SSHELFInd") 
# for( i in 1:11){eval(parse(text=paste(temp[i],"=dataf2@",temp[i],sep="")))} 
# eval(parse(text = paste("save(fullannot", paste(temp,collapse = ","), "file = 'fullannotInd_all.rda')", sep = "," ))) 

##############################################################################################################
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}

normalize450 <- function(meth450_data, method="BMIQ") {
  if (method == "BMIQ") {
    meth450_norm <- BMIQ(meth450_data)
  } else if (method == "SWAN") {
    meth450_norm <- SWAN(meth450_data)
  } else if (method == "dasen") {
    meth450_norm <- dasen(meth450_data)
  }
  meth450_norm
}

topTest <- function(res, diff=0.1, pval= 0.1,  pval.var= "adj.P.Val", diff.var="Beta.Difference"){
  
  res[which((abs(res[,diff.var])>=diff) & (res[,pval.var]<=pval)),]

}

topTest2 <- function(data.test, logFC=0.1, adjust.p = 0.1){
  
  data.test[which(abs(data.test[,"logFC"])>logFC & data.test[,"adj.P.Val"]<adjust.p),]
  
}

add.betaDiff<- function(res, b, clases, l1=NULL, l2=NULL){
 
   bd <- betaDiff(b, clases, l1=l1, l2=l2)
   for (l in levels(clases)){
     name <- paste('AVG',l)
     res[name] <- rowMeans(b[rownames(res),clases == l], na.rm=TRUE)
   }
   res$Beta.Difference <- bd[rownames(res)]
   res <- cbind(res,b[rownames(res),])
   res
}

betaDiff <- function(b, clases, l1=NULL, l2=NULL){
  l = levels(clases)
  if (is.null(l1)) l1 = l[1]
  if (is.null(l2)) l2 = l[2]
  
  difb = rowMeans(b[,clases == l1], na.rm=TRUE) - rowMeans(b[,clases == l2], na.rm=TRUE)
  difb
}  


grep_count <- function(x,meth450_data) {
  length(grep(x,meth450_data@featureData@data$UCSC_REFGENE_GROUP))
}

count_regions_1 <-function(meth450_data){
  regiones <- c('TSS1500','TSS200','Body','1stExon',"3'UTR","5'UTR")
  
  cuenta <- apply(as.array(regiones),MARGIN=1,FUN=grep_count,meth450_data=meth450_data)
  names(cuenta) <- regiones
  cuenta
}

count_regions_2 <-function(meth450_data){
  cuenta <- table(meth450_data@featureData@data$RELATION_TO_UCSC_CPG_ISLAND)
  cuenta[2:6]
}

count_regions <-function(meth450_data){
  
  c(count_regions_1(meth450_data),count_regions_2(meth450_data))
}

region_analysis <- function(meth450_data, meth450_data_selected){
  
  e <- count_regions(meth450_data) / nrow(meth450_data)
  o <- count_regions(meth450_data_selected) / nrow(meth450_data_selected)
  
  res <- log2(o/e)
  dtm <- data.frame(region=names(res),value=res)
  dtm[order(dtm$value),]
}

library(ggplot2)

region_plot<- function(dtm){
    
  dtm$colour <- ifelse(dtm$value < 0, "firebrick1",  "steelblue")
  
  dtm$hjust <- ifelse(dtm$value > 0, 1.3, -0.3)
  
  ggplot(dtm, aes(region, value, label = region, hjust = hjust)) +
    geom_text(aes(y = 0, colour = colour)) +
    geom_bar(stat = "identity", aes(fill = colour))
  
  last_plot() +
   coord_flip() + 
   labs(x = "", y = "") + 
   scale_x_discrete(breaks = NULL) +
   theme_bw() +
   theme(legend.position = "none")
  
}

removeNA <- function(x){
  x[apply(x, 1, function(y) !all(is.na(y))),]
} 

predictSex <-function (meth450_data_original,sample, thr=-2) {
  x <- assayData(meth450_data_original[,sample])$Intensity[meth450_data_original@featureData@data$CHR=="X"]
  
  y <- assayData(meth450_data_original[,sample])$Intensity[meth450_data_original@featureData@data$CHR=="Y"]
  
  delta <- log2(median(y)) - log2(median(x))
  
  if (delta < thr) {
    pred = "F"
  } else {
    pred = "M"
    
  }
  
  pred
  
}

convertBMIQ2IMA <- function(bmiq_betas,ima_betas){
  
  bmiq_cgs_names <- rownames(bmiq_betas);
  ima_cgs_names <- rownames(ima_betas);
  
  bmiq_samples_names <- colnames(bmiq_betas);
  ima_samples_names <- colnames(ima_betas);
  
  out <- ima_betas;
  out <- NA;
  
  if (length(bmiq_cgs_names) == length(ima_cgs_names) && sum(1*(bmiq_cgs_names == ima_cgs_names)) == length(bmiq_cgs_names)){
    print("Listado de cgs iguales; se puede continuar la conversion");
    

    if (length(bmiq_samples_names) == length(ima_samples_names) && sum(1*(bmiq_samples_names == ima_samples_names)) == length(bmiq_samples_names)){
      
      print("Orden de las muestras igual; se puede continar la conversion");
      
      out <- bmiq_betas;
    }else{
      
      print("Alguna muestra no coincide")
      print("BMIQ")
      print(bmiq_samples_names)
      print("IMA")
      print(ima_samples_names)
      
    }
  }
  
  return(out);
  
}

sitetest_local <-function (dataf=NULL, beta=NULL, 
                          grouplev= NULL, gcase = NULL, gcontrol = NULL, 
                          testmethod = c("wilcox", "limma", "pooled", "satterthwaite"), 
                          Padj = "BH", concov = "OFF", 
                          rawpcut = NULL, adjustpcut = NULL, betadiffcut = NULL, paired = FALSE) 
{
  
  if (!is.null(dataf)){
   beta = dataf@bmatrix
   group = dataf@groupinfo
   grouplev = group[, 2]
  }
  
  if (concov == "ON") {
    cat("Performing linear regression....\n")
    require(MASS)
    testout = apply(beta, 1, function(x) {
      temp = summary(lm(x ~ as.numeric(as.character(grouplev))))
      pvalue = temp$coefficients[2, c(1, 4)]
      return(pvalue)
    })
    adjustP = p.adjust(testout[2, ], method = Padj)
    out = cbind(testout[2, ], adjustP, testout[1, ])
    rownames(out) = rownames(beta)
    colnames(out) = c("P.Value", "adj.P.Value", "Coefficient")
  }
  else {
    caseind = which(grouplev %in% gcase)
    controlind = which(grouplev %in% gcontrol)
    if (paired == TRUE) {
      lev1 = caseind[order(group[caseind, 3])]
      lev2 = controlind[order(group[controlind, 3])]
    }
    else {
      lev1 = caseind
      lev2 = controlind
    }
    eset = beta[, c(lev1, lev2)]
    
    if (testmethod == "wilcox") {
      cat("Performing Wilcox testing ...\n")
      testout = apply(eset, 1, function(x) {
        wilcox.test(x[1:length(lev1)], x[(length(lev1) + 
                                            1):(length(lev1) + length(lev2))], paired = paired, na.action="na.exclude")$p.value
      })
    }
    if (testmethod == "limma") {
      require(limma)
      cat("Performing limma...\n")
      TS = as.factor(c(rep("T", length(lev1)), rep("C", 
                                                   length(lev2))))
      SS = as.factor(rep(1:length(lev1), 2))
      if (paired == FALSE) {
        design = model.matrix(~0 + TS)
        rownames(design) = colnames(eset)
        colnames(design) = c("C", "T")
        fit = lmFit(eset, design)
        cont.matrix = makeContrasts(comp = T - C, levels = design)
        fit2 = contrasts.fit(fit, cont.matrix)
        fit2 = eBayes(fit2)
        result1 = topTable(fit2, coef = 1, adjust.method = Padj, 
                           number = nrow(fit2))
      }
      else {
        design = model.matrix(~SS + TS)
        cat("Here is your design matrix\n")
        print(design)
        fit = lmFit(eset, design)
        fit2 = eBayes(fit)
        result1 = topTable(fit2, coef = "TST", adjust.method = Padj, 
                           number = nrow(fit2))
      }
      testout = result1[match(rownames(eset), result1[, 
                                                      1]), "P.Value"]
    }
    if (testmethod == "pooled") {
      cat("Performing pooled t.test...\n")
      testout = apply(eset, 1, function(x) {
        t.test(x[1:length(lev1)], x[(length(lev1) + 1):(length(lev1) + 
                                                          length(lev2))], var.equal = TRUE, paired = paired)$p.value
      })
    }
    if (testmethod == "satterthwaite") {
      cat("Performing satterthwaite t.test...\n")
      testout = apply(eset, 1, function(x) {
        t.test(x[1:length(lev1)], x[(length(lev1) + 1):(length(lev1) + 
                                                          length(lev2))], paired = paired)$p.value
      })
    }
    adjustP = p.adjust(testout, method = Padj)
    difb = apply(eset, 1, function(x) {
      median(x[1:length(lev1)],na.rm=TRUE) - median(x[(length(lev1) + 
                                              1):ncol(eset)],na.rm=TRUE)
    })
    out = cbind(testout, adjustP, difb, rowMedians(eset[, 1:length(lev1)],na.rm=TRUE), 
                rowMedians(eset[, (length(lev1) + 1):ncol(eset)],na.rm=TRUE))
    
    rownames(out) = rownames(eset)
    colnames(out) = c("P.Value", "adj.P.val", "Beta.Difference", 
                      paste("Median", paste(gcase, collapse = "_"), sep = "_"), 
                      paste("Median", paste(gcontrol, collapse = "_"), sep = "_"))
  }
  
  out <- out[order(out[,c("adj.P.val")]),]

  return(out)
}

testfunc_wilcox_mean_AL <- function (eset, concov = c("ON", "OFF"), testmethod = c("wilcox", "limma", "pooled", "satterthwaite", "paired"), Padj = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), groupinfo, gcase = "g1", gcontrol = "g2", paired = FALSE, porcentaje) 
{
    grouplev = groupinfo[, 2]
    caseind = which(grouplev %in% gcase)
    controlind = which(grouplev %in% gcontrol)
    if (paired == TRUE) {
      lev1 = caseind[order(groupinfo[caseind, 3])]
      lev2 = controlind[order(groupinfo[controlind, 3])]
    }
    else {
      lev1 = caseind
      lev2 = controlind
    }
    eset = eset[, c(lev1, lev2)]
    if (testmethod == "wilcox") {
      cat("Performing Wilcoxon test...\n")
      testout = apply(eset, 1, function(x) {
        wilcox.test(x[1:length(lev1)], x[(length(lev1) + 
                                            1):ncol(eset)], paired = paired)$p.value
      })
    }
    adjustP = p.adjust(testout, method = Padj)
    difb = apply(eset, 1, function(x) {
      mean(x[1:length(lev1)],na.rm=TRUE) - mean(x[(length(lev1) + 
                                          1):ncol(eset)],na.rm=TRUE)
    })
    
    betaMuestra = rowMeans(eset[, 1:length(lev1)],na.rm=TRUE);
    betaControl = rowMeans(eset[, (length(lev1) + 1):ncol(eset)],na.rm=TRUE)
    porcentajeMuestra = apply(eset[, 1:length(lev1)],1,function(x){100*sum(1*is.na(x))/(length(x))});
    porcentajeControl = apply(eset[, (length(lev1) + 1):ncol(eset)],1,function(x){100*sum(1*is.na(x))/(length(x))});
    out = cbind(testout, adjustP, difb, betaMuestra, betaControl,porcentajeMuestra,porcentajeControl)
    rownames(out) = rownames(eset)
    colnames(out) = c("P-Value", "Adjust Pval", "beta-Difference", 
                      paste("Mean", paste(gcase, collapse = "_"), sep = "_"), 
                      paste("Mean", paste(gcontrol, collapse = "_"), sep = "_"),
                      paste("Porcentaje", paste(gcase, collapse = "_"), sep = "_"),
                      paste("Porcentaje", paste(gcontrol, collapse = "_"), sep = "_"))
  
  return(out)
}

testfunc_wilcox_median_AL <- function (eset, concov = c("ON", "OFF"), testmethod = c("wilcox", "limma", "pooled", "satterthwaite", "paired"), Padj = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), groupinfo, gcase = "g1", gcontrol = "g2", paired = FALSE, porcentaje) 
{
  grouplev = groupinfo[, 2]
  caseind = which(grouplev %in% gcase)
  controlind = which(grouplev %in% gcontrol)
  if (paired == TRUE) {
    lev1 = caseind[order(groupinfo[caseind, 3])]
    lev2 = controlind[order(groupinfo[controlind, 3])]
  }
  else {
    lev1 = caseind
    lev2 = controlind
  }
  eset = eset[, c(lev1, lev2)]
  if (testmethod == "wilcox") {
    cat("Performing Wilcoxon test...\n")
    testout = apply(eset, 1, function(x) {
      wilcox.test(x[1:length(lev1)], x[(length(lev1) + 
                                          1):ncol(eset)], paired = paired)$p.value
    })
  }
  adjustP = p.adjust(testout, method = Padj)
  difb = apply(eset, 1, function(x) {
    median(x[1:length(lev1)],na.rm=TRUE) - median(x[(length(lev1) + 
                                                   1):ncol(eset)],na.rm=TRUE)
  })
  
  betaMuestra = rowMedians(eset[, 1:length(lev1)],na.rm=TRUE);
  betaControl = rowMedians(eset[, (length(lev1) + 1):ncol(eset)],na.rm=TRUE)
  porcentajeMuestra = apply(eset[, 1:length(lev1)],1,function(x){100*sum(1*is.na(x))/(length(x))});
  porcentajeControl = apply(eset[, (length(lev1) + 1):ncol(eset)],1,function(x){100*sum(1*is.na(x))/(length(x))});
  out = cbind(testout, adjustP, difb, betaMuestra, betaControl,porcentajeMuestra,porcentajeControl)
  rownames(out) = rownames(eset)
  colnames(out) = c("P-Value", "Adjust Pval", "beta-Difference", 
                    paste("Median", paste(gcase, collapse = "_"), sep = "_"), 
                    paste("Median", paste(gcontrol, collapse = "_"), sep = "_"),
                    paste("Porcentaje", paste(gcase, collapse = "_"), sep = "_"),
                    paste("Porcentaje", paste(gcontrol, collapse = "_"), sep = "_"))
  
  return(out)
}

indexregionfunc_AL <- function (indexlist, beta, indexmethod = c("mean", "median", "tbrm"), porcentaje = 200) 
{
  nr = length(indexlist$PID)
  temp2 = matrix(NA, nrow = nr, ncol = ncol(beta))
  temp3 = matrix(NA, nrow = nr, ncol = ncol(beta))
  rownames(temp2) = names(indexlist$SID)
  colnames(temp2) = colnames(beta)
  for (i in 1:nr) {  

    if (i%%(round(nr/10)) == 0){
      print ("10% más...");      
    }
    temp = beta[indexlist$PID[[i]], ]
    if (length(indexlist$PID[[i]]) == 1) {
      temp2[i, ] = temp;
      temp3[i, ] = 100*1*is.na(temp);
      temp2[i,temp3[i,] >= porcentaje] = NA;
    }
    else {
      if (indexmethod == "tbrm") {
        temp3[i, ] = apply(temp,2, function(x){100*sum(1*is.na(x))/(length(x))})
        temp2[i, ] = apply(temp, 2, eval(indexmethod))
        temp2[temp3[i,] >= porcentaje] = NA;
      }
      else {
        temp3[i, ] = apply(temp,2, function(x){100*sum(1*is.na(x))/(length(x))})   
        temp2[i, ] = apply(temp, 2, eval(indexmethod), na.rm = TRUE)
        temp2[i,temp3[i,] >= porcentaje] = NA;       
        
      }
    }
  }
  return(temp2)
}

indexregionfunc_tamanyoregion_before_AL <- function (indexlist, beta, indexmethod = c("mean", "median", "tbrm")) 
{
  nr = length(indexlist$PID)
  temp2 = matrix(NA, nrow = nr, ncol = ncol(beta))
  rownames(temp2) = names(indexlist$SID)
  colnames(temp2) = colnames(beta)
  for (i in 1:nr) {  
    if (i%%(round(nr/10)) == 0){
      print ("10% más...");      
    }
    temp = beta[indexlist$PID[[i]], ]
    print(length(row.names(temp)));
    temp2[i, ] = length(row.names(temp));

  }
  return(temp2)
}


indexregionfunc_porcentajeregion_after_AL <- function (indexlist, beta, indexmethod = c("mean", "median", "tbrm"), porcentaje) 
{
  nr = length(indexlist$PID)
  temp2 = matrix(NA, nrow = nr, ncol = ncol(beta))
  rownames(temp2) = names(indexlist$SID)
  colnames(temp2) = colnames(beta)
  for (i in 1:nr) {  
    if (i%%(round(nr/10)) == 0){
      print ("10% más...");      
    }
    temp = beta[indexlist$PID[[i]], ]
    if (length(indexlist$PID[[i]]) == 1) {
      temp2[i, ] = 100*1*is.na(temp);
      
    }
    else {
      if (indexmethod == "tbrm") {
        temp2[i, ] = apply(temp,2, function(x){100*sum(1*is.na(x))/(length(x))})
      }
      else {
        temp2[i, ] = apply(temp,2, function(x){100*sum(1*is.na(x))/(length(x))})
        
      }
    }
  }
  return(temp2)
}


porcentajeregionfunc_AL <- function (indexlist, beta, indexmethod = c("mean", "median", "tbrm")) 
{
  nr = 100;
  temp2 = matrix(NA, nrow = nr, ncol = ncol(beta))
  rownames(temp2) = names(indexlist$SID)[1:100]
  colnames(temp2) = colnames(beta)
    
  for (i in 1:nr) {
    temp = beta[indexlist$PID[[i]], ]
    if (length(indexlist$PID[[i]]) == 1) {
      
      temp2[i, ] = 100*1*is.na(temp);
    }
    else {
      if (indexmethod == "tbrm") {
        temp2[i, ] = apply(temp, 2, eval(indexmethod))
      }
      else {
        temp2[i,] = apply(temp,2, function(x){100*sum(1*is.na(x))/(length(x))})

      }
    }
  }
  return(temp2)
}

DMR.cpgs <- function(dmrcoutput,
                     annotation=c(array="IlluminaHumanMethylation450k",
                                  annotation="ilmn12.hg19"),
                     samps=NULL, toscale=FALSE, plotmedians=FALSE, ...)
{
  
  coords <- dmrcoutput$results$hg19coord[as.numeric(rownames(dmrs))]

  chr <- sub(":.*", "", coords)
  bookends <- sub(".*:", "", coords)
  startcpg <- as.integer(sub("-.*", "", bookends))
  stopcpg <- as.integer(sub(".*-", "", bookends))
  RSobject <- RatioSet(b, annotation=annotation)
  RSanno <- getAnnotation(RSobject)
  
  todas <-c()
  for (i in (1:length(startcpg))) {
    chri = chr[i]
    startcpgi = startcpg[i]
    stopcpgi = stopcpg[i]
    cpgs <- rownames(RSanno)[RSanno$chr %in% chri &
                               RSanno$pos >= startcpgi & RSanno$pos <= stopcpgi]
    todas <- c(todas,cpgs)
    
  }

  cpgs <- cpgs[order(RSanno[cpgs,"pos"])]
  cpgs
}


