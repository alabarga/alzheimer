
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

sapply(colnames(meth450_data_original), function(x){print predictSex(meth450_data_original,x)})
