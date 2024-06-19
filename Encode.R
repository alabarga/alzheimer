library(rtracklayer)
library(GenomicRanges)

fmt <- function(){
  f <- function(x) as.character(round(x,2))
  f
}

select_for_encode<- function(meth450_data_original, data.test.M, sel_cghs){
  all_cghs <- rownames(data.test.M)
  mapinfo <- meth450_data_original@featureData@data[all_cghs,c("CHR","MAPINFO")]
  seleccion_all <- cbind(Id=all_cghs,BetaDifference=data.test.M[,"Beta-Difference"],mapinfo)
  
  mapinfo <- meth450_data_original@featureData@data[sel_cghs,c("CHR","MAPINFO")]
  seleccion <- cbind(Id=sel_cghs,BetaDifference=data.test.M[sel_cghs,"Beta-Difference"],mapinfo)
  
  seleccion_ranges<- with(seleccion, GRanges(seqnames=paste("chr",CHR,sep=""), ranges = IRanges(as.integer(MAPINFO)-1, as.integer(MAPINFO)), strand="*", score = BetaDifference ))
  seleccion_all_ranges<- with(seleccion_all, GRanges(seqnames=paste("chr",CHR,sep=""), ranges = IRanges(as.integer(MAPINFO)-200, as.integer(MAPINFO)+200), strand="*", score = BetaDifference ))
  
  list(seleccion_ranges=seleccion_ranges, seleccion_all_ranges=seleccion_all_ranges)
}


histones_report <- function(ranges, broad_dir = ".", report_file = "histones_report.csv", type_histones = c(rep("H1 hesc",16),rep("Nha",13))) {
  
  seleccion_ranges <- ranges$seleccion_ranges
  seleccion_all_ranges <- ranges$seleccion_all_ranges
  
  files_histones <- list.files(path = broad_dir, pattern="total_pooled\\.broadPeak$",recursive=TRUE,include.dirs=TRUE)
  
  report <- data.frame(ratio_rel=numeric(), ratio_sel_denominador=numeric(),length_seleccion_ranges=numeric(),length_seleccion_all_ranges=numeric(), stringsAsFactors=FALSE)
  for(file_peaks_histones in files_histones){
    
    #Generamos el Rango con los picos del ENCODE
    #Directamente con el fichero broadPeak descargado. No se necesita modificar pero comprobar que el separador decimal es el . y no ,
    #http://www.genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeBroadHistone
    
    peaks_histones <- read.delim(file_peaks_histones,sep="\t")
    colnames(peaks_histones)<-c("chrom","chromStart","chromEnd","name","score","strand","signalValue","pValue","qValue")
    
    ## #Seleccionamos mejores
    ## filtro = 8
    ## #filtro = 0
    ## peaks_histones<-peaks_histones[peaks_histones$signalValue >= filtro,]
    peak_histones_ranges <- with(peaks_histones, GRanges(seqnames=chrom, ranges = IRanges(chromStart, chromEnd), strand="*", score = signalValue, signalValue= signalValue, pValue=-1, qValue=-1))
    
    overlaps_sel <- mean(peak_histones_ranges[countOverlaps(peak_histones_ranges,seleccion_ranges)>0]$signalValue)
    overlaps_all <- mean(peak_histones_ranges[countOverlaps(peak_histones_ranges,seleccion_all_ranges)>0]$signalValue)
    
    overlaps_sel_denominador <- sum(peak_histones_ranges[countOverlaps(peak_histones_ranges,seleccion_ranges)>0]$signalValue)/length(seleccion_ranges)
    overlaps_all_denominador <- sum(peak_histones_ranges[countOverlaps(peak_histones_ranges,seleccion_all_ranges)>0]$signalValue)/length(seleccion_all_ranges)
    
    ratio_sel <- overlaps_sel/overlaps_all
    
    ratio_sel_denominador <- overlaps_sel_denominador/overlaps_all_denominador
    
    report[nrow(report) + 1, 1:4] <- c(file_peaks_histones,ratio_sel,ratio_sel_denominador,length(seleccion_ranges),length(seleccion_all_ranges))
    
  }
  
  report$filename <- files_histones
  report$type <- type_histones
  report$name <- gsub("H1hesc|Nha", "", gsub("Pk.broadPeak/total_pooled.broadPeak","",gsub("broad/wgEncodeBroadHistone","",report$filename)))
  write.csv(report, file=report_file)
  
  ggplot(report, aes(x = name, y = ratio_rel,fill=name)) +
    geom_bar(stat="identity") +
    facet_wrap( ~ type) +
    scale_y_continuous(labels = fmt()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  report
}