##########################################
############## tDMR analysis #############
##########################################
##########################################
### Script written by Roderick Slieker ###
### Molecular Epidemiology             ###
### Leiden University Medical Center   ###
### r.c.slieker@lumc.nl                ###
##########################################

DMRfinder = function(data, mismatches=3, icd=1000, illumina=TRUE){
  
  #Check for annotation  
  if(illumina==TRUE){
    check  = length(grep("cg", data[,1]))> 1
    check2  = length(grep("0", data[,1]))> 1
  
    if(check==FALSE | check2==FALSE){
      print("Please check your data, is the format right?")
      }else{
        library("FDb.InfiniumMethylation.hg19")
        InfiniumMethylation <- features(FDb.InfiniumMethylation.hg19)
        probesselect= as.data.frame(data, stringsAsFactors = F)
        IM = InfiniumMethylation[names(InfiniumMethylation) %in% probesselect[,1],]
        IMx  =   cbind(as.character(seqnames(IM)),start(IM),end(IM))
        rownames(IMx) = names(IM)
        IMx = as.data.frame(IMx[probesselect[,1],],stringsAsFactors=F)
        probesselect$chr = IMx$V1
        probesselect$coord  = IMx$V2
        probesselect = probesselect[,c(3,4,2)]
        }
  }else{
    probesselect = as.data.frame(data)
  }
  
  pb = txtProgressBar(min=1,max=22,style=3)
  
  #DMR finder loop per chromosome
  DMR = function(x,data,mismatches,icd){
    setTxtProgressBar(pb,x)
    close(pb)
    tryCatch.W.E <- function(expr)
    {
      W <- NULL
      w.handler <- function(w){ # warning handler
        W <<- w
        invokeRestart("muffleWarning")
      }
      list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                       warning = w.handler),
           warning = W)
    }
    x = paste("chr",x,sep="")
    chr1 = probesselect[probesselect$chr==x,]
    order = order(as.numeric(chr1$coord))
    chr.sorted=chr1[order,]
    chr.final = data.frame(coord = as.numeric(chr.sorted$coord) , crit = as.numeric(chr.sorted[,3]))
    file = "temp.txt"
    sink(file)
    MAXIMUM_REGION_LENGTH = icd # constant, can be adjusted
    last_coordinate = length( chr.final$crit )
    next_coordinate = 0 # so that a region that has been called will be skipped
    
    for (i in 1:(last_coordinate-1)) {
      if ( i>=next_coordinate ) {
        if (chr.final$crit[ i ]==1) {
          start_location = chr.final$coord[ i ]
          last_visited_crit_loc = start_location
          sum_of_ones = 1
          number_of_items = 1
          
          # start crawling loop
          for (j in (i+1):last_coordinate ) {
            if (chr.final$coord[ j ] > (last_visited_crit_loc + MAXIMUM_REGION_LENGTH)) { break }
            if((number_of_items-sum_of_ones)>mismatches) { break }   #Number of mismatches
            number_of_items = number_of_items + 1
            if (chr.final$crit[j]==1) { 
              last_visited_crit_loc = chr.final$coord[ j ]
              sum_of_ones = sum_of_ones + 1 
            }
          }
          
          # now check if the result is good enough
          if (sum_of_ones>=3) {
            last_one=i+number_of_items-1
            for (k in (i+number_of_items-1):1) {
              if ( chr.final$crit[k] == 0 ) {
                last_one = last_one - 1
                number_of_items = number_of_items - 1
              }
              else {
                break
              }
            }
            cat(start_location,";",chr.final$coord[last_one],";",sum_of_ones/number_of_items,"\n")
            next_coordinate = last_one + 1
          }
        }
      }
    }
    sink()
    check = tryCatch.W.E(read.table("temp.txt" , sep=";"))
    check.for.error = class(check$value)[1] == "simpleError"
    if(check.for.error==TRUE){
      
    }else{
      
      dmr = as.matrix(check$value)
      dmr.chr = cbind(rep(x),dmr[,1:2])
      dmr.chr
    }
  }
  
  #Run DMR finder for all chromosomes
  dmr.func = sapply(seq(1:22),DMR,probesselect,mismatches,icd)
  dmr.allchr = do.call(rbind , dmr.func)
  dmrs.retGR = GRanges(seqnames=dmr.allchr[,1],IRanges(as.numeric(dmr.allchr[,2]),as.numeric(dmr.allchr[,3])))
  dmrs.retGR
}
