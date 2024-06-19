DMR.plot.region <- function(dmrcoutput, chr, startcpg, stopcpg, betas, phen.col,
                     annotation=c(array="IlluminaHumanMethylation450k",
                                  annotation="ilmn12.hg19"),
                     samps=NULL, toscale=FALSE, plotmedians=FALSE, ...)
{
  
 
  #coords <- dmrcoutput$results$hg19coord[dmr]
  #chr <- sub(":.*", "", coords)
  #bookends <- sub(".*:", "", coords)
  #startcpg <- as.integer(sub("-.*", "", bookends))
  #stopcpg <- as.integer(sub(".*-", "", bookends))
  
  RSobject <- RatioSet(betas, annotation=annotation)
  RSanno <- getAnnotation(RSobject)
  cpgs <- rownames(RSanno)[RSanno$chr %in% chr &
                             RSanno$pos >= startcpg & RSanno$pos <= stopcpg]
  cpgs <- cpgs[order(RSanno[cpgs,"pos"])]
  if(is.null(samps)){samps <- 1:ncol(betas)}
  betas <- betas[as.character(cpgs), samps]
  m <- match(as.character(cpgs), rownames(RSanno))
  clusta <- data.frame(
    gene=RSanno$UCSC_RefGene_Name[m],
    group=RSanno$UCSC_RefGene_Group[m], pos=RSanno$pos[m])
  if(!toscale) {
    plot(1:nrow(clusta), , type='n', xlab=paste(chr, "consecutive probes"),
         ylab="beta values", ylim=c(0,1.15), ...)
    for(i in 1:nrow(clusta)) {
      points(rep(i, ncol(betas)), betas[i,], col=phen.col, ...)
    } 
    
    if(plotmedians){
      medians <- matrix(0, nrow(clusta), ncol=length(unique(phen.col)))
      colnames(medians) <- unique(phen.col)
      for(j in colnames(medians)) {
        for (i in 1:nrow(clusta)){
          medians[i,j] <- median(betas[i,phen.col==j])
        }
        lines(1:nrow(clusta), medians[,j], col=j, lwd=2)
      }
    }
    
    abline(1, 0, col='gray')
    abline(0.5, 0, col='gray')
    abline(0, 0, col='gray')
    empty <- rep(-1, nrow(clusta))
    empty[grep("5'UTR", clusta$group)] <- 1
    points(1:length(empty), empty, adj=0.6, col='blue4', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("3'UTR", clusta$group)] <- 1.01
    points(1:length(empty), empty, adj=0.8, col='cyan', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("Body", clusta$group)] <- 1.02
    points(1:length(empty), empty, col='firebrick1', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("TSS1500", clusta$group)] <- 1.03
    points(1:length(empty), empty, adj=0.2, col='chartreuse', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("TSS200", clusta$group)] <- 1.04
    points(1:length(empty), empty, adj=0.4, col='forestgreen', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("1stExon", clusta$group)] <- 1.05
    points(1:length(empty), empty, adj=1, col='magenta', pch=15)
    empty <- rep('', nrow(clusta))
    geneidxs <- cumsum(rle(as.character(clusta$gene))$lengths)
    geneidxs <- ceiling(c(geneidxs[1]/2, geneidxs[2:length(geneidxs)] -
                            diff(geneidxs)/2))
    empty[geneidxs] <- as.character(clusta$gene[geneidxs])
    genevector <- empty
    text(1:nrow(clusta), 1.03, genevector, pos=3, offset=1)
  }
  if (toscale) {
    plot(clusta$pos, 1:nrow(clusta), type='n',
         xlab=paste(chr, "hg19 coords"), ylab="beta values",
         ylim=c(0,1.15), ...)
    for(i in 1:nrow(clusta)){
      points(rep(clusta$pos[i], ncol(betas)), betas[i,], col=phen.col, ...)
    } 
    
    if(plotmedians){
      medians <- matrix(0, nrow(clusta), ncol=length(unique(phen.col)))
      colnames(medians) <- unique(phen.col)
      for(j in colnames(medians)) {
        for (i in 1:nrow(clusta)){
          medians[i,j] <- median(betas[i,phen.col==j])
        }
        lines(clusta$pos, medians[,j], col=j, lwd=2)
      }
    }
    
    abline(1, 0, col='gray')
    abline(0.5, 0, col='gray')
    abline(0, 0, col='gray')
    empty <- rep(-1, nrow(clusta))
    empty[grep("5'UTR", clusta$group)] <- 1
    points(clusta$pos, empty, adj=0.6, col='blue4', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("3'UTR", clusta$group)] <- 1.01
    points(clusta$pos, empty, adj=0.8, col='cyan', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("Body", clusta$group)] <- 1.02
    points(clusta$pos, empty, col='firebrick1', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("TSS1500", clusta$group)] <- 1.03
    points(clusta$pos, empty, adj=0.2, col='chartreuse', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("TSS200", clusta$group)] <- 1.04
    points(clusta$pos, empty, adj=0.4, col='forestgreen', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("1stExon", clusta$group)] <- 1.05
    points(clusta$pos, empty, adj=1, col='magenta', pch=15)
    empty <- rep('', nrow(clusta))
    geneidxs <- cumsum(rle(as.character(clusta$gene))$lengths)
    geneidxs <- ceiling(c(geneidxs[1]/2, geneidxs[2:length(geneidxs)] -
                            diff(geneidxs)/2))
    empty[geneidxs] <- as.character(clusta$gene[geneidxs])
    genevector <- empty
    text(clusta$pos, 1.03, genevector, pos=3, offset=1)
  }
}

data <- read.table(paste("data","lincRNAs_transcripts.bed",sep="/"),header=F)
data <- data[,1:6]
colnames(data) <- c('chr','start','end','id','score','strand')
lincrna.bed <- with(data, GRanges(chr, IRanges(start, end), strand, score, id=id))


coords <- dmrcoutput$results$hg19coord
chr <- sub(":.*", "", coords)
bookends <- sub(".*:", "", coords)
startcpg <- as.integer(sub("-.*", "", bookends))
stopcpg <- as.integer(sub(".*-", "", bookends))
RSobject <- RatioSet(b, annotation=annotation)
RSanno <- getAnnotation(RSobject)

get_cpgs<-function(data){
  class(data)
  cpgs <- rownames(RSanno)[RSanno$chr %in% data[1] &
                             RSanno$pos >= as.integer(data[2]) & RSanno$pos <= as.integer(data[3])]
  cpgs <- cpgs[order(RSanno[cpgs,"pos"])]
  cpgs  
}

dmrs <- data.frame(chr=chr,startcpg=startcpg,stopcpg=stopcpg, stringsAsFactors=FALSE)

cpgs <- apply(dmrs,1, get_cpgs)
cpgs <- unlist(cpgs)

cpgs_in_dmr <- cbind(id=cpgs, annot[cpgs,c('CHR','MAPINFO')],betadiff=betaDiff(b[cpgs,],clases))
dmp.bed <- with(cpgs_in_dmr, GRanges(paste("chr",CHR,sep=""), IRanges(MAPINFO, MAPINFO+1), strand='*', score=0, id=id))

oo<- subsetByOverlaps(dmp.bed,lincrna.bed)

cpgs_total <- cbind(id=rownames(annot), annot[,c('CHR','MAPINFO')],betadiff=betaDiff(b,clases))
cpgs.bed <- with(cpgs_total, GRanges(paste("chr",CHR,sep=""), IRanges(MAPINFO, MAPINFO+1), strand='*', score=0, id=id))

oo<- subsetByOverlaps(cpgs.bed,lincrna.bed)

data <- read.table(paste("data","wgEncodeBroadHmmH1hescHMM.bed",sep="/"),header=F)
data <- data[,c(1:6,9)]
colnames(data) <- c('chr','start','end','type','score','strand','color')

table(data$type)
c <- ggplot(data, aes(factor(type), fill=factor(type))) + geom_bar() + coord_flip()
c

poised <- subset(data,type=='3_Poised_Promoter')

#1_Active_Promoter 10_Txn_Elongation       11_Weak_Txn      12_Repressed 13_Heterochrom/lo 
#12475             17027            118002             17589             91228 
#14_Repetitive/CNV 15_Repetitive/CNV   2_Weak_Promoter 3_Poised_Promoter 4_Strong_Enhancer 
#3716              2439             31965             11570              5264 
#5_Strong_Enhancer   6_Weak_Enhancer   7_Weak_Enhancer       8_Insulator  9_Txn_Transition 
#12812             84576            137484             60691             12223 

h1hesc.bed <- with(poised, GRanges(chr, IRanges(start, end), strand="*", score, id=paste(chr,":",start,"-",end,sep="")))

oo<- subsetByOverlaps(dmp.bed,h1hesc.bed)


#http://stattrek.com/online-calculator/hypergeometric.aspx
# Population size  
# Number of successes in population	
# Sample size	
# Number of successes in sample (x)	


n_cpgs <- length(cpgs.bed)
n_cpgs_in_dmr <- length(dmp.bed)
o_cpgs <- length(subsetByOverlaps(cpgs.bed,h1hesc.bed))
o_dmr <- length(subsetByOverlaps(dmp.bed,h1hesc.bed))

phyper(o_dmr-1, o_cpgs, n_cpgs - o_cpgs, n_cpgs_in_dmr, lower.tail=FALSE)

