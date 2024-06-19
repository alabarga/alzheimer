DMR.heat <- function(dmrcoutput, dmr, betas, phen.col,
                     annotation=c(array="IlluminaHumanMethylation450k",
                                  annotation="ilmn12.hg19"),
                     samps=NULL, toscale=FALSE, plotmedians=FALSE, ...)
{
  
  stopifnot(is(dmrcoutput, 'dmrcate.output')) 
  stopifnot(is.matrix(betas))
  stopifnot((length(dmr) == 1) && (dmr %in% 1:nrow(dmrcoutput$results)))
  stopifnot(all(samps %in% 1:ncol(betas)))
  stopifnot(ncol(betas) == length(phen.col))
  if(dmrcoutput$results$no.probes[dmr] < 2){
    stop("Region must have 2 or more CpGs")
  }
  coords <- dmrcoutput$results$hg19coord[dmr]
  chr <- sub(":.*", "", coords)
  bookends <- sub(".*:", "", coords)
  startcpg <- as.integer(sub("-.*", "", bookends))
  stopcpg <- as.integer(sub(".*-", "", bookends))
  RSobject <- RatioSet(betas, annotation=annotation)
  RSanno <- getAnnotation(RSobject)
  cpgs <- rownames(RSanno)[RSanno$chr %in% chr &
                             RSanno$pos >= startcpg & RSanno$pos <= stopcpg]
  cpgs <- cpgs[order(RSanno[cpgs,"pos"])]
  cpgs
}



library(BSgenome.Hsapiens.UCSC.hg19)
seqinfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19)

DMR.genome<- function(lista, annot, b, clases) {

gtr <- GenomeAxisTrack()
itr <- IdeogramTrack(genome="hg19", chromosome="chr2")

gr2 <- GRanges(seqnames=Rle(paste('chr',as.character(annot[lista,'CHR']),sep='')),
               IRanges(start=annot[lista,'MAPINFO'], width=1), seqinfo=seqinfo, mcols=b[lista,])

dTrack <-DataTrack(gr2, groups = clases, col= c("red","blue"),
                   type = c("a", "p"), legend=TRUE)


gen = "hg19"
chr = as.character(seqnames(gr2)[1])
start = min(start(gr2)) - 20
end = max(start(gr2)) + 20

library(biomaRt)

biomartTrack <- BiomartGeneRegionTrack(genome=gen, chromosome=chr, start=start, end=end, name="ENSEMBL", showId=TRUE)

martfunc <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
ensfunc <- 
  getBM(c("ensembl_gene_id","ensembl_transcript_id","exon_chrom_start",
          "exon_chrom_end","strand","gene_biotype","external_gene_name"),
           filters = c("chromosome_name", "start", "end"),
           values = list(chr, start, end), mart=martfunc)

biomartTrack <- AnnotationTrack(chr=chr,strand=ensfunc[,5],start=ensfunc[,3],end=ensfunc[,4],
                                  feature=ensfunc[,6],group=ensfunc[,1],id=ensfunc[,7],
                                  name = "genes")

biomartTrack <- BiomartGeneRegionTrack(genome=gen, chromosome=chr, start=start, end=end, name="genes", stacking='squish', showId=TRUE)

plotTracks(list(itr, gtr, biomartTrack, dTrack))
plotTracks(list(itr, gtr, biomartTrack, dTrack), from = start, to = end)

}

explore450 <- function(meth450_data, group){
  
  b <- betas(meth450_data)
  clases = phenoData(meth450_data)[[group]]
  plotMDS(b, labels=colnames(meth450_data), col=as.integer(clases))
}

metadata<- function(meth450_data,lista){
  
  metadatos <- (meth450_data@featureData@data[lista,c("CHR","UCSC_REFGENE_NAME","UCSC_REFGENE_GROUP")])
  
  metadatos <- metadatos[order(metadatos$CHR),]
  
  metadatos <- metadatos[metadatos$UCSC_REFGENE_NAME != "",]
  
  View(metadatos)
  
  metadatos
}

profilePlot<-function(b,lista,clases){
  
  
  sdata <- b[lista,order(clases)]
  sdata <- melt(sdata)
  
  sdata <- transform(sdata,X2=factor(Var2,levels=unique(Var2)))
  profile_plot <- ggplot(sdata, aes(x=Var2, y=value, group=Var1)) + geom_line()
  profile_plot
  
  #geom_line(aes(colour = fc)) + scale_colour_gradient(low="red")
  #geom_point(aes, size = 2.5) +
  #  geom_line(aes, size = 1)
  
  #require(plotly)
  #py <- plotly(username="alabarga", key="a935mlupe5")
  #py$ggplotly(profile_plot)
  
  #require(rCharts)
  #rPlot(value ~ as.factor(X2), color='Var1', data=sdata, type='line')
  
}


qcPlots <- function(betas, samples, clases = NULL){
  
  require(bioDist)
  eset = na.omit(betas)
  
  hc1 <- hclust(cor.dist(t(eset)), method = "average")
  
  if (nlevels(clases) > 2) {
    cols <- brewer.pal(nlevels(clases),'Set3')
  } else {
    cols <- c("red","blue")
  }
  
  plotColoredClusters(hc1, labs=samples, cols=cols[clases])  
  
  plot(hc1, samples, xlab = "Sample", main = "Clustering samples by all the CpG loci ", 
       lwd = 2, font.axis = 2, font.lab = 2)
  
  boxplot(betas, ylab = "beta Value")  
}

hcPlot <- function(betas, samples, clases = NULL){
  
  require(bioDist)
  eset = na.omit(betas)
  
  hc1 <- hclust(cor.dist(t(eset)), method = "average")
  
  if (nlevels(clases) > 2) {
    cols <- brewer.pal(nlevels(clases),'Set3')
  } else {
    cols <- c("red","blue")
  }
  
  plotColoredClusters(hc1, labs=samples, cols=cols[clases])  
  
}

pvalPlot<- function(pvalues){
  avgPval = apply(pvalues, 2, function(x) {
    sum(x >= 1e-05) * 100/length(x)
  })
  barplot(avgPval, ylab = "% of detect pvalue >1e-5")
}

b_heatmap<- function(betas,lista,clases)
{
  require(RColorBrewer)
  require(pheatmap)
  
  so <- sort(as.character(clases),index.return=T)
  mat <- betas[lista,]
  
  pheatmap(mat, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), cluster_rows=TRUE)
}

b_heatmap.2<- function(betas,lista,clases)
{
  require(RColorBrewer)
  require(pheatmap)
  
  col1 <- brewer.pal(8, "Set1")
  col.colors = col1[clases]
  so <- sort(as.character(clases),index.return=T)
  mat <- betas[lista,]
  
  heatmap.2(mat,  ColSideColors=col.colors)
}

b_heatmap.3<- function(betas,lista,clases)
{
  
  #install.packages("made4")
  require(made4)
  
  #install.packages('impute')
  require(impute)
  
  mat <- as.matrix(betas[lista,order(clases)])
  mat <- impute.knn(mat)$data
  
  col1 <- brewer.pal(8, "Set1")
  col.colors = col1[sort(clases)]
  #heatplot(mat, cols.default=FALSE, lowcol="white", highcol="red", ColSideColors=col.colors,lhei=c(1, 2))
  heatplot(mat, ColSideColors=col.colors,lhei=c(1, 2),dend = "row")
  
  par(lend = 1)           # square line ends for the color legend
  legend("topleft",      # location of the legend on the heatmap plot
         legend = levels(clases), # category labels
         col = col1[1:nlevels(clases)],  # color key
         lty= 1,             # line style
         lwd = 10            # line width
  )
  
}

compare.lists <- function(lista1,lista2) {
  draw.pairwise.venn(length(lista1),length(lista2),length(intersect(lista1,lista2)),category=c("lista1","lista2"))
  compare.lists <- intersect(lista1,lista2)
}


volcano450 <- function(res, adjust.p = 0.1, diff.thresh=0.1, pval.var= "adj.P.Val", diff.var="Beta.Difference", names= TRUE, lista = NULL) 
  {
  
  if (is.null(lista))  lista <- rownames(res) 
     
   mx <- max(abs(res[,diff.var]))
   my <- max(-log10(res[,"P.Value"]))
   
  res$color[(res[,diff.var] > diff.thresh) & (res[,pval.var] < adjust.p) & (rownames(res) %in% lista)] = "D+"
  res$color[(res[,diff.var] < (- diff.thresh)) & (res[,pval.var] < adjust.p) & (rownames(res) %in% lista)] = "D-"
  res$color[(abs(res[,diff.var]) < diff.thresh) | (res[,pval.var] > adjust.p) | !(rownames(res) %in% lista)]= "NoDM"
  
  res$size[(res[,diff.var] > diff.thresh) & (res[,pval.var] < adjust.p) & (rownames(res) %in% lista)] = 2
  res$size[(res[,diff.var] < (- diff.thresh)) & (res[,pval.var] < adjust.p) & (rownames(res) %in% lista)] = 2
  res$size[(abs(res[,diff.var]) < diff.thresh) | (res[,pval.var] > adjust.p) | !(rownames(res) %in% lista)]= 1.75
  
  
  if (diff.var=="Beta.Difference") {
  
  volcano = ggplot(data=res, aes(x=Beta.Difference, y=-log10(P.Value), colour=color, size=size)) +
    scale_colour_manual(values=c("#56B4E9", "#FF9999", "#888888")) +
    geom_point(alpha=0.4) +
    theme(legend.position = "none") +
    xlim(c(-mx, mx)) + ylim(c(0, my)) +
    xlab("beta difference") + ylab("-log10 p-value") 
  
   if (names) {
     
     volcano = volcano +
       geom_text(data=subset(res, as.logical(res$color == "D+")), aes(x=Beta.Difference, y=-log10(P.Value),
                                                                      label=geneName, size=1.2), colour="#FF9999") +
       geom_text(data=subset(res, as.logical(res$color == "D-")), aes(x=Beta.Difference, y=-log10(P.Value),
                                                                      label=geneName, size=1.2), colour="#56B4E9")
     
   }
    
  } else {
  volcano = ggplot(data=res, aes(x=logFC, y=-log10(adj.P.Val), colour=color)) +
    scale_colour_manual(values=c("#56B4E9", "#FF9999", "#888888")) +
    geom_point(alpha=0.4, size=1.75) +
    theme(legend.position = "none") +
    xlim(c(-mx, mx)) + ylim(c(0, my)) +
    xlab("beta difference") + ylab("-log10 p-value") 
  
  if (names) {
    
    volcano = volcano +
    geom_text(data=subset(res, as.logical(res$color == "D+")), aes(x=logFC, y=-log10(adj.P.Val),
                                                                   label=geneName, size=1.2), colour="#FF9999") +
    geom_text(data=subset(res, as.logical(res$color == "D-")), aes(x=logFC, y=-log10(adj.P.Val),
                                                                   label=geneName, size=1.2), colour="#56B4E9")
  }
  
  }
  plot(volcano)
  
  which(res$color %in% c("D+","D-"))
}


beta_FC <- function(res, adjust.p = 0.1, diff.thresh=0, pval.var= "adj.P.Val", diff.var="Beta.Difference") {
  
  mx <- max(abs(res[,"Beta.Difference"]))
  my <- max(abs(res[,"logFC"]))
  
  res$color[(res[,diff.var]) > diff.thresh & res[,pval.var] < adjust.p] = "D+"
  res$color[(res[,diff.var]) < - diff.thresh & res[,pval.var] < adjust.p] = "D-"
  res$color[abs(res[,diff.var]) < diff.thresh | res[,pval.var] > adjust.p] = "NoDM"
  
  volcano = ggplot(data=res, aes(x=Beta.Difference, y=logFC, colour=color)) +
      scale_colour_manual(values=c("#56B4E9", "#FF9999", "#888888")) +
      geom_point(alpha=0.4, size=1.75) +
      theme(legend.position = "none") +
      xlim(c(-mx, mx)) + ylim(c(-my, my)) +
      xlab("beta difference") + ylab("-log10 p-value") +
      geom_text(data=subset(res, as.logical(res$color == "D+")), aes(x=Beta.Difference, y=logFC,
                                                                     label=geneName, size=1.2), colour="#FF9999") +
      geom_text(data=subset(res, as.logical(res$color == "D-")), aes(x=Beta.Difference, y=logFC,
                                                                     label=geneName, size=1.2), colour="#56B4E9")
 
  plot(volcano)
  
}