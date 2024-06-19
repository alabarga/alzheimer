MethylDataFileName = "./data/multiple/GSE40360_GenomeStudio.txt"
SampleDescFileName = "./data/multiple/GSE40360_SampleDescription.csv"

group = "group"
gcase = c("EA1", "EA2", "EA3")   ### Specify the case group index in the sample.txt file (if "concov" is "ON")
gcontrol = c("EA_CONTROL")       ### Specify the control group index in the sample.txt file (if "concov" is "ON")

#### Filtramos muestras erroneas o que no pertenecen al estudio
remove_samples = c("46 ER","57 HC")

projectName = 'multiple'

clases <- as.factor(c("CONTROL","MS")[samps$disease_status])

require(biomaRt)
genes <- c("CYP2R1","CYP27B1", "CYP24A1", "VDR", "RXRB" )
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")

filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
View(attributes)
View(filters)

res <- getBM(attributes=c('hgnc_symbol', 'chromosome_name', 'strand', 'start_position', 'end_position'), filters = 'hgnc_symbol', values = genes, mart = ensembl)

res <- res[-c(3,5:9),]

library(GenomicRanges)


genes_gr <- GRanges(seqnames = paste("chr",as.character(res$chromosome_name),sep=""), ranges = IRanges(start = res$start_position - 5000, end = res$end_position + 5000), names = res$hgnc_symbol)

mi <- meth450_data_original@featureData@data
mi <- mi[,c("TargetID", "CHR","MAPINFO", "UCSC_REFGENE_NAME")]

mi <- mi[(mi$CHR %in% res$chromosome),]
mi <- mi[!is.na(mi$MAPINFO),]

probes_gr <- GRanges(seqnames = paste("chr",as.character(mi$CHR),sep=""), ranges = IRanges(start = mi$MAPINFO, end = mi$MAPINFO + 50), names = mi$TargetID)
overlaps<- findOverlaps(probes_gr, genes_gr, type="within")
mi[queryHits(overlaps),]

res[subjectHits(overlaps),]

lista_cpgs <- rownames(mi[queryHits(overlaps),])

library(Gviz)
data(cpgIslands)
atr <- AnnotationTrack(cpgIslands, name="CpG")
gtr <- GenomeAxisTrack()
itr <- IdeogramTrack(genome="mm9", chromosome="chr7")
data(geneModels)
grtr <- GeneRegionTrack(geneModels, name="Gene Model", showId=TRUE)
plotTracks(list(itr, gtr, atr, grtr))

for (x in unique(subjectHits(overlaps))) {
  cpgs_ix <- queryHits(overlaps)[subjectHits(overlaps) == x]
  cpgs <- probes_gr[cpgs_ix]
  region <- genes_gr[x]
  
  chr = as.character(seqnames(region))
  atr <- AnnotationTrack(cpgs, name="CpG")
  itr <- IdeogramTrack(genome="hg19", chromosome=chr)
  biomartTrack <- BiomartGeneRegionTrack(biomart=ensembl, genome="hg19", chromosome=chr, start=start(region), end=end(region), fontcolor="black",showId=TRUE,size=2)
  gtr <- GenomeAxisTrack()
  plotTracks(list(itr, gtr, atr, biomartTrack))
}
