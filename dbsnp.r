#ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv
#ftp://hgdownload.soe.ucsc.edu//apache/htdocs/goldenPath/hg19/database/snp138Common.txt.gz
#http://bit.ly/dbsnp137-schema
#ftp://hgdownload.soe.ucsc.edu//apache/htdocs/goldenPath/hg19/database/snp141Common.txt.gz
#gunzip -c snp141Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp141Common_small.txt.gz

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

processUCSCsnp <- function(snpfile) {
  require(GenomicRanges)
  cat("Reading file\n")
  df <- read.delim(gzfile(snpfile), header = FALSE,
                   stringsAsFactors = FALSE)
  names(df) <- c("chr", "start", "end", "name", "strand",
                 "refNCBI", "class", "alleleFreqs")
  print(table(df$chr))
  cat("Only keeping chrs 1-22, X, Y\n")
  df <- df[df$chr %in% paste0("chr", c(1:22, "X", "Y")),]
  print(table(df$class))
  cat("Only keeping class 'single'\n")
  df <- df[df$class == "single",]
  cat("Computing MAF\n")
  df$alleleFreqs <- sub(",$", "", df$alleleFreqs)
  sp <- strsplit(df$alleleFreqs, ",")
  minFreq <- sapply(sp, function(xx) min(as.numeric(xx)))
  cat("Instatiating object\n")
  grSNP <- GRanges(seqnames = df$chr, strand = df$strand, ranges = IRanges(start = df$start + 1, end = df$end),
                   MAF = minFreq, ref = df$refNCBI)
  names(grSNP) <- df$name
  grSNP
}

grSnp141CommonSingle <- processUCSCsnp("./data/snp141Common_small.txt.gz")
# gg <- grSnp138CommonSingle[(grSnp138CommonSingle$MAF <= 0.95) & (grSnp138CommonSingle$MAF >= 0.05),]
save(grSnp141CommonSingle, file=paste("./data","grSnp141CommonSingle.Rdata",sep="/"))

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

Locations <- anno[,c("chr","pos","strand")]
Locations <- as(Locations, "DataFrame")
Manifest <- anno[, c("Name", "AddressA", "AddressB",
                     "ProbeSeqA", "ProbeSeqB", "Type", "NextBase", "Color")]
Manifest <- as(Manifest, "DataFrame")

map <- cbind(Locations, Manifest)
map <- GRanges(seqnames = map$chr, ranges = IRanges(start = map$pos, width = 1),
               Strand = map$strand, Type = map$Type)

map <- minfi:::getProbePositionsDetailed(map)
names(map) <- rownames(Locations)

SNPs.141 <- minfi:::.doSnpOverlap(map, grSnp141CommonSingle)
save(SNPs.141, file=paste("./data","SNPs.141.Rdata",sep="/"))

map <- cbind(Locations, Manifest)
map <- GRanges(seqnames = map$chr, ranges = IRanges(start = map$pos -10, width = 20),
               Strand = map$strand, Type = map$Type)

SNP_10 = findOverlaps(map, grSnp141CommonSingle,ignore.strand=T)
save(SNP_10, file=paste("./data","SNP_10.Rdata",sep="/"))
