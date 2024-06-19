#  Load RCircos library
#  	_________________________________________________________________
#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

library(RCircos);


#	Load human cytoband data 
#  	_________________________________________________________________
#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

data(UCSC.HG19.Human.CytoBandIdeogram);
hg19.cyto <- UCSC.HG19.Human.CytoBandIdeogram;


#	Setup RCircos core components:
#  	_________________________________________________________________
#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# RCircos.Set.Core.Components(cyto.info, NULL, 5, 5);
RCircos.Set.Core.Components(cyto.info=hg19.cyto, chr.exclude=NULL, 
                            tracks.inside=8, tracks.outside=0 );

params <- RCircos.Get.Plot.Parameters();
params$track.height = 0.1
RCircos.Reset.Plot.Parameters(params)

#	Open the graphic device (here a pdf file)
#
#	png(file="RCircos.Layout.Demo.png", height=8, width=8, unit="in", 
#		type="cairo", res=300);
#  	_________________________________________________________________
#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

out.file <- "RCircos.Methyl450.Meditation.pdf";
pdf(file=out.file, height=8, width=8);

RCircos.Set.Plot.Area();
title("Methylation 450K analysis :: meditacion");

#	Draw chromosome ideogram
#  	_________________________________________________________________
#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cat("Draw chromosome ideogram ...\n");
RCircos.Chromosome.Ideogram.Plot();

#  Connectors in first track and gene names in the second track. 
#  	_________________________________________________________________
#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cat("Add Gene and connector tracks ...\n");

RCircos.Gene.Label.Data <- table_selected[,c("CHR","MAPINFO","MAPINFO","geneName")]
colnames(RCircos.Gene.Label.Data) <- c("Chromosome", "chromStart", "chromEnd", "Gene")
RCircos.Gene.Connector.Plot(genomic.data=RCircos.Gene.Label.Data, 
                            track.num=1, side="in");

RCircos.Gene.Name.Plot(gene.data=RCircos.Gene.Label.Data, name.col=4, 
                       track.num=2, side="in");

RCircos.Heatmap.Data <- table_selected[,c("CHR","MAPINFO","MAPINFO","Beta.Difference")]

RCircos.Heatmap.Plot(heatmap.data=RCircos.Heatmap.Data, data.col=4, 
                     track.num=4, side="in");

params <- RCircos.Get.Plot.Parameters();
params$track.height = 0.2
RCircos.Reset.Plot.Parameters(params)

RCircos.Scatter.Data <- table_selected[,c("CHR","MAPINFO","MAPINFO","Beta.Difference")]

#RCircos.Scatter.Data <- table_selected[,c("CHR","MAPINFO","MAPINFO")]
#RCircos.Scatter.Data$DATA <- table_selected$Beta.Difference

RCircos.Scatter.Plot(scatter.data=RCircos.Scatter.Data, data.col=4, 
                       track.num=3, side="in", by.fold=0.1);

#RCircos.Histogram.Plot(hist.data=RCircos.Scatter.Data, data.col=4, 
#                       track.num=3, side="in");

#locations <- RCircos.Track.Positions("in", 5)
#RCircos.Track.Outline(locations[1], locations[2],1)

#RCircos.Histogram.Data <- res[, c("CHR","MAPINFO","MAPINFO")]
#RCircos.Histogram.Data$DATA <- -log(res$adj.P.Val)
#RCircos.Histogram.Plot(hist.data=RCircos.Histogram.Data, data.col=4, 
#                       track.num=6, side="in");


#	Close the graphic device and clear memory
#  	_________________________________________________________________
#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

dev.off();
print("RCircos Layout Demo Done!");

#rm(list=ls(all=T));