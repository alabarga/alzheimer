	#	Load RCircos library
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	
	library(RCircos);


	#	Load human cytoband data and scatterplot data
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	data(RCircos.Scatter.Data);
	data(UCSC.HG19.Human.CytoBandIdeogram);
	cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;


	#	Setup RCircos core components:
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	RCircos.Set.Core.Components(cyto.info, NULL, 10, 0);


	#	Open the graphic device (here a pdf file)
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	
	# out.file <- "RCircos.Scatter.Plot.Demo.pdf";
	# pdf(file=out.file, height=8, width=8);

	RCircos.Set.Plot.Area();


	#	Draw chromosome ideogram
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	cat("Draw chromosome ideogram ...\n");

	RCircos.Chromosome.Ideogram.Plot();
	title("RCircos Scatter Plot Demo");
  
rcircos.params <- RCircos.Get.Plot.Parameters();
rcircos.params$track.height <- 0.3;
RCircos.Reset.Plot.Parameters(rcircos.params)


	#	Scatterplot 
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	scatter.data <- RCircos.Scatter.Data;
	scatter.colors <- rep("cyan", nrow(scatter.data));
	scatter.colors[which(scatter.data$seg.mean>=2)] <- "red";
	scatter.colors[which(scatter.data$seg.mean<=-2)] <- "blue";
	scatter.data["PlotColor"] <- scatter.colors;


	data.col <- 5;
	track.num <- 4; 
	RCircos.Scatter.Plot(scatter.data, data.col, 1, "in", 1);



	#	Close the graphic device
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	# dev.off();	print("RCircos Scatter Plot Demo Done!");

	# rm(list=ls(all=T));
