# Preprocessing and Quality Control Functions
# "The Nuts and Bolts of DNA Methylation Array Analysis Workshop"
# Devin C. Koestler Ph.D.
# July 11th 2013

# SUMMARY:
#######################################################################
# The code described below implements several diagnostic plots for 
# identifying poor performing samples
######################################################################

# install necessary packages 
#install.packages("RColorBrewer")
library(RColorBrewer)
#display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE)

# Multidimensional Scaling (MDS)plot 
MDSplot = function(betas, col.var = NULL) {
   
   d <- dist(betas) 
   fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
   x <- fit$points[,1]
   y <- fit$points[,2]
   
   if(is.null(col.var)) {
   par(mar = c(5,5,4,2))
   plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", 
        cex.axis = 1.5, cex.lab = 1.8, pch = 21, bg = "red", cex = 1.5)
   }
   else {
   	if (length(levels(as.factor(col.var))) >= 3) {	
   	cols = brewer.pal(length(levels(as.factor(col.var))), "Set3")
   	cols1 = cols	
   	names(cols) <- levels(as.factor(col.var))
   	}
   else {
   	cols = c("deepskyblue1", "darkorange1")
   	cols1 = cols
   	names(cols) <- levels(as.factor(col.var))
    }
   par(mar = c(5,5,4,9))
   plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS",
        cex.axis = 1.5, cex.lab = 1.8, bg = cols[as.factor(col.var)], pch = 21, cex = 1.5)
    text(x,y,rownames(betas))
    par(xpd = -2)
    legend("right", inset = c(-0.4,0), legend = levels(as.factor(col.var)), pch = 21, pt.bg = cols1, bty = "n", cex = 1.5)
  }
 }
 
# density plots of beta and M-values 

densplot = function(betas, col.var = NULL) {
	par(mar = c(5,5,4,2))
	
	if(is.null(col.var)) col.samp = rep("black", dim(betas)[1])
	else if(length(levels(as.factor(col.var))) >= 3) {
	    cols = brewer.pal(length(levels(as.factor(col.var))), "PuBuGn")
   	  cols1 = cols	
   	  names(cols) <- levels(as.factor(col.var))
   	  col.samp = cols[as.factor(col.var)]
   	}
    else if (length(levels(as.factor(col.var))) == 2) {
        cols = c("deepskyblue1", "darkorange1")
   	   cols1 = cols
   	   names(cols) <- levels(as.factor(col.var))
   	   col.samp = cols[as.factor(col.var)]
   }

	#density plot of the beta-values
	plot(x = c(0,0), y = c(0,0), xlim = c(0,1), ylim = c(0,5), cex.axis = 1.5, cex.lab = 1.8, xlab = expression(paste(beta, "-value")), 
	     col = "white", ylab = "Density")
	for(i in 1:dim(betas)[1]) {
        ymiss = is.na(betas[i,])
		lines(density(betas[i,!ymiss]), lwd = 1.5, col = col.samp[i])
	}
	
	if(!is.null(col.var)) {
		legend("topright", legend = levels(as.factor(col.var)), lwd = 1.5, col = cols1, bty = "n")
	}
}

# bar-plot depicting the percent of detection p-values greater than some threshold
dectPplot = function(detPvals, thr = 0.01) {
    percent = apply(detPvals, 2, function(w) mean(w >= thr))*100
    par(mar = c(5,5,4,2))
    barplot(percent, ylim = c(0, max(percent)), ylab = paste("DetP >", thr, "(%)"),
    cex.lab = 1.8, cex.axis = 1.5, xlab = "", xaxt = "n")
    ind = seq(0.7, by = 1.2, length.out = length(percent))
    axis(1, at = ind, labels = names(percent), las = 2)
}

# other relevant functions
logit2 = function(x) log2(x) - log2(1-x)
expit2 = function(x) 2^x/(1+2^x)



