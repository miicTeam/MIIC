#!/usr/bin/Rscript
## 2>&1 | tee logFiles/logFile.txt

#### Delete all variables
rm( list = ls ( all = TRUE ) )
options(warn=-1)
full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
                       error=function(e) # works when using R CMD
                         normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
rootDir = sub("gmPlot.R", "", full.fpath, fixed = TRUE)
require(getopt, quietly = TRUE)
require(igraph, quietly = TRUE, warn.conflicts = FALSE)
require(plotrix, quietly = TRUE)

cat( "\t# --------\n# -> START Plot...\n" )
#### Parse the command line
#### ----
##### Load libraries

# If testing
myTest = FALSE

#### FOR TEST
if( myTest == TRUE ){ setwd( "~/Projects/Projects_inference_VSign&Color/pipeline/inference/common" ) }
#### FOR TEST

#### Load supplementary functions
source( file.path( "..", "sharedLib", "shared.utils.lib.R" ) )
source( file.path( "lib", "gmPlot.lib.R" ) )

### FOR TEST ####
if( myTest == TRUE )
{
    plotArg.list.external = list()
    plotArg.list.external[["argOutDir"]] = "~/Project_miic/data/HEMATO/gottgens_2015/output/bypop_tfmq40_Binary_NewOrt/miic/test_colorSign/4FSGdata_TFandMQ_40_binary.tsv"
    plotArg.list.external[["argInData"]] = "~/Project_miic/data/HEMATO/gottgens_2015/input/unprocessed/data_by_pop_TFandMQ_40/4FSGdata_TFandMQ_40_binary.tsv"
    plotArg.list.external[["argLayout"]] = NULL
    plotArg.list.external[["argSummary"]] = "edgesList.miic.orientProba.summary.txt"
    plotArg.list.external[["argDoCurve"]] = FALSE
    plotArg.list.external[["argEntropy"]] = FALSE
    plotArg.list.external[["argLatent"]] = TRUE
}
### FOR TEST ####

#### Parse command line arguments
#### ----
plotArg.list = list();
if( !exists( "plotArg.list.external" ) )
{
	plotArg.list = plot.parseCommandLine()    
} else {
	
	plotArg.list = plotArg.list.external
}

#### Set Global variables
isCurved = plotArg.list[["argDoCurve"]]
myVariables = colnames( read.table(plotArg.list[["argInData"]],sep='\t',as.is=T,header=T, check.names=FALSE) )

setwd(plotArg.list[["argOutDir"]])

#### Define a log file
#### ----
logFilePath = "log.plot.txt"
if( file.exists( logFilePath ) ){ file.remove( logFilePath ) }
sink( file=logFilePath, append = TRUE, split = TRUE )

currdir = getwd()

setwd( rootDir )

#### Read the first line of the inputdata file to get all the properties saved in the correct order
gV <- new.env()
gV$allProperties =  unlist( strsplit( readLines( plotArg.list[["argInData"]], n = 1 ),split="\t" ) )
#### ----


#### Load the vertices positions

gV$myVerticesPos = NULL
if( length( plotArg.list[["argLayout"]] ) > 0 )
{
    if( file.exists( plotArg.list[["argLayout"]] ) )
    {

      
      
        gV$myVerticesPos <- as.data.frame( read.table( file = plotArg.list[["argLayout"]], header = FALSE, sep = '\t',
                                           stringsAsFactors = FALSE , check.names=FALSE) )
        if( ncol(gV$myVerticesPos) > 2 ) { gV$myVerticesPos <- gV$myVerticesPos[, -1] }
        gV$myVerticesPos <- as.matrix( gV$myVerticesPos )
    } else { warning("\tThe layout file does not exist") }
}
setwd(currdir)

#### ----
setwd(plotArg.list[["argOutDir"]])

#### Check the existence of the summary
if( file.exists( plotArg.list[["argSummary"]] ) == TRUE ){
  myPdfFile = sub( ".txt", ".plot_pCor.pdf", plotArg.list[["argSummary"]] )
  #plots for all technics on partial correlation
  mySummary = plot.loadSummary(plotArg.list[["argSummary"]])
  # create the color vector based on partial correlation values
  myLayout = layout_with_kk
  if( length( plotArg.list[["argLayout"]] ) > 0 ){
    myLayout = gV$myVerticesPos
  }
  else if (length(myVariables) >=40){
    myLayout = layout.circle

  }
  blue.gradient = rainbow(100, start = 3/6, end=4/6)
  red.gradient = rainbow(100, start=0, end=0.16)
  leg_colors = c(red.gradient,blue.gradient)


  if( length( na.omit( mySummary[,"partial_correlation"] ) ) <= nrow(mySummary) & length( na.omit( mySummary[,"partial_correlation"] ) ) > 0 ){
    myColors = pCor.edgeCol( mySummary, myVariables )
    # create the graph object
    myGraph = modif.Graph( mySummary, myVariables, myColors )
  
    # plot the Partial Correlation Graph
    pdf(myPdfFile)
    layout(t(1:2), widths=c(5,1))
    # Set margins and turn all axis labels horizontally (with `las=1`)
    par(mar=rep(.5, 4), oma=c(3,3,3,1), las=1)

    plot(myGraph, layout = myLayout )
    title( paste( "[", basename( plotArg.list[["argOutDir"]] ), "]" ), cex.main = 1 )  
    
    legend_image <- as.raster(matrix(leg_colors, ncol=1))
    plot(c(0,5),c(-1,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Partial correlation', cex.main = 0.7)
    text(x=1.5, y = seq(-1,1,l=5), labels = seq(-1,1,l=5), cex = 0.7)
    rasterImage(legend_image, 2.5, -1, 3.5,1)

    graphics.off()
  } else {
    leg_colors = red.gradient
  }

  # plot the Confidence Graph  only for miic on log_confidence
  if( length( na.omit( mySummary[["log_confidence"]] ) ) == nrow(mySummary) )
  { 
    myPdfFile = sub( ".txt", ".plot_confidence.pdf", plotArg.list[["argSummary"]] )
    myColors = conf.edgeCol( mySummary, myVariables )
    myGraph = modif.Graph( mySummary, myVariables, myColors )
    pdf(myPdfFile)
    layout(t(1:3), widths=c(5,1,1))
    #layout(t(1:2), widths=c(5,2))

     # Set margins and turn all axis labels horizontally (with `las=1`)
    par(mar=c(.5, .1,.5,.1), oma=c(3,3,3,1), las=1)
    plot(myGraph, layout = myLayout )
    title( paste( "[", basename( plotArg.list[["argOutDir"]] ), "]" ), cex.main = 1 )

    # Legend
    # Positive correlations
    legend_image <- as.raster(matrix(red.gradient, ncol=1))
    plot(c(0, 4),c(0.2,0.8),type = 'n', axes = F,xlab = '', ylab = '')
    par(adj=1)
    if( length( na.omit( mySummary[,"sign"] ) ) == nrow(mySummary) ){
      title('Confidence\npcor+', cex.main = 1, line=-4)
    

      rasterImage(legend_image, 3.3, 0.25, 3.8, 0.75)
      # Negative correlations
      legend_image <- as.raster(matrix(rev(blue.gradient), ncol=1))
      plot(c(0, 4),c(0.2,0.8),type = 'n', axes = F,xlab = '', ylab = '')
      par(adj=0.5)
      title("(NI' = -log Pxy)\npcor-", cex.main = 1, line=-4)
      text(x=rep(.5,3), y = c(0.26,0.5,0.74) , labels = c("< 1", "50", "> 100"), cex = 1)
      rasterImage(legend_image, 1.5, 0.25, 2, 0.75)

    } else {
      # legend_image <- as.raster(matrix(leg_colors, ncol=1))
      # plot(c(0,5),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Confidence', cex.main = 0.7)
      # text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5), cex = 0.7)
      # rasterImage(legend_image, 3.3, 0, 3.8, 1)

      graphics::title("Confidence   \n(NI' = -log Pxy)", cex.main = 1, line=-3)
      graphics::rasterImage(legend_image, 2.5, 0.25, 3, 0.75)
      graphics::plot(c(0, 4),c(0.2,0.8),type = 'n', axes = F,xlab = '', ylab = '')
      graphics::par(adj=0.3)
      graphics::text(x=rep(.5,3), y = c(0.26,0.5,0.74), labels = c("< 1", "50", "> 100") , cex = 1)
    }

    graphics.off()
  }

  # Terminate log
  sink(file=NULL)
  
  cat( "\t# --------\n# -> END Plot...\n" )
}
