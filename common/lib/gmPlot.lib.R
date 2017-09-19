plot.parseCommandLine <- function()
{
    #### Specify the input arguments
    mySpec = matrix( c( 'inputData', 'i', 1, "character"	# File path to the input data
                  , 'outDirPath', 'o', 1, "character"	    # Path to the working directory (skeleton output)	
                  , 'layout', 'l', 1, "character"           # File path to the vertices layout
                  , 'summary', 's', 1, "character"          # File path to the graph summary
                  , 'isCurved', 'c', 0, "logical"           # Curve edges
                  , 'isLatent', 'h', 0, "logical" )          # Are the hidden variables allowed
                  , byrow = TRUE, ncol = 4 );

    #### Read the input arguments
    opt = getopt(mySpec)
    
    #### Variable for the arguments
    myInputData = character(0)
    myOutDirPath = character(0)
    myLayout = character(0)
    mySummary = character(0)
    myLatent = FALSE
    myDoCurve = FALSE
    myIsVerbose = FALSE
    
    #### Check the input arguments
    if( !is.null( opt$inputData ) )
    { myInputData <- opt$inputData } else { stop("The input data file is required (-i)") }
    
    if( !is.null( opt$outDirPath ) )
    { myOutDirPath <- opt$outDirPath } else { stop("The output dir path is required (-o)") }
    
    if( !is.null( opt$layout ) )
    { myLayout <- opt$layout } else { cat("\t# Reminder -> No file path for the nodes layout given (-l)\n") }
    
    if( !is.null( opt$summary ) )
    { mySummary <- opt$summary } else { stop("The learnt graphical model summary is required (-s)") }
    
    if( !is.null(opt$isCurved) ) {myDoCurve <- TRUE }
    if( !is.null(opt$isLatent) ) {myLatent <- TRUE }

    if( !is.null(opt$isVerbose) ) {myIsVerbose <- TRUE }

    return( list( argInData = myInputData, argOutDir = myOutDirPath
                    , argLayout = myLayout, argSummary = mySummary
                    , argDoCurve = myDoCurve, argVerbose=myIsVerbose, argLatent = myLatent  ))
}



plot.loadSummary <- function( mySummaryFilePath )
{
    #### Load the summary of the edges
    mySummary <- read.table( file = mySummaryFilePath, header = TRUE, row.names = 1, sep = '\t', stringsAsFactors = FALSE , check.names = F)
    rownames(mySummary) = c()

    #### Ignore the TN edges
    myTypesToIgnore = c()
    myTypesToIgnore = c("TN","N","FN")
    myLinesToIgnore = which( mySummary[,"type"] %in% myTypesToIgnore )
    if( length( myLinesToIgnore ) > 0 ) { mySummary = mySummary[-myLinesToIgnore,]} 
    return(mySummary) 
}

plot.createDefaultGraph <- function( mySummary, myAllGenes )
{
    #### Replace names by numbers to create the graph
    inf.edgesList.nbr <- apply( mySummary[, c("x","y")], MARGIN = c(1,2), function(x) { x <- which( myAllGenes == x ) } )

    #### Create an unoriented igraph with all the nodes
    inf.graph <- graph( t( inf.edgesList.nbr ), length( myAllGenes ), directed = TRUE )

    #### Set the vertices options
    V(inf.graph)$label <- myAllGenes
    V(inf.graph)$shape <- "circle"
    V(inf.graph)$color <- "lightblue"
    V(inf.graph)$label.family <- "Helvetica"
    V(inf.graph)$label.cex <- 0.6
    V(inf.graph)$size <- 10

    #### Set the general edges options
    E(inf.graph)$arrow.size <- 0.5
    E(inf.graph)$arrow.width <- 3
    E(inf.graph)$width <- 3
    E(inf.graph)$curved <- FALSE
    E(inf.graph)$color <- "red2"
    E(inf.graph)$lty <- "solid"
    E(inf.graph)$arrow.mode <- 0
    return(inf.graph)
}

plot.setOrientation <- function( mySummary, inf.graph )
{
    #### Set the options specific for forward oriented
    ort.fwd.idx <- which( ( mySummary[, "infOrt"] %in% c(2,4) ) | ( mySummary[, "type"] == 'FN' & mySummary[, "trueOrt"] == 2 ) )
    if( length( ort.fwd.idx ) > 0 ){ E(inf.graph)[ort.fwd.idx]$arrow.mode <- 2 }

    #### Set the options specific for backward oriented
    ort.bck.idx <- which( ( mySummary[, "infOrt"] %in% c(-2,-4) ) | ( mySummary[, "type"] == 'FN' & mySummary[, "trueOrt"] == (-2) ) )
    if( length( ort.bck.idx ) > 0 ) { E(inf.graph)[ort.bck.idx]$arrow.mode <- 1 }
    
    #### Set the options specific for bidirectional orientations
    bidir.idx <- which( mySummary[, "infOrt"] == 6 )
    if( length( bidir.idx ) > 0 )
    {
        E(inf.graph)[bidir.idx]$arrow.mode <- 3
    }
    return(inf.graph)

}

littlefunc <- function(edge1,edge2)
{
    return(paste(edge1, collapse=",") == paste(rev(edge2), collapse=","))
}

# ---- Function to plot graphes with edges matching to their partial correlation
pCor.edgeCol <- function(summary, features)
{
    # Define the color gradients
    blue.gradient = rainbow(100, start = 3/6, end=4/6)
    red.gradient = rev(rainbow(100, start=0, end=0.16))

    myEdgesColor = rep(NA, nrow(summary)) # Set the color vector for the edges
    max.pcor.neg = min(summary[which(summary[,"sign"] == "-"),"partial_correlation"]) # get the maximum negative pcor
    max.pcor.pos = max(summary[which(summary[,"sign"] == "+"),"partial_correlation"]) # get the maximum positive pcor
    for(edge in 1:nrow(summary)) # loop on all the edges present in the network
    {   
        if(! is.na(summary[edge, "sign"]) )     
        {
        # Set the correct tmp.max.pcor
        if(sign(summary[edge, "partial_correlation"]) == -1){tmp.max.pcor = max.pcor.neg}
        else {tmp.max.pcor = max.pcor.pos}
        # Compute the ratio between the tmp.max.pcor and the edge's pcor, and use it as an index to get a color
        edge.pCor.ind = abs(summary[edge, "partial_correlation"])/abs(tmp.max.pcor)
        edge.colIndex = round(edge.pCor.ind * 100)
        if( edge.colIndex == 0 ) { edge.colIndex = 1}
        ### Get the sign of the link to look at the correct color gradient

            if(summary[edge, "sign"] == "+")
            {   
                myEdgesColor[edge] = red.gradient[edge.colIndex]
            }
            else
            { 
                myEdgesColor[edge] = blue.gradient[edge.colIndex]
            }
        }
        else { myEdgesColor[edge] = "grey48" }
        
    }

    return(myEdgesColor)
}

# ---- Function to plot graphes with edges matching to their mutual information (log_confidence column in summary)
conf.edgeCol <- function(summary, features)
{
    # Define the color gradients
    blue.gradient = rainbow(100, start = 3/6, end=4/6)
    red.gradient =  rev(rainbow(100, start=0, end=0.16))
    myEdgesColor = rep(NA, nrow(summary)) # Set the color vector for the edges
    max.conf = 100 # get the maximum
    min.conf = 1 # get the minimum

    for(edge in 1:nrow(summary)) # loop on all the edges present in the network
    {
        # Set the correct tmp.max.pcor
        # Compute the ratio between the tmp.max.pcor and the edge's pcor, and use it as an index to get a color
        edge.colIndex = round(summary[edge, "log_confidence"])
        if( edge.colIndex < min.conf  ) { edge.colIndex = 1 }
        else if( edge.colIndex > max.conf ){ edge.colIndex = 100 }
        ### Get the sign of the link to look at the correct color gradient
        if(! is.na(summary[edge, "sign"]) )     
        {
            if(summary[edge, "sign"] == "+")
            {   
                myEdgesColor[edge] = red.gradient[edge.colIndex]
            }
            else
            { 
                myEdgesColor[edge] = blue.gradient[edge.colIndex]
            }
        }
        else { myEdgesColor[edge] = "grey48" }
        
    }

    return(myEdgesColor)
}

# ---- Function which returns a graph object from edges colors and node sizes eventually
modif.Graph <- function(summary, features, edgeColors, nodeSizes = 10, nodeColors = 'lightblue')
{
    mygraph = plot.createDefaultGraph(summary, features)
    E(mygraph)$color = edgeColors
    V(mygraph)$color = nodeColors
    mygraph = plot.setOrientation(summary, mygraph)
    V(mygraph)$size = nodeSizes
    return(mygraph)
}
