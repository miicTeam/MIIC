#!/usr/bin/Rscript

computeSHD2SampleStat <- function( myTrueGraphFilePath, myFinalSummaryFilePath, skeletonStats, myInfMethod, myAllProperties, myVb = FALSE, plotDir = NULL )
{
	#### Scores to return	
	myErrWarKeys = c( "infCpdagErrWar", "infVstructErrWar", "infSkelErrWar", "trueCpdagErrWar", "trueVstructErrWar", "trueSkelErrWar" )
	myKeys = c( "distance__pattern", "distance__cpdag", "distance__inferred", myErrWarKeys )
	myRetList = vector( "list", length(myKeys) )
	names( myRetList ) = myKeys
	myRetList = lapply(myRetList, function(x) x = NA)
	#### Prepare a list to store the graphs to compare
	myKeys = c( "true", "inferred" )
	graphToCompare = vector( "list", length(myKeys) )
	names( graphToCompare ) = myKeys
	myKeys = c( "pattern", "cpdag", "graph", "vstructs", "skeleton", "dataFilePath" )
	graphToCompare = lapply( graphToCompare, function( x ) { tmpList  = vector( "list", length(myKeys) ); names(tmpList) = myKeys; x = tmpList } )

	graphToCompare[["true"]][["dataFilePath"]] = myTrueGraphFilePath
	graphToCompare[["inferred"]][["dataFilePath"]] = myFinalSummaryFilePath

	#### TRUE ---------------------------------------------------
	# - Load the true model either from a .bif file or an edge list
	if( grepl( ".bif", myTrueGraphFilePath ) == TRUE )
	{
	    trueGraph.bif = read.bif( file = graphToCompare[["true"]][["dataFilePath"]], debug = myVb )
	    myAllProperties = names(trueGraph.bif)
	    # - Make a model string of the graph from the bn.fit object
	    trueGraph.modelString = modelstring( x = trueGraph.bif )
	    # - Make a bn object from the model string
	    graphToCompare[["true"]][["graph"]] = model2network( string = trueGraph.modelString, ordering = myAllProperties, debug = myVb )

	} else {
	    # - Read the edge list
	    trueGraph.edgelist = read.table( file = graphToCompare[["true"]][["dataFilePath"]], as.is = TRUE, header = FALSE, sep = '\t', check.names=F )
	    trueGraph.edgelist = trueGraph.edgelist[,c(1,2)]    # the third column, if any, is for the essentiality
	    #### Define an empty true graph and add the edges from the list
    	graphToCompare[["true"]][["graph"]] = bnlearn::empty.graph( myAllProperties )
    	
    	# Set the edges and their direction (we suppose here that all edges are oriented as col1->col2)
    	for( iTrueEdge in seq_len(nrow(trueGraph.edgelist)) )
	    {
		
			myX = as.character(trueGraph.edgelist[iTrueEdge, 1])
			myY = as.character(trueGraph.edgelist[iTrueEdge, 2])
            graphToCompare[["true"]][["graph"]] = bnlearn::set.arc( x = graphToCompare[["true"]][["graph"]], from = myX, to= myY, check.cycles = FALSE, debug = myVb)
    	}
	}
	
	#### Plot the true graph
	if( !is.null(plotDir) )
	{
		graphviz.plot( x = graphToCompare[["true"]][["graph"]], main = paste( "True Graph", basename(graphToCompare[["true"]][["dataFilePath"]]), sep=", " ) )
		dev.copy2pdf( file = file.path( plotDir, "trueGraph.pdf" ) )
		dev.off()
	}

	#### INFERRED ----------------------------------------------
	#### Load the summary of the edges
	infEdges.summary <- read.table( file = graphToCompare[["inferred"]][["dataFilePath"]], header = TRUE, row.names = 1, sep = '\t', stringsAsFactors = FALSE , check.names=F )
	infEdges.summary.save <- infEdges.summary
	#### Ignore the FN and TN edges
	myLinesToIgnore <- which( infEdges.summary[,"type"] == 'TN' | infEdges.summary[,"type"] == 'FN' )
	if( length( myLinesToIgnore ) > 0 ) { infEdges.summary = infEdges.summary[-myLinesToIgnore,]} 

	#### Define first an empty graph for the inferred one
	graphToCompare[["inferred"]][["graph"]] = empty.graph( myAllProperties )
	
	# Set the edges and their direction if any
	for( iInfEdge in seq_len(nrow(infEdges.summary)) )
	{
		#### Get the direction
		# - if '2', then x --> y, - if '-2', then y --> x, - if '1', set an undirected edge
		myDirection = infEdges.summary[iInfEdge, "infOrt"]
		
		if( abs( myDirection ) > 0 )
		{
			myX = as.character(infEdges.summary[iInfEdge, "x"])
			myY = as.character(infEdges.summary[iInfEdge, "y"])

			if( myDirection == 2 )
			{
				graphToCompare[["inferred"]][["graph"]] = bnlearn::set.arc( x = graphToCompare[["inferred"]][["graph"]], from = myX, to= myY, check.cycles = FALSE, debug = myVb)
			
			} else if( myDirection == (-2) )	{

				graphToCompare[["inferred"]][["graph"]] = bnlearn::set.arc( x = graphToCompare[["inferred"]][["graph"]], from = myY, to= myX, check.cycles = FALSE, debug = myVb)
			
			} else if( myDirection == 1 ) {
			#} else if( myDirection %in% c(1,6) ) {
	
				graphToCompare[["inferred"]][["graph"]] = bnlearn::set.edge( x = graphToCompare[["inferred"]][["graph"]], from = myX, to= myY, check.cycles = FALSE, debug = myVb)

			}
		} 
	}

	#### Plot the inferred graph
	if( !is.null(plotDir) )
	{
		graphviz.plot( x = graphToCompare[["inferred"]][["graph"]], main = paste( "Inferred Graph", basename(graphToCompare[["inferred"]][["dataFilePath"]]), sep=", " ) )
		dev.copy2pdf( file = file.path( plotDir, "inferredGraph.pdf" ) )
		dev.off()
	}

	#### Get CPDAG and PATTERN graph
	#### ------------------------------------------------------------
	for( strType in names( graphToCompare ) )
	{
		# - the cpdag
		# ------------------------------------------------------------
		myKey = ifelse( strType == "true", "trueCpdagErrWar", "infCpdagErrWar" )
		myTry = tryCatch( bnlearn::cpdag( graphToCompare[[strType]][["graph"]], moral = TRUE,  debug = myVb ), error = function( e ) { e }, 	warning = function( w ) { w } )

		if( !( inherits( myTry, "simpleError" ) | inherits( myTry, "simpleWarning" ) )  )
		{ 
			# if no error or if warning, get the result
			graphToCompare[[strType]][["cpdag"]] = myTry

		} else if( inherits( myTry, "simpleError" ) ) {
			
			# if there is an error
			graphToCompare[[strType]][["cpdag"]] = NULL
			myRetList[[myKey]] = gsub( "\n", "", as.character(myTry) )

		} else if( inherits( myTry, "simpleWarning" ) ) {

			# if it is only a warning
			myRetList[[myKey]] = gsub( "\n", "", as.character(myTry) )
			graphToCompare[[strType]][["cpdag"]] = bnlearn::cpdag( graphToCompare[[strType]][["graph"]], moral = TRUE,  debug = myVb )
			assign("last.warning", NULL, envir = baseenv())
		}

		# - the skeleton
		# ------------------------------------------------------------
		myKey = ifelse( strType == "true", "trueSkelErrWar", "infSkelErrWar" )
		myTry = tryCatch( bnlearn::skeleton( graphToCompare[[strType]][["graph"]] ), 
											error = function( e ) { e }, 	warning = function( w ) { w } )

		if( !( inherits( myTry, "simpleError" ) | inherits( myTry, "simpleWarning" ) )  )
		{ 
			# if no error or if warning, get the result
			graphToCompare[[strType]][["skeleton"]] = myTry

		} else if( inherits( myTry, "simpleError" ) ) {

			# if there is an error
			graphToCompare[[strType]][["skeleton"]] = NULL
			myRetList[[myKey]] = gsub( "\n", "", as.character(myTry) )

		} else if( inherits( myTry, "simpleWarning" ) ) {

			# if it is only a warning
			myRetList[[myKey]] = gsub( "\n", "", as.character(myTry) )
			graphToCompare[[strType]][["skeleton"]] = bnlearn::skeleton( graphToCompare[[strType]][["graph"]] )
			assign("last.warning", NULL, envir = baseenv())
		}

		# - the vStructures
		myKey = ifelse( strType == "true", "trueVstructErrWar", "infVstructErrWar" )
		myTry = tryCatch( bnlearn::vstructs(graphToCompare[[strType]][["graph"]], moral = FALSE, debug = myVb), 
											error = function( e ) { e }, 	warning = function( w ) { w } )

		if( !( inherits( myTry, "simpleError" ) | inherits( myTry, "simpleWarning" ) )  )
		{ 
			# if no error or if warning, get the result
			graphToCompare[[strType]][["vstructs"]] = myTry

		} else if( inherits( myTry, "simpleError" ) ) {

			# if there is an error
			graphToCompare[[strType]][["vstructs"]] = NULL
			myRetList[[myKey]] = gsub( "\n", "", as.character(myTry) )

		} else if( inherits( myTry, "simpleWarning" ) ) {

			# if it is only a warning
			myRetList[[myKey]] = gsub( "\n", "", as.character(myTry) )
			graphToCompare[[strType]][["vstructs"]] = bnlearn::vstructs(graphToCompare[[strType]][["graph"]], moral = FALSE, debug = myVb)
			assign("last.warning", NULL, envir = baseenv())
		}

		# if(!is.null(graphToCompare[[strType]][["cpdag"]]) )
		# {
		# 	cat("STRTYPE",strType," CLASS OF MY CPDAG!!!! ::::", class(graphToCompare[[strType]][["cpdag"]]), "\n")
		# 	pdf("/home/louis/miic_memoryAndThreads/common/testCPDAG.pdf")
		# 	graphviz.plot( x = graphToCompare[[strType]][["cpdag"]], main = paste( strType, "CPDAG", basename(graphToCompare[[strType]][["dataFilePath"]]), sep=", "))		
			
		# 	graphics.off()
		# }

		#### Create the pattern --> skeleton + vstructs
		if( !is.null(graphToCompare[[strType]][["skeleton"]]) )
		{
			# - Go through the table of vstructures and set the orientation X --> Z and Y --> Z on the skeleton
			graphToCompare[[strType]][["pattern"]] = graphToCompare[[strType]][["skeleton"]]
			if( !is.null( graphToCompare[[strType]][["vstructs"]] ) )
			{
				for( iVstruct in seq_len(nrow(graphToCompare[[strType]][["vstructs"]])) )
				{
					# Sink of the vStructure
					strTo = graphToCompare[[strType]][["vstructs"]][iVstruct, "Z"]
	
					# X --> Z & Y --> Z
					for( nBase in c("X", "Y") )
					{
						strFrom = graphToCompare[[strType]][["vstructs"]][iVstruct, nBase]
						graphToCompare[[strType]][["pattern"]] = set.arc( x = graphToCompare[[strType]][["pattern"]], from = strFrom, to= strTo, check.cycles = FALSE, debug = myVb)	
					}
				}
			}

			if( !is.null(plotDir) )
			{
				graphviz.plot( x = graphToCompare[[strType]][["pattern"]], main = paste( strType, "CPDAG", basename(graphToCompare[[strType]][["dataFilePath"]]), sep=", "))		
				dev.copy2pdf( file = file.path( plotDir, paste( strType, "Graph.pattern.pdf", sep="" ) ) )
				dev.off()
			}
		}
	}


	if( !is.null(graphToCompare[["true"]][["cpdag"]]) & !is.null(graphToCompare[["inferred"]][["graph"]]) )
	{
		infEdges.summary.save = computeSHD2( myTrueGraph = graphToCompare[["true"]][["cpdag"]], myInfGraph = graphToCompare[["inferred"]][["graph"]], infEdges.summary.save )

	} else { cat( "# -----> CPDAG SHD2 = NA\n") 
		TP_skel = length(which(infEdges.summary.save[,"type"] == "TP"))
		FP_skel = length(which(infEdges.summary.save[,"type"] == "FP"))
		FN_skel = length(which(infEdges.summary.save[,"type"] == "FN"))

		Prec_skel = TP_skel / (FP_skel + TP_skel)
		Rec_skel = TP_skel / (FN_skel + TP_skel)
		Fscore_skel = 2 * Prec_skel * Rec_skel / ( Prec_skel + Rec_skel)
		# Print the results
		line="\tPrecision\tRecall\tFscore\n"
		line=paste(line,"Skeleton\t",Prec_skel,"\t",Rec_skel,"\t",Fscore_skel,"\n",sep="")
	}


	return(infEdges.summary.save)
}

computeSHD2 <- function( myTrueGraph, myInfGraph, infEdges.summary.save, myVb = FALSE )
{

	#### Compute the distance between the two graph (shd2)
	#### ------------------------------------------------------------
	# - For each true oriented edge
	# -- if the edge exists in inferred and is oriented in opposite direction or not oriented --> +1
	# - For each true non-oriented edge
	# -- if the edge exists in inferred and is oriented --> +1
	shd2 = 0

	if( myVb == TRUE  ){ cat( "# Oriented true edges:\n# ----\n" ) }
	for(iOriented in seq_len(nrow(directed.arcs(myTrueGraph))))
	{
		# Get the oriented edge
		myEdge = directed.arcs(myTrueGraph)[iOriented, ]

		# if the edge has not been inferred (FN) or has been inferred with the true direction --> OK
		# else shd2 +1
		if( !( ( amat(myInfGraph)[myEdge["from"], myEdge["to"]] == 0 & amat(myInfGraph)[myEdge["to"], myEdge["from"]] == 0 ) | 
				( amat(myInfGraph)[myEdge["from"], myEdge["to"]] == 1 & amat(myInfGraph)[myEdge["to"], myEdge["from"]] == 0 ) ) )
		{
			shd2 = (shd2 + 1)
			infEdges.summary.save[which((infEdges.summary.save[,"x"] == myEdge["from"] & infEdges.summary.save[,"y"] == myEdge["to"]) | (infEdges.summary.save[,"x"] == myEdge["to"] & infEdges.summary.save[,"y"] == myEdge["from"]) ), "isOrtOk"] ="N"

		}
		else
		{
			infEdges.summary.save[which(infEdges.summary.save[,"x"] == myEdge["from"] & infEdges.summary.save[,"y"] == myEdge["to"]| (infEdges.summary.save[,"x"] == myEdge["to"] & infEdges.summary.save[,"y"] == myEdge["from"]) ), "isOrtOk"] ="Y"
		}
	}

	# The non oriented edge are given twice in the table
	if( myVb == TRUE ){ cat( "\n# Non-oriented true edges:\n# ----\n" ) }
	# - use a vector to indicate which edge has already been processed
	isDone = rep( 0, nrow(undirected.arcs(myTrueGraph)) )
	for(iNotOriented in seq_len(nrow(undirected.arcs(myTrueGraph))))
	{
		if( isDone[iNotOriented] == 1 ) { next }

		# Get the oriented edge
		myEdge = undirected.arcs(myTrueGraph)[iNotOriented, ]

		# Set the same edge given in the 'opposite direction' as done
		myIdenticalEdge = which( undirected.arcs(myTrueGraph)[, "from"] == myEdge["to"] &
															undirected.arcs(myTrueGraph)[, "to"] == myEdge["from"] )
		if( length(myIdenticalEdge) == 1 ) { isDone[myIdenticalEdge] = 1}
		else { stop( "# Only one non Oriented (", paste(myEdge, collapse = "--"), ") in the table!\n", sep = "" )  }

		# if the edge has not been inferred (FN) or has been inferred with no direction --> OK
		# else shd2 +1
		if( !( ( amat(myInfGraph)[myEdge["from"], myEdge["to"]] == 0 & amat(myInfGraph)[myEdge["to"], myEdge["from"]] == 0 ) | 
				( amat(myInfGraph)[myEdge["from"], myEdge["to"]] == 1 & amat(myInfGraph)[myEdge["to"], myEdge["from"]] == 1 ) ) )
		{
			shd2 = (shd2 + 1)
			infEdges.summary.save[which(infEdges.summary.save[,"x"] == myEdge["from"] & infEdges.summary.save[,"y"] == myEdge["to"]| (infEdges.summary.save[,"x"] == myEdge["to"] & infEdges.summary.save[,"y"] == myEdge["from"]) ), "isOrtOk"] ="N"

		}
		else
		{
			infEdges.summary.save[which(infEdges.summary.save[,"x"] == myEdge["from"] & infEdges.summary.save[,"y"] == myEdge["to"]| (infEdges.summary.save[,"x"] == myEdge["to"] & infEdges.summary.save[,"y"] == myEdge["from"]) ), "isOrtOk"] ="Y"			
		}
	}

	# Compute the various scores for the comparison
	TP_skel = length(which(infEdges.summary.save[,"type"] == "TP"))
	FP_skel = length(which(infEdges.summary.save[,"type"] == "FP"))
	FN_skel = length(which(infEdges.summary.save[,"type"] == "FN"))
	TP_ort = TP_skel - shd2
	FP_ort = FP_skel + shd2

	cat("SHD2 :::: ", shd2,"\n")
	cat("TP_skel :::: ", TP_skel,"\n")
	cat("FP_skel :::: ", FP_skel,"\n")
	cat("TP_ort :::: ", TP_ort,"\n")
	cat("FP_ort :::: ", FP_ort,"\n")

	Prec_skel = TP_skel / (FP_skel + TP_skel)
	Rec_skel = TP_skel / (FN_skel + TP_skel)
	Fscore_skel = 2 * Prec_skel * Rec_skel / ( Prec_skel + Rec_skel)

	Prec_ort = TP_ort / (FP_ort + TP_ort)
	Rec_ort = TP_ort / (FN_skel + TP_ort)
	Fscore_ort = 2 * Prec_ort * Rec_ort / ( Prec_ort + Rec_ort)

	# Print the results
	line="\tPrecision\tRecall\tFscore\n"
	line=paste(line,"Skeleton\t",Prec_skel,"\t",Rec_skel,"\t",Fscore_skel,"\n",sep="")
	line=paste(line,"CPDAG\t",Prec_ort,"\t",Rec_ort,"\t",Fscore_ort,"\n",sep="")
    writeLines(line, "stats.txt")
    
	return(infEdges.summary.save)

}

