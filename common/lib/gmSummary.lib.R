fromStringToNumberArrowType <- function(val){
    ret = 0
    if(val == "arrow")
        ret = 6
    else if(val == "T")
        ret = 15

    return(ret)
}
summary.parseCommandLine <- function()
{
    #### Specify the input arguments
    mySpec = matrix( c( 'inputData', 'd', 1, "character"	# File path to the input data
				  , 'outDirPath', 'o', 1, "character"        # Path to the working directory (skeleton output)
				  , 'stateOrderFile', 's', 1, "character"	# File path to the description of each variable state order
                      , 'trueEdges', 't', 1, "character"        # File path to the true egdes
                      , 'layoutFile', 'l', 1, "character"
                      , 'inferredEdges', 'i', 1, "character"    # File path to the inferred egdes
                      , 'adjMat', 'a', 1, "character"           # File path to the adjacency matrix
                      , 'isContinuous', 'c', 1, "integer"		# data continuous or not.
				              , 'writeNetwork', 'n', 0, "logical"
                      , 'isVerbose', 'v', 0, "logical")         # Verbose
                      , byrow = TRUE, ncol = 4 );

    #### Read the input arguments
    opt = getopt(mySpec)

    #### Variable for the arguments
    myInputDataFile = character(0)
    myOutDirPath = character(0)
    myStateOrderFile = character(0)
    myLayoutFile = character(0)
    myTrueEdgesFilePath = character(0)
    myInfEdgesFilePath = character(0)
    myAdjMatFilePath = character(0)
    myIsCnt = 0
    myIsVerbose = FALSE
    myWriteNetwork = FALSE

    #### Check the input arguments
    # ----
    # ----
    if( !is.null( opt$inputData ) )
    {myInputDataFile <- opt$inputData } else { stop("The input data file is required (-d)") }

    if( !is.null( opt$outDirPath ) )
    {myOutDirPath <- opt$outDirPath } else { stop("The output dir path is required (-o)") }

	# Variables states order file
    if( !is.null( opt$stateOrderFile ) ) { myStateOrderFile <- opt$stateOrderFile } else { cat("\t# Reminder -> No state order file provided (-s)\n") }
  	
  	# Layout file (used for graphml)
    if( !is.null( opt$layoutFile ) ) { myLayoutFile <- opt$layoutFile } else { cat("\t# Reminder -> No layout file provided (-l)\n") }

    if( !is.null( opt$inferredEdges ) )
    { myInfEdgesFilePath <- opt$inferredEdges } else { stop("The filename of the inferred edges is required (-e)") }    

    # True edges file
    if( !is.null( opt$trueEdges ) )
    { myTrueEdgesFilePath <- opt$trueEdges } else { cat("\t# Reminder -> No file path for the true edges given (-t)\n") }
    
    # Set if yser wants the graphml or xgmml
    if( !is.null( opt$writeNetwork ) ){ 
        myWriteNetwork <- TRUE
    } else { cat("\t# Reminder -> Network will not be written in graphml format (-n)\n") }


    if( !is.null( opt$adjMat ) )
    { myAdjMatFilePath <- opt$adjMat } else { stop("The filename of the adjacency matrix is required (-a)") }
       
    if( !is.null(opt$isVerbose) ) {myIsVerbose <- TRUE }

    if( !is.null(opt$isContinuous) ) {myIsCnt <- opt$isContinuous }

    return( list( argInData = myInputDataFile, argOutDir = myOutDirPath,
    			stateOrder = myStateOrderFile,
    			argWriteNetwork = myWriteNetwork,
     			argTrueEdgesFile = myTrueEdgesFilePath,
     			argLayoutFile = myLayoutFile, 
     			argInfEdgesFile = myInfEdgesFilePath,
     			argVerbose = myIsVerbose,
     			argCnt = myIsCnt,
    			argAdjMat = myAdjMatFilePath ))
}



### AJOUT LV #### 
### Functions to get the sign of the edges


findUnassignedNodes <- function(myGv, myStateOrderTable){
	unassignedNodes = c()

	if(!is.null(myStateOrderTable) ){
		if(ncol(myStateOrderTable) == 2)
			cols2look = colnames(myGv$data)
		else
			cols2look = myStateOrderTable[which(myStateOrderTable[,"var_type"] == 0),1]

		for(col in cols2look){
			mySafety = sort(as.character(unique(myGv$data[which(!is.na(myGv$data[,col])), col])))

			if(!is.numeric(myGv$data[,col]) & ncol(myStateOrderTable) == 2){
				# if stateorder is provided by user (it must have levels order), same code for 3 columns
				if("levels_increasing_order" %in% colnames(myStateOrderTable)){
					myCatStr = unlist(strsplit(myStateOrderTable[col,"levels_increasing_order"], ","))
					if(identical(mySafety, sort(myCatStr)) == FALSE) {
						unassignedNodes = c(unassignedNodes, col)
					}
					else{
						myGv$data[which(!is.na(myGv$data[,col])), col] =  as.integer(factor(myGv$data[which(!is.na(myGv$data[,col])),col], levels = myCatStr))
						### Convert the NA Values as well
						myGv$data[,col] = as.integer(myGv$data[,col])
					}
				} else {
					unassignedNodes = c(unassignedNodes, col)
				}				
			}
			else if(ncol(myStateOrderTable) == 3 ){
				myCatStr = unlist(strsplit(myStateOrderTable[col,"levels_increasing_order"], ","))
				print(myCatStr)
				if(identical(mySafety, sort(myCatStr)) == FALSE) {
					unassignedNodes = c(unassignedNodes, col)
				}
				else{
					myGv$data[which(!is.na(myGv$data[,col])), col] =  as.integer(factor(myGv$data[which(!is.na(myGv$data[,col])),col], levels = myCatStr))
					### Convert the NA Values as well
					myGv$data[,col] = as.integer(myGv$data[,col])
				}
				
			}


		}
	} else{
		cols2look = colnames(myGv$data)
		for(col in cols2look){	
			if(!is.numeric(myGv$data[,col]) ){
				unassignedNodes = c(unassignedNodes, col)	
				myGv$data[, col] <- factor(myGv$data[, col])
				myGv$data[, col] <- as.numeric(myGv$data[,col])
			}
		}
	}

	return( unassignedNodes )
}

# To compute sign with partial correlation, based on the previously determined ui set
computeSign.cont.pcor <- function(edgeList, myGv, outDirPath, originalPath)
{

	edgeList.signed <- cbind(edgeList, c(rep(NA,nrow(edgeList))), c(rep(NA,nrow(edgeList))))
	colnames(edgeList.signed) <- c(colnames(edgeList),"Sign", "Diff")
	tmp.core = NA
	# Find nodes for which we can't evaluate correlations (categorical variables without order reference)
	unassignedNodes = c()
	myStateOrderTable = NULL
	if( length( stateOrderFile ) > 0 )
	{

    	if( file.exists( stateOrderFile ) ){
    		myStateOrderTable = read.table(stateOrderFile, header=T, as.is=T, sep='\t')
    		rownames(myStateOrderTable) = myStateOrderTable[,"var_names"]
		}

    }			
	unassignedNodes = findUnassignedNodes(myGv, myStateOrderTable)
	dataTest = myGv$data+1

	write.table(dataTest, paste(outDirPath, "/datasetOrdered.txt",sep=""), sep="\t", col.names = T, row.names = F,quote = F)
	write(unassignedNodes, paste(outDirPath, "/unassignedNodes.txt",sep=""), sep="\t")
  
	setwd(originalPath)
	command = paste(paste(outDirPath, "/datasetOrdered.txt",sep=""), 
	                paste(outDirPath, "/unassignedNodes.txt",sep=""), 
	                paste(outDirPath, "/adjacencyMatrix.miic.orientProba.txt",sep=""),
	                paste(outDirPath, "/mediationAnalysis.txt",sep=""))

	# system2("../extendedAnalysis/extendedMediation", command)
  	setwd(outDirPath)
  
	edgeTabIdx = which(edgeList[,"confidence"] > 0)
	# Loop on all edges
	for(edge in edgeTabIdx)
	{

		edgeList[edge,"Sign"] <- NA; edgeList[edge,"Diff"] <- NA 			

		### get the Ui set of separation if it exists
		x <- myGv$data[,edgeList[edge,"x"]] ; y <- myGv$data[,edgeList[edge,"y"]]
		x.name <- edgeList[edge,"x"]; y.name <- edgeList[edge,"y"]
		# Check if nodes are in unassignedNodes
		if(!x.name %in% unassignedNodes & !y.name %in% unassignedNodes)
		{
			if(!is.na(edgeList[edge,"ui"]))
			{	
				if(edgeList[edge,"ui"] != "") # partial correlation only if ui exists
				{
					ui.name = unlist(strsplit(edgeList[edge,"ui"], ","))
					ui = myGv$data[ , ui.name ]
					
					# Remove NA samples on x,y and ui
					NaVal =  sort(unique(which(is.na(cbind(x,y,ui)), arr.ind = T)[,1]))
				
					if(length(NaVal) != 0)
					{
						if(length(ui.name) > 1) { ui = ui[-NaVal,] }
						else { ui = ui[-NaVal] }

						tmp.pcor = try(pcor.test(x[-NaVal],y[-NaVal], ui, method = "pearson")$"estimate", silent = T) 
						if (class(tmp.pcor) == "try-error") #--- Some of the ui do not have a consistent pcor
						{
							ui.tokeep = c()
							for( n in ui.name)
							{
								if(class(try(pcor.test(x[-NaVal],y[-NaVal], ui[,n], method = "pearson")$"estimate", silent = T)) != "try-error") 
								{
									ui.tokeep = c(ui.tokeep, n)
								}
							}
							ui = ui[,ui.tokeep]	
							tmp.pcor = pcor.test(x[-NaVal],y[-NaVal], ui, method = "pearson")$"estimate" 
						}
					}
					else 
					{
						tmp.pcor = pcor.test(x,y, ui, method = "pearson")$"estimate"
					}
				} 
				else # if no ui, we simply compute the spearman's correlation
				{
					NaVal = sort(unique(which(is.na(cbind(x,y)),arr.ind = T)[,1]))
					if(length(NaVal) != 0){ tmp.pcor = cor(x[-NaVal],y[-NaVal], method = "pearson") }
					else {tmp.pcor = cor(x,y, method = "pearson")}
				}
			}
			else # if no ui, we simply compute the spearman's correlation
			{
				NaVal = sort(unique(which(is.na(cbind(x,y)),arr.ind = T)[,1]))
				if(length(NaVal) != 0){ tmp.pcor = cor(x[-NaVal],y[-NaVal], method = "pearson") }
				else {tmp.pcor = cor(x,y, method = "pearson")}
			}

			# Test
			if(sign(tmp.pcor) == 1){ edgeList.signed[edge,"Sign"] <- "+" }
			else{ edgeList.signed[edge,"Sign"] <- "-" } 
			edgeList.signed[edge,"Diff"] <- (tmp.pcor)
		}
			
	}
	return(edgeList.signed[,c(ncol(edgeList.signed)-1, ncol(edgeList.signed))])
}

