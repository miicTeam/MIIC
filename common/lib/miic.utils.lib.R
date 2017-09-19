discoNet.parseCommandLine <- function(userDir, rootDir)
{
  setwd(userDir)
    #### Specify the input arguments
    mySpec = matrix( c(  'param', 'p', 1, "character"	    # String of arguments for a specific method [param1:val1,param2:val2]
						# For all
						, 'method', 'm', 1, "character"		# Discovery method: miic, aracne, pc, hc (default: miic)
						, 'stepsToDo', 'd', 1, "character"	# Steps that should be performed: 1.skel, 2.ort, 3.sumry., 4.plot [all by dflt, else -d 1,3 etc...]

						, 'inputData', 'i', 1, "character"	# File path to the discrete input data
						, 'outDirPath', 'o', 1, "character" # Output dir path
						, 'trueEdges', 't', 1, "character"	# File path to the true egdes (for summary)
						, 'weights', 'w', 1, "character"	# File path to the sample weights
						, 'blackBox', 'x', 1, "character"	# File path to the blackbox file
						, 'shfPvalue', 'c', 1, "character"	# shf:100,aph:0.01
						, 'stateOrder', 's', 1, "character" # File path to the each variable state definition (mainly used to sign the edges)
						, 'layout', 'l', 1, "character"     # File path to the vertices layout (for plot)
						, 'writeNetwork', 'n', 0, "character"     # Write XGMML?
						, 'isBent', 'b', 0, "logical" )       # Curved edges (for plot, default is FALSE)	
   					, byrow = TRUE, ncol = 4 );

    #### Read the input arguments
    opt = getopt(mySpec)
    #### Variable for the arguments
	# Specific miic
	myEffN = -1
	myCplx = "nml"	        # (1) NML by default
	myEta = 1	            # default at 1
	myShuffle = 0
	myUnorient = FALSE
	myLatent = FALSE
	myTplReuse = TRUE
	myK23 = TRUE
	myDegeneracy = FALSE
	myPropag = TRUE
	myNoInitEta = FALSE
	myHvs = 0
	nThreads = 1
	mySeed = 0
	
	# Specific aracne
	myEpsilon = 0
	myEstimator = "mi.empirical"
	
	# Specific bayes hc
	myScore = "bic"
	myRestart = 20
	
	# Specific cb pc
	myIndepTest = "disCItest" # possible: gaussCItest, dsepTest, disCItest, binCItest; default: gaussCItest
	myAlpha = 0.01
	myU2pd = "relaxed"          # possible: "relaxed", "rand", "retry"
	myAmbg = "majority"         # possible: "conservative", "majority"
	myConflict = "ignore"        # possible: "solve", "ignore"
	mySkelMeth = "stable"	    # possible: "stable", "original", "stable.fast"
	
	# Specific fci
	myFciType = "normal"
	myMajRule = FALSE
	myConserv = FALSE
	
	# Specific fci
	myIndepTest = "disCItest" # possible: gaussCItest, dsepTest, disCItest, binCItest; default: gaussCItest
	myAlpha = 0.01
	myU2pd = "relaxed"          # possible: "relaxed", "rand", "retry"
	myConflict = "ignore"        # possible: "solve", "ignore"
	mySkelMeth = "stable"	    # possible: "stable", "original", "stable.fast"
	
	myFciType = "normal"
	myMajRule = FALSE
	myConserv = FALSE
	


	# For all
	myMethod = "miic"
	mySteps.vect = c(1, 2, 3, 4)
	myInputData = character(0)
	myOutDirPath = character(0)
	

	#myTplDirCycle = "accept"

	myTrueEdgesFilePath = character(0)
	myBlackBox = character(0)
	myConfidenceRatio = character(0)
	myCutShf = character(0)
	myStateOrder = character(0)
	myLayout = character(0)
	myDoCurve = FALSE
	myEntropy= FALSE
	myIsVerbose = FALSE
	myWriteNetwork = FALSE

  	#### COMMON ARGUMENTS
  	# ----
	# The method used for the reconstruction
    if( !is.null( opt$method ) ){ myMethod <- opt$method } else { cat("\t# Reminder -> The default learning method is miic (-m)\n")}

	# The steps that should be performed
    if( !is.null( opt$stepsToDo ) ){ mySteps.vect <- as.numeric( unlist( strsplit(opt$stepsToDo, split=",") ) ) }

	# INPUTData and OUTPUT dir path are compulsory
    if( !is.null( opt$inputData ) ){ myInputData <-  normalizePath(opt$inputData) } else { stop("The input data file is required (-i)\n") }
	if( !is.null( opt$outDirPath ) ){ myOutDirPath <-  suppressWarnings(normalizePath(opt$outDirPath)) } else { stop("The output directory path is required (-o)\n") }
	if( !is.null( opt$trueEdges ) ){ myTrueEdgesFilePath <- normalizePath(opt$trueEdges) } else { cat("\t# Reminder -> No information about the true edges (-t)\n") }

	# Get the variables states orders if given
	if( !is.null( opt$stateOrder ) ){ myStateOrder <- normalizePath(opt$stateOrder) } else { cat("\t# Reminder -> No information about the variables states order (-s)\n") }

	# Get the Black box file if provided
    if( !is.null( opt$blackBox ) ){ myBlackBox <- normalizePath(opt$blackBox) } else { cat("\t# Reminder -> No restrictions about the edges inference (blackbox) (-bb)\n") }

    if( !is.null( opt$shfPvalue ) )
    { 
		tmp.param <- opt$shfPvalue

		# Remove extra white spaces
		tmp.param = gsub( "\\s+", "", tmp.param )

		# Parameters are given as: "param1:val1,param2:val2"
		tmp.param.split = sapply( unlist( strsplit(tmp.param, split=',') ), strsplit, split = ':' )

        # Check that all argument has its value
        tmp.checkPairParamValue = unlist( lapply(tmp.param.split, length) )
        names(tmp.checkPairParamValue) = c()
        if( all.equal( tmp.checkPairParamValue, rep( 2, length(tmp.param.split) ) ) != TRUE )
        { stop("Some values are missing for the complementary parameters (-c)\n") }
		names( tmp.param.split ) = sapply(tmp.param.split, "[[", 1)
		tmp.param.split = sapply(tmp.param.split, "[[", 2)

		if( "ccr" %in% names(tmp.param.split)){ myConfidenceRatio = as.numeric( tmp.param.split[["ccr"]] ) }
		if( "csh" %in% names(tmp.param.split)){ myCutShf = as.numeric( tmp.param.split[["csh"]] ) }

    } 
    else if( myMethod == "miic" )
    { 
    	cat("\t# Reminder -> No Confidence ratio filter (-c)\n") 
    }

	# Get the layout of the nodes if defined
	if( !is.null( opt$layout ) ){ myLayout <- normalizePath(opt$layout) } else { cat("\t# Reminder -> No information about the nodes layout (-l)\n") }
	
	# Set if yser wants the graphml or xgmml
	if( !is.null( opt$writeNetwork ) ){ 
	    myWriteNetwork <- TRUE
  } else { cat("\t# Reminder -> Network will not be written in graphml format (-n)\n") }
  
	# Should the plotted edges be curved?
	if( !is.null(opt$isCurved) ) {myDoCurve <- TRUE }

	# Should the confidence gradient be weighetd by entropies?
	if( !is.null( opt$isEntropy ) ) { myEntropy <- TRUE }

	# Give details?
    if( !is.null(opt$isVerbose) ) { myIsVerbose <- TRUE }
	# ----
	#### Check the input arguments specific to a method
    # ----

	if( !is.null( opt$param ) )
	{

		tmp.param <- opt$param

		# Remove extra white spaces
		tmp.param = gsub( "\\s+", "", tmp.param )

		# Parameters are given as: "param1:val1,param2:val2"
		tmp.param.split = sapply( unlist( strsplit(tmp.param, split=',') ), strsplit, split = ':' )

        # Check that all argument has its value
        tmp.checkPairParamValue = unlist( lapply(tmp.param.split, length) )
        names(tmp.checkPairParamValue) = c()
        if( all.equal( tmp.checkPairParamValue, rep( 2, length(tmp.param.split) ) ) != TRUE )
        { stop("Some values are missing for the complementary parameters (-p)\n") }
		names( tmp.param.split ) = sapply(tmp.param.split, "[[", 1)
		tmp.param.split = sapply(tmp.param.split, "[[", 2)

		# ---- miic ----
		if( myMethod == "miic" )
		{

			if( "efn" %in% names(tmp.param.split)){ myEffN = as.numeric( tmp.param.split[["efn"]] ) }
			if( "shf" %in% names(tmp.param.split)){ myShuffle = as.numeric( tmp.param.split[["shf"]] ) }
			if( "cpx" %in% names(tmp.param.split)){ myCplx = tmp.param.split[["cpx"]] }
			if( "eta" %in% names(tmp.param.split)){ myEta = as.numeric( tmp.param.split[["eta"]] ) }
			if( "urt" %in% names(tmp.param.split)){ if( tmp.param.split[["urt"]] == "yes" ){ myUnorient <- TRUE } }
			if( "lat" %in% names(tmp.param.split)){ if( tmp.param.split[["lat"]] == "yes" ){ myLatent <- TRUE } }
			if( "reu" %in% names(tmp.param.split)){ if( tmp.param.split[["reu"]] == "no" ){ myTplReuse <- FALSE } }
			if( "k23" %in% names(tmp.param.split)){ if( tmp.param.split[["k23"]] == "no" ){ myK23 <- FALSE } }
			if( "deg" %in% names(tmp.param.split)){ if( tmp.param.split[["deg"]] == "yes" ){ myDegeneracy <- TRUE } }
			if( "prg" %in% names(tmp.param.split)){ if( tmp.param.split[["prg"]] == "no" ){ myPropag <- FALSE } }
			if( "nie" %in% names(tmp.param.split)){ if( tmp.param.split[["nie"]] == "yes" ){ myNoInitEta <- TRUE } }
			if( "hvs" %in% names(tmp.param.split)){ if( tmp.param.split[["hvs"]] == "yes" ){ myHvs <- 1 } }
			if( "hvs" %in% names(tmp.param.split)){ if( tmp.param.split[["hvs"]] == "no" ){ myHvs <- 0 } }
			if( "thr" %in% names(tmp.param.split)){ nThreads = tmp.param.split[["thr"]] }
			if( "sed" %in% names(tmp.param.split)){ mySeed = as.numeric( tmp.param.split[["sed"]]) }

		} else if( myMethod == "fci" ) {
		  
		  if( "aph" %in% names(tmp.param.split)){ myAlpha = as.numeric( tmp.param.split[["aph"]] ) }
		  if( "cit" %in% names(tmp.param.split)){ myIndepTest = tmp.param.split[["cit"]]}
		  if( "typ" %in% names(tmp.param.split)){ myFciType = tmp.param.split[["typ"]] }
		  if( "skm" %in% names(tmp.param.split)){ mySkelMeth = tmp.param.split[["skm"]] }
		  if( "maj" %in% names(tmp.param.split)){ myMajRule = tmp.param.split[["maj"]] }
		  if( "con" %in% names(tmp.param.split)){ myConserv = tmp.param.split[["con"]] }
		} 
	}
	setwd(rootDir)
	return( list( argMethod=myMethod, argSteps=mySteps.vect, argUnorient=myUnorient
	              , argInData=myInputData, argOutDir=myOutDirPath, argThreads=nThreads, argSeed = mySeed
	              , argEffN=myEffN, argCplx=myCplx, argEta=myEta, argShuffle = myShuffle, argHvs = myHvs 
	              , argWriteNetwork = myWriteNetwork, argBlackBox=myBlackBox, argStateOrder = myStateOrder
	              , argConfidenceRatio = myConfidenceRatio, argCutShf = myCutShf
	              , argLatent = myLatent, argTplReuse = myTplReuse, argK23 = myK23, argDegeneracy = myDegeneracy, argPropag = myPropag, argNoInitEta = myNoInitEta
	              , argTrueEdgesFile = myTrueEdgesFilePath, argLayout = myLayout, argDoCurve = myDoCurve, argEntropy=myEntropy
	              , argEstimator = myEstimator, argEpsilon = myEpsilon
	              , argScore = myScore, argRestart = myRestart
	              , argAlpha = myAlpha, argCItest = myIndepTest, argU2pd = myU2pd, argSkelMeth = mySkelMeth, argAmbg = myAmbg, argConflict = myConflict
	              , argFciType = myFciType, argMajRule = myMajRule, argConserv = myConserv
	              , argVerbose=myIsVerbose ) )
}
 
cleanOutput <- function(outDirPath){
	file_list_logs = list.files(outDirPath, pattern = "log.", full.names=TRUE)
	file_list_edges = list.files(outDirPath, pattern = "edges.txt", full.names=TRUE)
	file_list_arrow = list.files(outDirPath, pattern = "Arrowhead.txt", full.names=TRUE)
	file_list = c(file_list_logs, file_list_edges, file_list_arrow)
	invisible(lapply(file_list, file.remove))
}



isContinuous <- function(tmpArray){
	nbLevels = unique(tmpArray)
	nbSamples = length(tmpArray)
	if(length(nbLevels) >= 0.8 * nbSamples | length(nbLevels) >= 100 ){return(1) }
	else{ return(0) }
}

checkInput <- function(dataFile, method, data){
	errCode = "0"
	isCnt_test = 0
	

	str = readChar(dataFile, file.info(dataFile)$size)
	str_names = unlist(strsplit(str,"\n"))[1]
	if(grepl("#", str) | grepl("&", str) | grepl("<", str_names) | grepl(">", str_names) | grepl("\"", str) | grepl("'", str)){
		errCode = "111"
	}
	else{

		isCnt_test = sum(apply( data[-1,],2,isContinuous ))
  		   
		if(isCnt_test > 0){ errCode = "118"}
		if (method == "miic" & length(unique(data[1,])) != ncol(data)) { errCode = "117"}
	}
	
	return(errCode)
}


checkTrueEdges <- function(edgesFile, inputData.df){
	errCode = "0"
	if( !file.exists(edgesFile) ){errCode="020"}
	else{
		str = readChar(edgesFile, file.info(edgesFile)$size)
		if(grepl("#", str) | grepl("&", str) | grepl("<", str) | grepl(">", str) | grepl("\"", str) | grepl("'", str)){
			errCode = "121"
		}
		else {
			data = try( read.table(edgesFile, header=F, as.is=F, sep='\t', check.names = F) )
			if(class(data) == "try-error"){errCode = "021"}
			else if(ncol(data) != 2 & ncol(data) != 3){errCode = "023"}	
			else if(length(which(!c(as.vector(data[,1]), as.vector(data[,2])) %in% colnames(inputData.df))) > 0){

				print(colnames(inputData.df))

				errCode = "025"}
		}
	}
		
	return(errCode)
}

checkLayout <- function(layoutFile){
	errCode = "0"
	if( !file.exists(layoutFile) ){errCode="030"}
	else{
		str = readChar(layoutFile, file.info(layoutFile)$size)
		if(grepl("#", str) | grepl("&", str) | grepl("<", str) | grepl(">", str) | grepl("\"", str) | grepl("'", str)){
			errCode = "131"
		}
		else{
			data = try(read.table(layoutFile, header=F, as.is=T, sep='\t', check.names = F) )
			if( class(data) == "try-error" ){errCode = "031"}
			else{
				if( ncol(data) != 2 & ncol(data) != 3  ){errCode = "033"}
	      else {
	  			if(ncol(data) == 2){
	  			  if(!is.numeric(data[,1]) | !is.numeric(data[,2]) ){errCode = "034"}
	  			} else if(ncol(data) == 3 & (!is.numeric(data[,2]) | !is.numeric(data[,3]) ) ) {errCode = "038"}
	      }
			}
		}
	}
	return(errCode)
}

checkStateOrder <- function(stateOrderFile,inData){
	errCode = "0"
	if( !file.exists(stateOrderFile) ){errCode="040"}
	else{
		str = readChar(stateOrderFile, file.info(stateOrderFile)$size)
		str_names = unlist(strsplit(str,"\n"))[1]
		if(grepl("#", str) | grepl("&", str) | grepl("<", str_names) | grepl(">", str_names) | grepl("\"", str) | grepl("'", str)){
			errCode = "141"
		} else{
			data = try(read.table(stateOrderFile, header=T, as.is=T, sep='\t', check.names = F))
			if( class(data) == "try-error" ){errCode = "041"}
			else {
				rownames(data) = data[,"var_names"]
				myVariables = rownames(data)

				if( ncol(data) != 2  & ncol(data) != 3){errCode = "043"}
			}
		}
	}
	return(errCode)
}

errorCodeToString <- function(error_code){
	errorList1 = list(
				'0' = "Warning:", 
				'1' = "Fatal error:"
			);

	errorList2 = list(
				'0' = 'Unknown Error',
				'1' = "inputFile",
				'2' = "trueEdgeFile", 
				'3' = "layoutFile",
				'4' = "stateOrderFile"
			);

	errorList3 = list(
		'0' = "does not exist", 
		'1' = "is not readable, check the file format. Special characters like #,&amp;,<,> are not allowed.",
		'2' = "has rownames, please remove them",
		'3' = "should have exactly two columns",
		'4' = "should be numerical",
		'5' = "does not correspond to the input data matrix. Please check node names.",
		'6' = "problem as states in stateOrderFile are not consistent with dataset",
		'7' = "has duplicated column names",
		'8' = "contains one or several variables with too many different levels (might be continuous data)",
		'9' = 'occured.'
	);
	error_string = unlist(strsplit(error_code, ""))
	return(paste(errorList1[[error_string[1]]], errorList2[[error_string[2]]], errorList3[[error_string[3]]], sep = " "))
}