#!/usr/bin/Rscript

#### Delete all variables
rm( list = ls ( all = TRUE ) )

startTime.whole <- proc.time()
cat( "\t# --------\n# -> START Summary...\n" )
options(warn=-1)
# If testing
myTest = FALSE

##### Load libraries
require(getopt, quietly = TRUE, warn.conflicts = FALSE)            # Parse command line options
require(bnlearn, quietly = TRUE, warn.conflicts = FALSE)
require(tools, quietly = TRUE, warn.conflicts = FALSE)

library(ppcor)


### FOR TEST ####
if( myTest == TRUE )
{ setwd( "~/DiscoNet/common" ) }
### FOR TEST ####

#### Load supplementary functions
source( file.path( "..", "sharedLib", "shared.utils.lib.R" ) )
source( file.path( "lib", "gmSummary.lib.R" ) )
source( file.path( "lib", "gmStatistics.skeleton.lib.R" ) )
source( file.path( "lib", "gmStatistics.orient.lib.R" ) )

### FOR TEST ####
if( myTest == TRUE )
{
    summaryArg.list.external = list()
    summaryArg.list.external[["argInData"]] = "~/Project_3off2/pipeline_test_lat/data_benchmark/child/input/rawData/child_1000/child_1000_0001.txt"
    summaryArg.list.external[["argOutDir"]] = "~/Project_3off2/pipeline_test_lat/data_benchmark/child/output/test1"
    summaryArg.list.external[["stateOrder"]] = NULL
    summaryArg.list.external[["argTrueEdgesFile"]] = "~/Project_3off2/pipeline_test_lat/data_benchmark/child/input/graphData/edges_name.txt"
    summaryArg.list.external[["argInfEdgesFile"]] = "edgesList.miic.txt"
    summaryArg.list.external[["argAdjMat"]] = "adjacencyMatrix.miic.orientProba.txt"
    summaryArg.list.external[["argCnt"]] = 0

    summaryArg.list.external[["argVerbose"]] = TRUE
}
### FOR TEST ####

#### Parse command line arguments
#### ----
summaryArg.list = list();
if( !exists( "summaryArg.list.external" ) )
{
	summaryArg.list = summary.parseCommandLine()
	if( summaryArg.list[["argVerbose"]] == TRUE )
	{
	    cat( "\t# --------> No Arguments from any variable inside the environment\n" )
	}
	
} else {
	
	summaryArg.list = summaryArg.list.external
    if( summaryArg.list[["argVerbose"]] == TRUE )
	{
	    cat( "\t# --------> Arguments from a variable inside the environment\n" )
	}

}

#### Set Global variables
inputDataFile = summaryArg.list[["argInData"]]
outDirPath = summaryArg.list[["argOutDir"]]
stateOrderFile = summaryArg.list[["stateOrder"]]
trueEdgesFile = summaryArg.list[["argTrueEdgesFile"]]
layoutFile = summaryArg.list[["argLayoutFile"]]
infEdgesFile = summaryArg.list[["argInfEdgesFile"]]
adjMatFile = summaryArg.list[["argAdjMat"]]
isCnt = summaryArg.list[["argCnt"]]
isVerbose = summaryArg.list[["argVerbose"]]
confRatio = as.numeric(summaryArg.list[["argConfRatio"]])

#read dataset
inputData.df <- read.table( file = inputDataFile, header = TRUE, stringsAsFactors = FALSE, sep = "\t", 
                            na.strings = c('', NA), check.names=FALSE )

#### Load the list of the true edges
true.edgesList = NULL
if( length(trueEdgesFile) > 0 )
{
    if( nchar(trueEdgesFile) > 0 )
    { true.edgesList <- read.table( file = trueEdgesFile, header = FALSE, stringsAsFactor = FALSE, sep='\t', check.names=FALSE) 
    }
}

originalPath = getwd()
setwd( outDirPath )

#### Define a log file
logFilePath = "log.summary.txt"
if( file.exists( logFilePath ) ){ file.remove( logFilePath ) }
sink( file=logFilePath, append = TRUE, split = TRUE )

#### Set an env variable
gV <- new.env()

#### Load fist the adjacency matrix to get all the variables name
adjMat <- read.table( file = adjMatFile, header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = "\t", check.names=FALSE )
gV$allProperties = colnames(adjMat)


if(!is.null(true.edgesList)){
  #### Set as rowname the key of each edge
  rownames( true.edgesList ) <- apply( true.edgesList, MARGIN = 1 , function(myRow){ binVectToStr( myRow[1:2], gV ) } ) 
}

#### Load also the  inferred edges with their mutual information
allMutInfo.df <- read.table(file= infEdgesFile,header= TRUE,row.names= 1,stringsAsFactors= FALSE,sep = "\t", check.names=FALSE)
  

#### Check if the orientations have been learned with assuming latent variables
gV$isLatent = FALSE
if( max(adjMat) > 2 | min(adjMat) < -2 ){ gV$isLatent = TRUE }



#### Convert the adj mat it to upper triangular matrix
adjMat[ lower.tri( adjMat, diag = TRUE ) ] = 0

#### Get all inferred edges and gather their idx and val into one list
edgesNbr = 0
edges.list = list( nOrt = list( idx = data.frame( nrow = numeric(0), ncol = numeric(0) ), val = numeric(0), key = character(0) )
                    , ort = list( idx = data.frame( nrow = numeric(0), ncol = numeric(0) ), val = numeric(0), key = character(0) )
                    , ph = list( idx = data.frame( nrow = numeric(0), ncol = numeric(0) ), val = numeric(0), key = character(0) ) )

#### ----
# -- not oriented
if( isVerbose == TRUE ){ cat("\t# -- Get not oriented edges idx & val\n") }
edges.list[["nOrt"]][["idx"]] <- which( adjMat == 1, arr.ind = TRUE )
edges.list[["nOrt"]][["val"]] <- rep( 1, nrow(edges.list[["nOrt"]][["idx"]]) )
edgesNbr = edgesNbr + length( edges.list[["nOrt"]][["val"]] )

# -- oriented
if( isVerbose == TRUE ){ cat("\t# -- Get oriented edges idx & val\n") }
#### get the number of oriented edges
ortN = ( nrow( which( abs( adjMat ) >= 2, arr.ind = TRUE ) ) )

if( ortN > 0 )
{
    #### Initialize the dimension array of index with the nbr of oriented edges
    edges.list[["ort"]][["idx"]] <- data.frame( nrow = numeric(ortN), ncol = numeric(ortN) )
    edges.list[["ort"]][["val"]] <- rep(NA, ortN)

    countRow = 0
    for( iOrt in c(2,-2,4,-4,6) )
    {
        #### Get the array of index
        tmp.ortArr.idx <- which( adjMat == iOrt, arr.ind = TRUE )
        
        #### Get the number of edges
        tmp.ortN = nrow( tmp.ortArr.idx )
        tmp.ortArr.val <- rep( iOrt, (tmp.ortN) )
        
        #### Insert in the initalised matrix
        if( tmp.ortN > 0 )
        {
            edges.list[["ort"]][["idx"]][(countRow+1):(countRow+tmp.ortN),] <- tmp.ortArr.idx
            edges.list[["ort"]][["val"]][(countRow+1):(countRow+tmp.ortN)] <- tmp.ortArr.val
            
            countRow = (countRow + tmp.ortN)
        }
    }
}
edgesNbr = edgesNbr + length( edges.list[["ort"]][["val"]] )

# -- phantom
if( isVerbose == TRUE ){ cat("\t# -- Get phantom edges idx & val\n") }
edges.list[["ph"]][["idx"]] <- which( adjMat == 0, arr.ind = TRUE )

if( nrow( edges.list[["ph"]][["idx"]] ) > 0 )
{
	edges.list[["ph"]][["idx"]] <- edges.list[["ph"]][["idx"]][which( edges.list[["ph"]][["idx"]][,2] > edges.list[["ph"]][["idx"]][,1] ),, drop=FALSE]
	edges.list[["ph"]][["val"]] <- rep( 0, nrow( edges.list[["ph"]][["idx"]] ) )
}
edgesNbr = edgesNbr + length( edges.list[["ph"]][["val"]] )

#### Initialize the output data frame 
### LV // adding 2 columns: signs & partial_correlation
outputSummary.df <- data.frame( x = character(edgesNbr), y = character(edgesNbr)
                                , type = character(edgesNbr), ui = character(edgesNbr), info = numeric(edgesNbr)
								, cplx = numeric(edgesNbr), Nxy_ui = numeric(edgesNbr)
								, confidence = numeric(edgesNbr)
                                , infOrt = numeric(edgesNbr), trueOrt = numeric(edgesNbr)
								, isOrt = character(edgesNbr), isOrtOk = character(edgesNbr)
                                , essential = character(edgesNbr)
						  , sign = character(edgesNbr)
						  , partial_correlation = character(edgesNbr)
                                , stringsAsFactors = FALSE )

#### ----
countRow <- 0

for( iEdgeCategory in names( edges.list ) )
{
    if( isVerbose == TRUE ){ cat("\t# ---->", iEdgeCategory, "\n" ) }
    if( nrow( edges.list[[iEdgeCategory]][["idx"]] ) > 0 )
    {        
        #### Replace numbers by names
        edges.list[[iEdgeCategory]][["idx"]] = t(apply( edges.list[[iEdgeCategory]][["idx"]], MARGIN = c(1), FUN = function(x) { x = gV$allProperties[x] } ))
        
        #### Compute the key of each edge
        edges.list[[iEdgeCategory]][["key"]] = t(apply( edges.list[[iEdgeCategory]][["idx"]], MARGIN = c(1), FUN = function(x) { x = binVectToStr( x, gV ) } ))
        
        #### Set the bounds for insertion
        currentN = nrow(edges.list[[iEdgeCategory]][["idx"]])
        lastRow = ( countRow + currentN )                

        #### Insert the new data
        outputSummary.df[c((countRow+1):lastRow), ] <- data.frame( 
                                x = as.character( as.vector( edges.list[[iEdgeCategory]][["idx"]][,1] ) )
                                , y = as.character( as.vector( edges.list[[iEdgeCategory]][["idx"]][,2] ) )
                                , type = rep( NA, currentN )
								, ui = rep( NA, currentN )
                                , info = rep( NA, currentN )
								, cplx = rep( NA, currentN )
								, Nxy_ui = rep( NA, currentN )
								, confidence = rep( NA, currentN )
                                , infOrt = edges.list[[iEdgeCategory]][["val"]]
                                , trueOrt = rep( NA, currentN )
                                , isOrt = rep( NA, currentN ), isOrtOk = rep( NA, currentN )
                                , essential = rep( NA, currentN )
						  , sign = rep( NA, currentN )
						  , partial_correlation = rep( NA, currentN )
                                , stringsAsFactors = FALSE )
                                
        #### Update the total row count
        countRow = lastRow
    }
}

#### Set the key of each edge as rowname
rownames( outputSummary.df ) <- do.call(c, list(edges.list$nOrt$key, edges.list$ort$key, edges.list$ph$key))

  
#### Get the true links and check the TP, FP, FN, TN
if( !is.null(true.edgesList) )
{
    #### Set the name of the columns for the true edges list
    colnames( true.edgesList )[1:2] = c( "x", "y" )
    if( ncol(true.edgesList) == 3 ){ colnames( true.edgesList )[3] = "essential" }
       
    #### Distinguish the TP from the FP among the non phantom links
    allEdge.idx <- which( abs( outputSummary.df[,"infOrt"] ) > 0 )

    # [TP]    
    TP.type.idx = c()
    if( length(allEdge.idx) > 0 )
    {
        #### First, initialise all inferred edge as FP
        outputSummary.df[allEdge.idx, "type"] <- 'FP' 
    
        #### Then, find the TP
        #### ie., among the inferred edges, which are in the true edges list
        TP.type.idx <- which( abs( outputSummary.df[,"infOrt"] ) > 0 & ( rownames( outputSummary.df ) %in% rownames( true.edgesList ) ) )

        if( length(TP.type.idx) > 0 ) { outputSummary.df[TP.type.idx, "type"] <- 'TP' }
    }
    #### ----
    
    #### Distinguish the TN from the FN among the phantom links
    allPh.idx <- which( outputSummary.df[,"infOrt"] == 0 )
    
    # [FN]
    FN.type.idx = c()
    if( length(allPh.idx) > 0 ) 
    {
        #### First, initialise all inferred phantom as TN
        outputSummary.df[allPh.idx, "type"] <- 'TN'
    
        #### Then, find the FN
        #### ie., among the phantom, which have a key found in the true edges list
        FN.type.idx <- which( rownames( true.edgesList ) %in% rownames(outputSummary.df)[allPh.idx] )
    
        if( length(FN.type.idx) > 0 ) { outputSummary.df[rownames( true.edgesList )[FN.type.idx], "type"] <- 'FN' }
    }
    #### ----

    #### Add the true orientations into the summary
    #### ----
    if( length(FN.type.idx) > 0 )
    {
        #### For the FN (forward = 2, backward = -2)
        FN.ort <- sapply( rownames( true.edgesList )[FN.type.idx]
            , FUN=function(myKey)
            {
                ifelse( isTRUE(all.equal(outputSummary.df[myKey, c("x", "y")], true.edgesList[myKey, c("x", "y")] )), 2, -2 )
            } )
        if( length( names( FN.ort ) ) > 0)
        {
            outputSummary.df[names(FN.ort), "trueOrt"] = FN.ort
        }
    }
    
    if( length(TP.type.idx) > 0 )
    {
        #### For the TP (forward = 2, backward = -2)
        TP.ort <- sapply( rownames( outputSummary.df )[TP.type.idx]
            , FUN=function(myKey)
            {
                ifelse( isTRUE(all.equal(outputSummary.df[myKey, c("x", "y")], true.edgesList[myKey, c("x", "y")] )), 2, -2 )
            } )
        if( length( names( TP.ort ) ) > 0)
        {
            outputSummary.df[names(TP.ort), "trueOrt"] = TP.ort
        }
    }
    #### ----

    #### Get the information about the correct orientation
    #### ----
	# -- oriented TP edge
	ortEdge.idx = which( ( abs( outputSummary.df[, "infOrt"] ) > 1 ) & ( outputSummary.df[, "type"] == 'TP' ) )
	if( length( ortEdge.idx ) > 0 ) { outputSummary.df[ortEdge.idx, "isOrt"] = 'Y' }

    # -- compatible orientation
    compatibleOrt.idx = which( ( outputSummary.df[, "type"] == 'TP') &
                               ( ( sign( outputSummary.df[, "infOrt"] ) == sign( outputSummary.df[, "trueOrt"] ) ) |
                                 ( outputSummary.df[, "infOrt"] %in% c(1,6) ) ) )
    if( length( compatibleOrt.idx ) > 0 ) { outputSummary.df[compatibleOrt.idx, "isOrtOk"] = 'Y' }

    # -- not compatible orientation
    notCompatibleOrt.idx = which( ( outputSummary.df[, "type"] == 'TP' ) & ( is.na( outputSummary.df[, "isOrtOk"] ) ) )
    if( length( notCompatibleOrt.idx ) > 0 ) { outputSummary.df[notCompatibleOrt.idx, "isOrtOk"] = 'N' }
    #### ----

    #### Add the essentiality info if any
    if( "essential" %in% colnames(true.edgesList) )
    { outputSummary.df[rownames(true.edgesList), "essential"] <- true.edgesList[,"essential"] }

 
} else {
    
    #### First, initialise all inferred edge as N
    outputSummary.df[, "type"] <- 'N' 
    
    #### Then, get the index of oriented edges and set them to P
    remainingEdges.idx <- which( abs(outputSummary.df[, "infOrt"]) > 0 )
    if( length(remainingEdges.idx) > 0 )
    { outputSummary.df[rownames(outputSummary.df[remainingEdges.idx,]), "type"] <- 'P' }
}
#### Use the mutual information values
#### ----
#### Add the mutual information values, the cplx values and the partial_correlationerence (ie, confidence) if the cplx exisits...
if( "cplx" %in% colnames( allMutInfo.df ) )
{
	outputSummary.df[rownames(outputSummary.df), c( "info", "cplx")] <- allMutInfo.df[rownames(outputSummary.df), c("Ixy_ui", "cplx")]
	outputSummary.df[, "confidence"] = as.numeric(outputSummary.df[, "info"] - outputSummary.df[, "cplx"] )

} else {
	outputSummary.df[rownames(outputSummary.df), c( "info")] <- allMutInfo.df[rownames(outputSummary.df), c("Ixy_ui")]
}

#### Add the Nxy_ui if it exists
if( "Nxy_ui" %in% colnames( allMutInfo.df ) )
{ outputSummary.df[rownames(outputSummary.df), "Nxy_ui"] <- allMutInfo.df[rownames(outputSummary.df), "Nxy_ui"] }


#### Add the {ui}
if( "ui.vect" %in% colnames( allMutInfo.df ) )
{ 

    outputSummary.df[rownames(outputSummary.df), "ui"] <- allMutInfo.df[rownames(outputSummary.df), "ui.vect"] 
}

#### Order by decreasing value of confidence
outputSummary.df <- outputSummary.df[order(as.numeric(outputSummary.df[,"confidence"]), decreasing = TRUE),]
### Compute the sign of each edge and fill the two last columns (sign and partial_correlation)
cat("\t## ---- Computing the sign of the edges ---- \n")
#### Remove the lines that are all 'NA'
allNAs.idx = which( rowSums( is.na( inputData.df ) ) == ncol( inputData.df ) )
if( length( allNAs.idx ) > 0 ) { inputData.df = inputData.df[-allNAs.idx,] }

gV$data = inputData.df
cat("\t# Signs are calculated using partial correlation coefficient.\n")
# Order data according to stateOrderfile
# dataToStateOrder( gV, stateOrderFile, outDirPath )
# Calculate sign using partial correlation
outputSummaryFile = sub( ".txt$", ".summary.txt", infEdgesFile )



# if(correlationToEvaluate)
outputSummary.df[,c("sign","partial_correlation")] <- computeSign.cont.pcor(outputSummary.df, gV, outDirPath, originalPath)

# outputSummary.df[,c("sign","partial_correlation")] <- computeSign.pcor(outputSummary.df, gV)


#### Save the summary


colnames(outputSummary.df) = gsub("ui", "ai", colnames(outputSummary.df))
colnames(outputSummary.df) = gsub("confidence", "log_confidence", colnames(outputSummary.df))

write.table( outputSummary.df, file = outputSummaryFile, col.names = TRUE, row.names = TRUE, sep = '\t', quote = FALSE )

#### Compute complementary stats
#### ----
#### Open the final summary and count the number of TP, FP, TN, FN, TPnort, precision, recall, fscore
retListSkeleton = list()
retListSkeleton = computeSingleSampleStat( myFinalSummaryFilePath = outputSummaryFile )

#### SHD2 for ORIENTATIONS
if( length(summaryArg.list[["argTrueEdgesFile"]]) > 0 ){
  infMethod = "miic"
  retListOrient = list()
  outputSummary.df = computeSHD2SampleStat( myTrueGraphFilePath = summaryArg.list[["argTrueEdgesFile"]], myFinalSummaryFilePath = outputSummaryFile, skeletonStats = retListSkeleton, myInfMethod = infMethod, myAllProperties = gV$allProperties )
  cat("\t# Summary File will contain information about the comparison with the true graph provided (-t enabled)\n")
}
getridOff = c("isOrt","essential")
write.table( outputSummary.df[,!colnames(outputSummary.df)%in%getridOff], file = outputSummaryFile, col.names = TRUE, row.names = TRUE, sep = '\t', quote = FALSE )

sink(file=NULL)

spentTime.whole <- (proc.time() - startTime.whole)
if( isVerbose == TRUE ){ cat( "\t# -> STOP miic summary elapsed time:", (spentTime.whole[["elapsed"]]/60), "min\n\n" ) }



## NS
tmp_argOutDir = getwd()
tmp_pvalFile = file.path(tmp_argOutDir, "confRatios.txt")
if( file.exists(tmp_pvalFile) ){
    tmp_sum = outputSummary.df[,!colnames(outputSummary.df)%in%getridOff]
    conf_col = rep(1, nrow( tmp_sum ) )
    isCut = rep(NA, nrow(tmp_sum))
    tmp_sum = cbind(tmp_sum,conf_col,isCut)

    tmp_sum = cbind(tmp_sum,conf_col)
    colnames(tmp_sum) = c('x','y','type','ai','info','cplx','Nxy_ai','log_confidence','infOrt','trueOrt', 'isOrtOk', 'sign','partial_correlation','confidence_ratio')
    tmp_pval = read.table(tmp_pvalFile, header=T, as.is=T, sep="\t", check.names = F)
    tmp_pval[,"confidence_ratio"] = as.numeric(tmp_pval[,"confidence_ratio"])

    for(r in 1:nrow(tmp_pval))
    {
      tmp_sum[which(tmp_sum[,"x"] == tmp_pval[r,"x"] & tmp_sum[,"y"] == tmp_pval[r,"y"]),'confidence_ratio'] =tmp_pval[r,"confidence_ratio"]
      if(tmp_pval[r,"confidence_ratio"] < confRatio){
        tmp_sum[which(tmp_sum[,"x"]==tmp_pval[r,"x"] & tmp_sum[,"y"]==tmp_pval[r,"y"]),'isCut']='N'
      }
      else{
        tmp_sum[which(tmp_sum[,"x"]==tmp_pval[r,"x"] & tmp_sum[,"y"]==tmp_pval[r,"y"]),'isCut']='Y'
      }
    }
    tmp_sum = tmp_sum[,c('x','y','type','ai','info','cplx','Nxy_ai','log_confidence','confidence_ratio','infOrt','trueOrt', 'isOrtOk', 'sign','partial_correlation','isCut')]
    write.table(tmp_sum,outputSummaryFile, col.names=T, row.names=T, sep='\t',quote=F)
    outputSummary.df=tmp_sum
    rm(tmp_sum); rm(tmp_pval)
  }

##################################### NETWORK IN GRAPHML


line = "<graphml>\n"

#attributes part nodes

line = paste(line,"\t<key id=\"weight\" for=\"node\" attr.name=\"weight\" attr.type=\"double\">\n",sep="")
line = paste(line,"\t\t<default>0.2</default>\n",sep="")
line = paste(line,"\t</key>\n",sep="")
line = paste(line,"\t<key id=\"label\" for=\"node\" attr.name=\"label\" attr.type=\"string\"/>\n",sep="")

#attributes part edges

line = paste(line,"\t<key id=\"weight\" for=\"edge\" attr.name=\"weight\" attr.type=\"double\"/>\n",sep="")

line = paste(line,"\t<key id=\"label\" for=\"edge\" attr.name=\"label\" attr.type=\"string\"/>\n",sep="")

line = paste(line,"\t<key id=\"sourceArrowShape\" for=\"edge\" attr.name=\"sourceArrowShape\" attr.type=\"string\"/>\n",sep="")
line = paste(line,"\t<key id=\"targetArrowShape\" for=\"edge\" attr.name=\"targetArrowShape\" attr.type=\"string\"/>\n",sep="")
line = paste(line,"\t<key id=\"upstream\" for=\"edge\" attr.name=\"upstream\" attr.type=\"string\"/>\n",sep="")
line = paste(line,"\t<key id=\"info\" for=\"edge\" attr.name=\"info\" attr.type=\"double\"/>\n",sep="")
line = paste(line,"\t<key id=\"complexity\" for=\"edge\" attr.name=\"complexity\" attr.type=\"double\"/>\n",sep="")
line = paste(line,"\t<key id=\"nSamples\" for=\"edge\" attr.name=\"nSamples\" attr.type=\"int\"/>\n",sep="")
line = paste(line,"\t<key id=\"log_confidence\" for=\"edge\" attr.name=\"log_confidence\" attr.type=\"double\"/>\n",sep="")
if(confRatio != 1){
    line = paste(line,"\t<key id=\"confidenceRatio\" for=\"edge\" attr.name=\"confidenceRatio\" attr.type=\"double\"/>\n",sep="")
}
line = paste(line,"\t<key id=\"sign\" for=\"edge\" attr.name=\"sign\" attr.type=\"string\"/>\n",sep="")
line = paste(line,"\t<key id=\"partialCorrelation\" for=\"edge\" attr.name=\"partialCorrelation\" attr.type=\"double\"/>\n",sep="")
line = paste(line,"\t<key id=\"edgeType\" for=\"edge\" attr.name=\"edgeType\" attr.type=\"int\"/>\n",sep="")

line = paste(line,"\n",sep="")


line = paste(line,"\t<graph edgedefault=\"directed\">\n",sep="")

#cicle on nodes
for(node in colnames(inputData.df)){  
  line = paste(line,"\t\t<node id=\"", node, "\">\n",sep="")
  line = paste(line,"\t\t\t<data key=\"label\">",node, "</data>\n",sep="")
  line = paste(line,"\t\t\t<data key=\"weight\">0.5</data>\n",sep="")

  line = paste(line,"\t\t</node>\n",sep="")
}


line = paste(line,"\n",sep="")

indexes = which(outputSummary.df["type"] == "P" | outputSummary.df["type"] == "TP" | outputSummary.df["type"] == "FP")

#cicle on edges
for(index in indexes){
  if(!is.na(outputSummary.df[index, "log_confidence"]))
    weigth =  outputSummary.df[index, "log_confidence"] 
  else
    weigth = (outputSummary.df[index, "partial_correlation"]) 
  
  if(outputSummary.df[index, "infOrt"] == 1){
    
    
    line = paste(line,"\t\t<edge target=\"", outputSummary.df[index,2], "\" source=\"", outputSummary.df[index,1], "\" directed=\"false\">\n",sep="")

    line = paste(line,"\t\t\t<data key=\"label\">",  outputSummary.df[index,2], "---", outputSummary.df[index,1], "</data>\n",sep="")
    line = paste(line,"\t\t\t<data key=\"edgeType\">1</data>\n",sep="")

  } else   if(outputSummary.df[index, "infOrt"] == 2){
        
        line = paste(line,"\t\t<edge target=\"", outputSummary.df[index,2], "\" source=\"", outputSummary.df[index,1], "\" directed=\"true\">\n",sep="")
        line = paste(line,"\t\t\t<data key=\"sourceArrowShape\">none</data>\n",sep="")

        if(is.na(outputSummary.df[index, "partial_correlation"])){
            line = paste(line,"\t\t\t<data key=\"targetArrowShape\">arrow</data>\n",sep="")
            line = paste(line,"\t\t\t<data key=\"label\">",  outputSummary.df[index,2], "--&gt;", outputSummary.df[index,1], "</data>\n",sep="")
        }
        else {
            if(outputSummary.df[index, "partial_correlation"] > 0){
              line = paste(line,"\t\t\t<data key=\"targetArrowShape\">arrow</data>\n",sep="")
              line = paste(line,"\t\t\t<data key=\"label\">",  outputSummary.df[index,2], "--&gt;", outputSummary.df[index,1], "</data>\n",sep="")
            } else {
              line = paste(line,"\t\t\t<data key=\"targetArrowShape\">T</data>\n",sep="")
              line = paste(line,"\t\t\t<data key=\"label\">",  outputSummary.df[index,2], "--|", outputSummary.df[index,1], "</data>\n",sep="")
            }
        }
        line = paste(line,"\t\t\t<data key=\"edgeType\">2</data>\n",sep="")
         
   }  else   if(outputSummary.df[index, "infOrt"] == -2){
        
        line = paste(line,"\t\t<edge target=\"", outputSummary.df[index,1], "\" source=\"", outputSummary.df[index,2], "\" directed=\"true\">\n",sep="")
        line = paste(line,"\t\t\t<data key=\"sourceArrowShape\">none</data>\n",sep="")

        if(is.na(outputSummary.df[index, "partial_correlation"])){
            line = paste(line,"\t\t\t<data key=\"targetArrowShape\">arrow</data>\n",sep="")
            line = paste(line,"\t\t\t<data key=\"label\">",  outputSummary.df[index,1], "--&gt;", outputSummary.df[index,2], "</data>\n",sep="")
        }
        else {
            if(outputSummary.df[index, "partial_correlation"] > 0){
              line = paste(line,"\t\t\t<data key=\"targetArrowShape\">arrow</data>\n",sep="")
              line = paste(line,"\t\t\t<data key=\"label\">",  outputSummary.df[index,1], "--&gt;", outputSummary.df[index,2], "</data>\n",sep="")
            } else {
              line = paste(line,"\t\t\t<data key=\"targetArrowShape\">T</data>\n",sep="")
              line = paste(line,"\t\t\t<data key=\"label\">",  outputSummary.df[index,1], "--|", outputSummary.df[index,2], "</data>\n",sep="")
            }
        }
        line = paste(line,"\t\t\t<data key=\"edgeType\">2</data>\n",sep="")

  } else if(outputSummary.df[index, "infOrt"] == 6){
    line = paste(line,"\t\t<edge target=\"", outputSummary.df[index,2], "\" source=\"", outputSummary.df[index,1], "\" directed=\"true\">\n",sep="")

    if(is.na(outputSummary.df[index, "partial_correlation"])){
            line = paste(line,"\t\t\t<data key=\"sourceArrowShape\">arrow</data>\n",sep="")
            line = paste(line,"\t\t\t<data key=\"targetArrowShape\">arrow</data>\n",sep="")
            line = paste(line,"\t\t\t<data key=\"label\">",  outputSummary.df[index,2], "&lt;-&gt;", outputSummary.df[index,1], "</data>\n",sep="")
    } else {
        if(outputSummary.df[index, "partial_correlation"] > 0){
          line = paste(line,"\t\t\t<data key=\"sourceArrowShape\">arrow</data>\n",sep="")
          line = paste(line,"\t\t\t<data key=\"targetArrowShape\">arrow</data>\n",sep="")
          line = paste(line,"\t\t\t<data key=\"label\">",  outputSummary.df[index,2], "&lt;-&gt;", outputSummary.df[index,1], "</data>\n",sep="")
        } else {
          line = paste(line,"\t\t\t<data key=\"sourceArrowShape\">T</data>\n",sep="")
          line = paste(line,"\t\t\t<data key=\"targetArrowShape\">T</data>\n",sep="")
          line = paste(line,"\t\t\t<data key=\"label\">",  outputSummary.df[index,2], "|-|", outputSummary.df[index,1], "</data>\n",sep="")
        }
    }
    line = paste(line,"\t\t\t<data key=\"edgeType\">6</data>\n",sep="")
  }

  if(!all(is.na(outputSummary.df[, "log_confidence"]))) {
    if(outputSummary.df[index, "log_confidence"] <= 1)
      value = 1
    else if(outputSummary.df[index, "log_confidence"] >= 20)
      value = 8
    else 
      value = outputSummary.df[index, "log_confidence"] * 8 / 20
  } else {
      value = (abs(outputSummary.df[index, "partial_correlation"]) + 1) * 4
  }
  
  line = paste(line,"\t\t\t<data key=\"weight\">", value ,"</data>\n",sep="")

  line = paste(line,"\t\t\t<data key=\"upstream\">", outputSummary.df[index, "ai"] ,"</data>\n",sep="")
  line = paste(line,"\t\t\t<data key=\"info\">", outputSummary.df[index, "info"] ,"</data>\n",sep="")
  line = paste(line,"\t\t\t<data key=\"complexity\">", outputSummary.df[index, "cplx"] ,"</data>\n",sep="")
  line = paste(line,"\t\t\t<data key=\"nSamples\">", outputSummary.df[index, "Nxy_ui"] ,"</data>\n",sep="")
  if(confRatio != 1)
    line = paste(line,"\t\t\t<data key=\"confidenceRatio\">", outputSummary.df[index, "confidence_ratio"] ,"</data>\n",sep="")
  line = paste(line,"\t\t\t<data key=\"log_confidence\">", outputSummary.df[index, "log_confidence"] ,"</data>\n",sep="")
  line = paste(line,"\t\t\t<data key=\"sign\">", outputSummary.df[index, "sign"] ,"</data>\n",sep="")
  line = paste(line,"\t\t\t<data key=\"partialCorrelation\">", outputSummary.df[index, "partial_correlation"] ,"</data>\n",sep="")

  line = paste(line,"\t\t</edge>\n",sep="")

}
line = paste(line,"\t</graph>\n",sep="")

line = paste(line,"</graphml>\n",sep="")

writeLines(line, "network.graphml" )

############    XGMML
#if we do not have the layout we save a graphml file else a xgmml file with the node position
if(length(summaryArg.list[["argLayoutFile"]]) > 0){
    
    layout = read.table(summaryArg.list[["argLayoutFile"]], header=F, as.is=T, sep="\t", check.names = F)
    if(ncol(layout) == 2){
        xcol=1
        ycol=2
        rownames(layout) = colnames(inputData.df)
    } else {
        xcol=2
        ycol=3
        rownames(layout) = layout[,1]
    }
    
    line = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
    line = paste(line, "<graph label=\"graph\"",sep = "")

    line = paste(line, " xmlns:dc=\"http://purl.org/dc/elements/1.1/\"",sep = "")
    line = paste(line, " xmlns:xlink=\"http://www.w3.org/1999/xlink\"",sep = "")
    line = paste(line, " xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"",sep = "")
    line = paste(line, " xmlns:cy=\"http://www.cytoscape.org\"",sep = "")
    line = paste(line, " xmlns=\"http://www.cs.rpi.edu/XGMML\"",sep = "")

    line = paste(line, " directed=\"1\">\n",sep = "")
    #cicle on nodes
    for(node in colnames(inputData.df) ){
      line = paste(line,"\t\t<node label=\"", node ,"\" id=\"", node, "\">\n",sep="")
      line = paste(line,"\t\t\t<att name=\"size\" type=\"integer\" value=\"32\"/>\n",sep="")
      x = layout[node,xcol]*10
      y = -layout[node,ycol]*10
      line = paste(line,"\t\t\t<graphics fill=\"#f5f5f5\" x=\"",x,"\" y=\"",y,"\" cy:nodeLabelFont=\"Arial-0-11\" labelanchor=\"c\" type=\"ELLIPSE\" cy:nodeTransparency=\"0.8\" h=\"32\" width=\"1\" outline=\"#666666\" w=\"32\"/>\n",sep="")
      line = paste(line,"\t\t</node>\n",sep="")
    }


    line = paste(line,"\n",sep="")

    indexes = which(outputSummary.df["type"] == "P" | outputSummary.df["type"] == "TP" | outputSummary.df["type"] == "FP")

    #cycle on edges
    for(index in indexes) {
        sourceArrowNum = 0
        targetArrowNum = 0
      if(outputSummary.df[index, "infOrt"] == 1){
        line = paste(line,"\t\t<edge label=\"", outputSummary.df[index,2], "---", outputSummary.df[index,1], 
                     "\" target=\"", outputSummary.df[index,2], "\" source=\"", outputSummary.df[index,1], "\">\n",sep="")
        line = paste(line,"\t\t\t<att name=\"edgeType\" type=\"integer\" value=\"1\"/>\n",sep="")

       
        
      } else if(outputSummary.df[index, "infOrt"] == 2){

        if(is.na(outputSummary.df[index, "partial_correlation"])){
            value="arrow"
            varchar=intToUtf8(187)
            label= paste(outputSummary.df[index,1], "--&gt;", outputSummary.df[index,2], sep="")

        } else {
            if(outputSummary.df[index, "partial_correlation"] > 0){
              value="arrow"
              varchar=intToUtf8(187)
              label= paste(outputSummary.df[index,1], "--&gt;", outputSummary.df[index,2], sep="")
            } else {
              value="T"
              varchar="|"
              label= paste(outputSummary.df[index,1], "--|", outputSummary.df[index,2], sep="")
            }
        }
        line = paste(line,"\t\t<edge label=\"", label, 
                     "\" target=\"", outputSummary.df[index,2], "\" source=\"", outputSummary.df[index,1], "\">\n",sep="")
        line = paste(line,"\t\t\t<att name=\"targetArrowShape\" type=\"string\" value=\"", value ,"\"/>\n",sep="")
        line = paste(line,"\t\t\t<att name=\"sourceArrowShape\" type=\"string\" value=\"none\"/>\n",sep="")
        line = paste(line,"\t\t\t<att name=\"edgeType\" type=\"integer\" value=\"2\"/>\n",sep="")

        sourceArrowNum = 0
        targetArrowNum = fromStringToNumberArrowType(value)

      }  else if(outputSummary.df[index, "infOrt"] == -2){
        if(is.na(outputSummary.df[index, "partial_correlation"])){
            value="arrow"
            varchar=intToUtf8(187)
            label= paste(outputSummary.df[index,2], "--&gt;", outputSummary.df[index,1], sep="")
        } else {
            if(outputSummary.df[index, "partial_correlation"] > 0){
              value="arrow"
              varchar=intToUtf8(187)
              label= paste(outputSummary.df[index,2], "--&gt;", outputSummary.df[index,1], sep="")
            } else {
              value="T"
              varchar="|"
              label= paste(outputSummary.df[index,2], "--|", outputSummary.df[index,1], sep="")
            }
        }
        line = paste(line,"\t\t<edge label=\"", label, 
                     "\" target=\"", outputSummary.df[index,1], "\" source=\"", outputSummary.df[index,2], "\">\n",sep="")
        line = paste(line,"\t\t\t<att name=\"targetArrowShape\" type=\"string\" value=\"", value ,"\"/>\n",sep="")
        line = paste(line,"\t\t\t<att name=\"sourceArrowShape\" type=\"string\" value=\"none\"/>\n",sep="")
        line = paste(line,"\t\t\t<att name=\"edgeType\" type=\"integer\" value=\"2\"/>\n",sep="")
        
        sourceArrowNum = 0
        targetArrowNum = fromStringToNumberArrowType(value)

      } else if(outputSummary.df[index, "infOrt"] == 6){
        if(is.na(outputSummary.df[index, "partial_correlation"])){
            value="arrow"
            varchar=intToUtf8(187)
            label= paste(outputSummary.df[index,2], "&lt;-&gt;", outputSummary.df[index,1], sep="")
        } else {
            if(outputSummary.df[index, "partial_correlation"] > 0){
              value="arrow"
              varchar=intToUtf8(187)
              label= paste(outputSummary.df[index,2], "&lt;-&gt;", outputSummary.df[index,1], sep="")
            } else {
              value="T"
              varchar="|"
              label= paste(outputSummary.df[index,2], "|-|", outputSummary.df[index,1], sep="")
          }
        }
        line = paste(line,"\t\t<edge label=\"", label, 
                     "\" target=\"", outputSummary.df[index,1], "\" source=\"", outputSummary.df[index,2], "\">\n",sep="")
        line = paste(line,"\t\t\t<att name=\"targetArrowShape\" type=\"string\" value=\"", value, "\"/>\n",sep="")
        line = paste(line,"\t\t\t<att name=\"sourceArrowShape\" type=\"string\" value=\"", value, "\"/>\n",sep="")
        line = paste(line,"\t\t\t<att name=\"edgeType\" type=\"integer\" value=\"6\"/>\n",sep="")

        sourceArrowNum = fromStringToNumberArrowType(value)
        targetArrowNum = fromStringToNumberArrowType(value)
      }
      
      if(outputSummary.df[index, "log_confidence"] <= 1)
        value = 1
      else if(outputSummary.df[index, "log_confidence"] >= 20)
        value = 8
      else 
        value = outputSummary.df[index, "log_confidence"] * 8 / 20
      
      line = paste(line,"\t\t\t<att name=\"weight\" type=\"double\" value=\"", value ,"\"/>\n",sep="")
      line = paste(line,"\t\t\t<att name=\"upstream\" type=\"string\" value=\"", outputSummary.df[index, "ai"] ,"\"/>\n",sep="")
      line = paste(line,"\t\t\t<att name=\"info\" type=\"double\" value=\"", outputSummary.df[index, "info"] ,"\"/>\n",sep="")
      line = paste(line,"\t\t\t<att name=\"complexity\" type=\"double\" value=\"", outputSummary.df[index, "cplx"] ,"\"/>\n",sep="")
      line = paste(line,"\t\t\t<att name=\"nSamples\" type=\"integer\" value=\"", outputSummary.df[index, "Nxy_ai"] ,"\"/>\n",sep="")
      if(confRatio != 1)
        line = paste(line,"\t\t\t<att name=\"confidenceRatio\" type=\"double\" value=\"", outputSummary.df[index, "confidence_ratio"] ,"\"/>\n",sep="")
      line = paste(line,"\t\t\t<att name=\"log_confidence\" type=\"double\" value=\"", outputSummary.df[index, "log_confidence"] ,"\"/>\n",sep="")
      line = paste(line,"\t\t\t<att name=\"sign\" type=\"string\" value=\"", outputSummary.df[index, "sign"] ,"\"/>\n",sep="")
      line = paste(line,"\t\t\t<att name=\"partialCorrelation\" type=\"double\" value=\"", outputSummary.df[index, "partial_correlation"] ,"\"/>\n",sep="")
      if(confRatio != 1)
        line = paste(line,"\t\t\t<att name=\"isCut\" type=\"double\" value=\"", outputSummary.df[index, "isCut"] ,"\"/>\n",sep="")

      if(outputSummary.df[index, "sign"] == "+")
        fillColor = "#ff3300"
      else if(outputSummary.df[index, "sign"] == "-")
        fillColor = "#aaaaff"
      else
        fillColor = "#808080"
      line = paste(line,"\t\t\t<graphics cy:sourceArrowColor=\"#000000\" cy:targetArrowColor=\"#000000\" width=\"", value, 
                "\" cy:edgeLineType=\"SOLID\" cy:targetArrow=\"",targetArrowNum,"\" cy:sourceArrow=\"",sourceArrowNum,"\" fill=\"", fillColor,"\"/>\n", sep="")



      
      line = paste(line,"\t\t</edge>\n",sep="")
    }
    line = paste(line,"\t</graph>\n",sep="")

    # name = basename(file_path_sans_ext(outDirPath))

    writeLines(line, paste("network.xgmml" ,sep=""))
}


##################################### NETWORK IN GRAPHML END

rm(inputData.df)

cat( "\t# --------\n# -> END Summary...\n" )

