#!/usr/bin/Rscript
#### Delete all variables
rm( list = ls ( all = TRUE ) )
library(tools)
getos <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
userDir = getwd()
isTest = FALSE

osVersion = getos()
slash = "/"
if(osVersion == "windows")
  slash = "\\"

##### Load libraries
suppressMessages(require(getopt))             # Parse command line options

full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
                       error=function(e) # works when using R CMD
                         normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
rootDir = sub("miic.R", "", full.fpath, fixed = TRUE)

setwd( rootDir )

### FOR TEST ####
if( isTest == TRUE ){ setwd( rootDir ) }
### FOR TEST ####

#### Load supplementary functions
source( file.path( "lib", "miic.utils.lib.R" ) )

### FOR TEST ####
if( isTest == TRUE )
{
  discoNetArg.list.external = list()
  discoNetArg.list.external[["argMethod"]] = "miic";
  discoNetArg.list.external[["argSteps"]] = c(1,5,3,4);
  discoNetArg.list.external[["argInData"]] = "";
  discoNetArg.list.external[["argOutDir"]] = "";
  discoNetArg.list.external[["argBlackBox"]] =  ""
  discoNetArg.list.external[["argEffN"]] = -1;
  discoNetArg.list.external[["argCplx"]] = "nml";
  discoNetArg.list.external[["argEta"]] = 1;
  discoNetArg.list.external[["argDegeneracy"]] = FALSE;
  discoNetArg.list.external[["argLatent"]] = "no";
  discoNetArg.list.external[["argHvs"]] = 0;
  discoNetArg.list.external[["argScore"]] = "bic";
  discoNetArg.list.external[["argRestart"]] = 10;

  discoNetArg.list.external[["argTrueEdgesFile"]] = "";
  discoNetArg.list.external[["argStateOrder"]] = "";
  discoNetArg.list.external[["argLayout"]] = "";

}
### FOR TEST ####
stateOrderAdded = FALSE
cat( "# --------\n# -> START miic...\n" )

#### Parse command line arguments
#### ----
discoNetArg.list = list();
if( !exists( "discoNetArg.list.external" ) )
{
  cat(userDir)
  discoNetArg.list = discoNet.parseCommandLine(userDir, rootDir)
  #check if the output directory is correct

  if(slash == "\\")
    regexp = "[\\\\]"
  else
    regexp = slash

  if(unlist(strsplit(discoNetArg.list[["argOutDir"]], ""))[length(unlist(strsplit(discoNetArg.list[["argOutDir"]], "")))] == slash)
    discoNetArg.list[["argOutDir"]] = paste(unlist(strsplit(discoNetArg.list[["argOutDir"]], ""))[1:(length(unlist(strsplit(discoNetArg.list[["argOutDir"]], "")))-1)], collapse="")

  tmp = paste(unlist(strsplit(discoNetArg.list[["argOutDir"]], regexp))[1:(length(unlist(strsplit(discoNetArg.list[["argOutDir"]], regexp)))-1)], collapse=slash)
  folderName = paste(unlist(strsplit(discoNetArg.list[["argOutDir"]], regexp))[length(unlist(strsplit(discoNetArg.list[["argOutDir"]], regexp)))], collapse="")

  if(!file.exists(tmp) | (length(unlist(strsplit(discoNetArg.list[["argOutDir"]], regexp))) == 1)){
    discoNetArg.list[["argOutDir"]] = paste(userDir,folderName, sep = regexp)
    cat("\t# --W--> Output folder will be created in\n\t",userDir,"\n")
  } else {
  discoNetArg.list[["argOutDir"]] = paste(file_path_as_absolute(tmp),folderName, sep = slash)
  }
} else {

  discoNetArg.list = discoNetArg.list.external
}

# Check that the input file exists
if(!file.exists(discoNetArg.list[["argInData"]])){
  cat("# -- EXIT miic -- The input file does not exist!\n")
  quit(save = "no", status = 0)
}

# Check that the method can be processed
if( !( discoNetArg.list[["argMethod"]] %in% c( "miic") ) ){
  # stop( "# Unkown method: ", discoNetArg.list[["argMethod"]] )
  cat("# -- EXIT miic -- Unkown method:", discoNetArg.list[["argMethod"]], "\n")
  quit(save = "no", status = 0)
}

# When searching for the best eta (ie shuffle>0), only the skeleton step should be performed
if( discoNetArg.list[["argShuffle"]] > 0 ){ discoNetArg.list[["argSteps"]] = c(1) }

# 'generic' name of the ouput files
inferredEdgesFileName = "edgesList.txt"
adjMatFileName = "adjacencyMatrix.txt"
summaryFileName = inferredEdgesFileName

currCmd = ""
# Learn the skeleton
# --------
# ---- miic ----




###################################################################  check datasets

if(!file.exists(discoNetArg.list[["argInData"]])){errCode="110";stop(errorCodeToString(err_code))}
data = try(as.matrix(read.table(discoNetArg.list[["argInData"]], header = FALSE, stringsAsFactors = FALSE, sep = "\t", check.names = F, na.strings = c('', NA) )))
if(class(data) == "try-error"){errCode = "121";stop(errorCodeToString(err_code))}

inputData = read.table(discoNetArg.list[["argInData"]], header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = F, na.strings = c('', NA) )


err_code = checkInput(discoNetArg.list[["argInData"]], discoNetArg.list[["argMethod"]], inputData)
if(err_code != "0"){
  stop(errorCodeToString(err_code))
}

  


out =  paste('"',discoNetArg.list[["argOutDir"]],'"',sep="")
modeOrient=NA
if( 5 %in% discoNetArg.list[["argSteps"]] )
  modeOrient=5
if( 2 %in% discoNetArg.list[["argSteps"]] )
  modeOrient=2

if( discoNetArg.list[["argMethod"]] == "miic" )
{
  executable = "./miic"

  setwd( file.path( "..", discoNetArg.list[["argMethod"]] ) )

  if( dir.exists(discoNetArg.list[["argOutDir"]]) & (1 %in% discoNetArg.list[["argSteps"]]) ) { stop( "# The output directory already exists." ) }
  if(!dir.exists(discoNetArg.list[["argOutDir"]]))
    dir.create(discoNetArg.list[["argOutDir"]])

  if(length(discoNetArg.list[["argCutShf"]]) > 0 & length(discoNetArg.list[["argConfidenceRatio"]]) > 0 ){
    out_shf_dir = file.path( discoNetArg.list[["argOutDir"]],paste ('shuffle',discoNetArg.list[["argCutShf"]], sep="_"),paste(
                    "filtered_network",discoNetArg.list[["argConfidenceRatio"]], sep='_' ) ,fsep = slash )

    if(!dir.exists(out_shf_dir))
      dir.create( out_shf_dir , recursive=T) 
    else
      stop("The folder of the confidence cut already exists.")
  }
 
  ## START CALLS

  if( 1 %in% discoNetArg.list[["argSteps"]] )
  {
    if( 2 %in% discoNetArg.list[["argSteps"]] | 5 %in% discoNetArg.list[["argSteps"]]){
      currCmd = paste( executable, " -i ", discoNetArg.list[["argInData"]], " -o ", out
                       , " -d 1,", modeOrient, " -n ", discoNetArg.list[["argEffN"]], " -c ",
                       discoNetArg.list[["argCplx"]], sep = '' )
      if( discoNetArg.list[["argPropag"]] == FALSE ){ currCmd = paste( currCmd, ' -p 0', sep='' ) }
    }
    else {
      currCmd = paste( executable, " -i ", discoNetArg.list[["argInData"]], " -o ", out
                       , " -d 1,", " -n ", discoNetArg.list[["argEffN"]], " -c ", discoNetArg.list[["argCplx"]], sep = '' )
    }

    if( length(discoNetArg.list[["argBlackBox"]]) > 0  ){ currCmd = paste( currCmd, ' -b ', discoNetArg.list[["argBlackBox"]], sep='' ) }
    if( discoNetArg.list[["argLatent"]] == TRUE ){ currCmd = paste( currCmd, ' -l', sep='' ) }
    if( discoNetArg.list[["argTplReuse"]] == FALSE ){ currCmd = paste( currCmd, ' -r 0', sep='' ) }
    if( discoNetArg.list[["argK23"]] == FALSE ){ currCmd = paste( currCmd, ' -k 0', sep='' ) }
    if( discoNetArg.list[["argDegeneracy"]] == TRUE ){ currCmd = paste( currCmd, ' -g', sep='' ) }
    if( discoNetArg.list[["argNoInitEta"]] == TRUE ){ currCmd = paste( currCmd, ' -f', sep='' ) }
    if( discoNetArg.list[["argThreads"]] > 1) {currCmd = paste( currCmd, ' -z ', discoNetArg.list[["argThreads"]], sep='') }
    if( length(discoNetArg.list[["argWeights"]]) > 0  ){ currCmd = paste( currCmd, ' -q ', discoNetArg.list[["argWeights"]], sep='' ) }

    if( length(discoNetArg.list[["argCutShf"]]) > 0 & length(discoNetArg.list[["argConfidenceRatio"]]) > 0 ){
      currCmd = paste( currCmd, ' -e ', discoNetArg.list[["argConfidenceRatio"]], sep='')
      currCmd = paste( currCmd, ' -s ', discoNetArg.list[["argCutShf"]], sep='')
    }


    currCmd = paste( currCmd, ' -x ', discoNetArg.list[["argSeed"]], sep='' )

    # Recall parameters
    cat('\t# --------\n')

    cat('# -> START Skeleton & Orientation Learning ...\n')
    cat('\t# Input Data file --> ', discoNetArg.list[["argInData"]],'\n')
    cat('\t# Output directory --> ', discoNetArg.list[["argOutDir"]],'\n')

    #print(currCmd)

    # Call the script
    system( currCmd)#, intern=TRUE )
    cat('# -> END Skeleton & Orientation Learning\n')


  }
  # this is only for the confidence tool execution
  else
  {

    if(!is.na(modeOrient)){
      if( 2 %in% discoNetArg.list[["argSteps"]] ){
        currCmd = paste( executable, " -i ", discoNetArg.list[["argInData"]], " -o ", out
                         , " -d 2", " -n ", discoNetArg.list[["argEffN"]], " -c ", discoNetArg.list[["argCplx"]], sep = '' )
      } else if( 5 %in% discoNetArg.list[["argSteps"]] ){
        currCmd = paste( executable, " -i ", discoNetArg.list[["argInData"]], " -o ", out
                         , " -d 5", " -n ", discoNetArg.list[["argEffN"]], " -c ", discoNetArg.list[["argCplx"]]
                         ,," -a ", discoNetArg.list[["argHvs"]], sep = '' )
      }

      # Build the command line
      if( 2 %in% discoNetArg.list[["argSteps"]] | 5 %in% discoNetArg.list[["argSteps"]]){
        currCmd = paste( currCmd, ' -x ', discoNetArg.list[["argSeed"]], sep='' )

        if( discoNetArg.list[["argLatent"]] == TRUE ){ currCmd = paste( currCmd, ' -l ', sep='' ) }
        if( discoNetArg.list[["argK23"]] == FALSE ){ currCmd = paste( currCmd, ' -k 0', sep='' ) }
        if( discoNetArg.list[["argDegeneracy"]] == TRUE ){ currCmd = paste( currCmd, ' -g ', sep='' ) }
        if( discoNetArg.list[["argPropag"]] == FALSE ){ currCmd = paste( currCmd, ' -p 0', sep='' ) }
        if( discoNetArg.list[["argThreads"]] > 1) {currCmd = paste( currCmd, ' -z ', discoNetArg.list[["argThreads"]], sep='') }
        if( length(discoNetArg.list[["argWeights"]]) > 0  ){ currCmd = paste( currCmd, ' -q ', discoNetArg.list[["argWeights"]], sep='' ) }


        if( length(discoNetArg.list[["argCutShf"]]) > 0 & length(discoNetArg.list[["argConfidenceRatio"]]) > 0 ){
          currCmd = paste( currCmd, ' -e ', discoNetArg.list[["argConfidenceRatio"]], sep='')
          currCmd = paste( currCmd, ' -s ', discoNetArg.list[["argCutShf"]], sep='')
        }

        #print(currCmd)

        # Call the script
        system( currCmd, intern=TRUE )
      }
    }
  }
}

setwd( rootDir )

# ----- End Other Techniques ----

# Summarize the results
# --------
if( 3 %in% discoNetArg.list[["argSteps"]] )
{
  if( length( discoNetArg.list[["argTrueEdgesFile"]] ) > 0 ){
    err_code = checkTrueEdges(discoNetArg.list[["argTrueEdgesFile"]], inputData)
    if(err_code != "0"){
      print(errorCodeToString(err_code))
      print("WARNING: True edges file will be ignored!")
      discoNetArg.list[["argTrueEdgesFile"]] = character(0)
    }
  }

  #STATE ORDER FILE
  if( length( discoNetArg.list[["argStateOrder"]] ) > 0 ){
    err_code = checkStateOrder(discoNetArg.list[["argStateOrder"]], inputData)
    if(err_code != "0"){
      print(errorCodeToString(err_code))
      print("WARNING: Cathegory order file will be ignored!")
      discoNetArg.list[["argStateOrder"]] = character(0)
    }
  }


  if( length(discoNetArg.list[["argCutShf"]]) > 0 & length(discoNetArg.list[["argConfidenceRatio"]]) > 0 ){

    vec = c(discoNetArg.list[["argOutDir"]], paste(discoNetArg.list[["argOutDir"]], slash, "shuffle_", discoNetArg.list[["argCutShf"]], slash,
          "filtered_network_", discoNetArg.list[["argConfidenceRatio"]], sep=""))
  } else {
    vec = discoNetArg.list[["argOutDir"]]
  }
   
  for(val in vec){

    out =  paste('"',val,'"',sep="")

    tmp.inferredEdgesFileName = gsub( ".txt", paste( '', discoNetArg.list[["argMethod"]], "txt", sep='.' ), inferredEdgesFileName )

    #### If the orientation with probabilities is required, choose the outputs .orientProba.
    if( 2 %in% discoNetArg.list[["argSteps"]] )
    {
      tmp.adjMatFileName = gsub( ".txt", paste( '', discoNetArg.list[["argMethod"]], "orientProba","txt", sep='.' ), adjMatFileName )
    } else if ( 5 %in% discoNetArg.list[["argSteps"]] ){
      tmp.adjMatFileName = gsub( ".txt", paste( '', discoNetArg.list[["argMethod"]], "orient","txt", sep='.' ), adjMatFileName )
    } else {
      tmp.adjMatFileName = gsub( ".txt", paste( '', discoNetArg.list[["argMethod"]], "orientProba.txt", sep='.' ), adjMatFileName )
    }
    

    if(!file.exists(file.path(val,tmp.inferredEdgesFileName)) | !file.exists(file.path(val,tmp.adjMatFileName))){
      print(file.path(val,tmp.inferredEdgesFileName))
      print(file.path(val,tmp.adjMatFileName))
      stop("ERROR: summary can't be called on not existing network")
    }

    # Build the command line
    currCmd = paste( "Rscript gmSummary.R -o ", out , " -d ", discoNetArg.list[["argInData"]], " -i ", tmp.inferredEdgesFileName, " -a ", tmp.adjMatFileName, sep = '' )

    if( length( discoNetArg.list[["argTrueEdgesFile"]] ) > 0 )
    { currCmd = paste( currCmd, ' -t ', discoNetArg.list[["argTrueEdgesFile"]], sep='' ) }
    if( length( discoNetArg.list[["argLayout"]] ) > 0 )
    { currCmd = paste( currCmd, ' -l ', discoNetArg.list[["argLayout"]], sep='' ) }
    if( length( discoNetArg.list[["argStateOrder"]] ) > 0 )
    { currCmd = paste( currCmd, ' -s ', discoNetArg.list[["argStateOrder"]], sep='' ) }

    if(discoNetArg.list[["argWriteNetwork"]])
      currCmd = paste( currCmd, ' -n ', sep='' )
    print(paste("COMMAND SUMMARY:" , currCmd))
    system( currCmd )
  }


  if( length(discoNetArg.list[["argCutShf"]]) > 0 & length(discoNetArg.list[["argConfidenceRatio"]]) > 0 ){

    tmp_argOutDir = paste(discoNetArg.list[["argOutDir"]], slash ,paste("shuffle",discoNetArg.list[["argCutShf"]],sep="_")
                                 , slash ,paste("filtered_network",discoNetArg.list[["argConfidenceRatio"]],sep="_"), sep = "")
    tmp_sumFile = file.path(tmp_argOutDir, "edgesList.miic.summary.txt")
    tmp_pvalFile = file.path(tmp_argOutDir, "confRatios.txt")
    tmp_sum = read.table(tmp_sumFile, header=T, as.is=T, sep="\t", check.names = F)
    conf_col = rep(1, nrow( tmp_sum ) )
    isCut = rep(NA, nrow(tmp_sum))
    tmp_sum = cbind(tmp_sum,conf_col,isCut)

    tmp_sum = cbind(tmp_sum,conf_col)
    colnames(tmp_sum) = c('x','y','type','ai','info','cplx','Nxy_ai','log_confidence','infOrt','trueOrt', 'isOrtOk', 'sign','partial_correlation','confidence_ratio')
    tmp_pval = read.table(tmp_pvalFile, header=T, as.is=T, sep="\t", check.names = F)
    tmp_pval[,"confidence_ratio"] = as.numeric(tmp_pval[,"confidence_ratio"])
    #tmp_pval = tmp_pval[which(tmp_pval[,3] <= discoNetArg.list[["argConfidenceRatio"]]), 3]
    for(r in 1:nrow(tmp_pval))
    {
      tmp_sum[which(tmp_sum[,"x"] == tmp_pval[r,"x"] & tmp_sum[,"y"] == tmp_pval[r,"y"]),'confidence_ratio'] =tmp_pval[r,"confidence_ratio"]
      if(tmp_pval[r,"confidence_ratio"] < discoNetArg.list[["argConfidenceRatio"]]){
        tmp_sum[which(tmp_sum[,"x"]==tmp_pval[r,"x"] & tmp_sum[,"y"]==tmp_pval[r,"y"]),'isCut']='N'
      }
      else{
        tmp_sum[which(tmp_sum[,"x"]==tmp_pval[r,"x"] & tmp_sum[,"y"]==tmp_pval[r,"y"]),'isCut']='Y'
      }
    }
    tmp_sum = tmp_sum[,c('x','y','type','ai','info','cplx','Nxy_ai','log_confidence','confidence_ratio','infOrt','trueOrt', 'isOrtOk', 'sign','partial_correlation','isCut')]
    write.table(tmp_sum,tmp_sumFile, col.names=T, row.names=T, sep='\t',quote=F)
    rm(tmp_sum); rm(tmp_pval)
  }
}

# Plot the graph
# --------
if( 4 %in% discoNetArg.list[["argSteps"]] )
{

  #LAYOUT FILE
  if( length( discoNetArg.list[["argLayout"]] ) > 0 ){
    err_code = checkLayout(discoNetArg.list[["argLayout"]])
    if(err_code != "0"){
     print(errorCodeToString(err_code))
      print("WARNING: Layout file will be ignored!")
      discoNetArg.list[["argLayout"]] = character(0)
    }
  }

  if( length(discoNetArg.list[["argCutShf"]]) > 0 & length(discoNetArg.list[["argConfidenceRatio"]]) > 0 ){

    vec = c(discoNetArg.list[["argOutDir"]], paste(discoNetArg.list[["argOutDir"]], slash, "shuffle_", discoNetArg.list[["argCutShf"]], slash,
          "filtered_network_", discoNetArg.list[["argConfidenceRatio"]], sep=""))
  } else {
    vec = discoNetArg.list[["argOutDir"]]
  }

  for(val in vec){

    out =  paste('"',val,'"',sep="")


    # Update the name of the ouput files from previous steps

    tmp.summaryFileName = gsub( ".txt", paste( '', discoNetArg.list[["argMethod"]],"summary","txt", sep='.' ), summaryFileName )

    if(!file.exists(file.path(discoNetArg.list[["argOutDir"]],tmp.summaryFileName)))
      stop("ERROR: plot can't be called on not existing network")



    # Build the command line
    currCmd = paste( "Rscript gmPlot.R -i ", discoNetArg.list[["argInData"]], " -o ", out, " -s ", tmp.summaryFileName, sep = '' )

    if( length( discoNetArg.list[["argLayout"]] ) > 0 ){ currCmd = paste( currCmd, ' -l ', discoNetArg.list[["argLayout"]], sep='' ) }
    if( discoNetArg.list[["argLatent"]] == TRUE ){ currCmd = paste( currCmd, ' -h ', sep='' ) }
    #print(paste("COMMAND PLOT:" , currCmd))
    # Call the script
    system( currCmd )
  }
}
cleanOutput(discoNetArg.list[["argOutDir"]])


## Plot the edges strength in histogramm
edgeStr_file = file.path( discoNetArg.list[["argOutDir"]], "edgesScores.txt" )

if(file.exists(edgeStr_file))
{
  file.remove(edgeStr_file)
}