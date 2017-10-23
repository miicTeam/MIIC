#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <algorithm>
#include <ctime>
#include <unistd.h>
#include <string.h>

#include "structure.h"
#include "utilities.h"
#include "log.h"
#include "skeletonInitialization.h"
#include "skeletonIteration.h"


using namespace std;

bool skeleton(Environment& environment, string slash, clock_t startTimeWhole, int argc, char* argv[]){

	std::stringstream ss;
	ss.str("");
	ss << environment.outDir << slash << "log.miic.skeleton.txt";

	char *addr = new char[ss.str().length()+1];
	strcpy(addr,ss.str().c_str());
	Log* pLog = new Log(addr);


	int** matrixOriginal = NULL; // used to save original values 
	clock_t startTime;
	stringstream log_shuffle;
	log_shuffle.str("");

	ExecutionTime execTime;

	string outputFilePath;
	string filePath;

	//// Define a log file
	if( environment.isVerbose == true ){cout << "# --------\n# Define a log file...\n";}


	// set the environment
	setEnvironment(environment);

	if(environment.blackbox_name.compare("") != 0)
		readBlackbox(environment, slash);
						

	// print the environment
 	printEnvironment(environment, pLog);
 	
	startTime = std::clock();


	// write progress
	ss.str("");
	ss << environment.outDir << slash << "progress.txt";
	ofstream output;
	string filename = ss.str();
	output.open(filename.c_str());
	output << "Skeleton initialization";
	output.close();
	
	skeletonInitialization(environment);

	// ----
	long double spentTime = (std::clock() - startTime) / (double)(CLOCKS_PER_SEC / 1000) /1000;
	environment.execTime.init = spentTime;
	if( environment.isVerbose == true ){ cout << "\n# ----> First contributing node elapsed time:" << spentTime << "sec\n\n"; }
	

	if( environment.numNoMore == 0 && environment.numSearchMore == 0 ) { 
		if( environment.isVerbose == true ){ cout << "# ------| Only phantom edges found.\n"; }
	} else if( environment.numSearchMore > 0 ) {
	
		//// Search for other Contributing node(s) (possible only for the edges still in 'searchMore', ie. 2)
		if( environment.isVerbose == true ){ cout << "\n# ---- Other Contributing node(s) ----\n\n"; }
		startTime = std::clock();
			
		stringstream ss;
		ss.str("");
		//// Save the execTime
		ss << environment.outDir << slash << "progress.txt";
		ofstream output;
		string filename = ss.str();
		output.open(filename.c_str());
		output << "Skeleton iteration";
		output.close();
		skeletonIteration(environment);
		long double spentTime = (std::clock() - startTime) / (double)(CLOCKS_PER_SEC / 1000) /1000;
		environment.execTime.iter = spentTime;
		environment.execTime.initIter = environment.execTime.init + environment.execTime.iter;
	}



	//// Save the ordered edges and make an adjacency matrix
	if(environment.isVerbose){ cout << "\n# ---- Save the ordered edges and make an adjacency matrix ----\n\n"; }

	if( environment.numNoMore > 0 )
	{
		ss.str("");
		//// Make a table from the list
		ss << environment.outDir << slash << "edgesList.miic.txt";
		saveEdgesListAsTable(environment, ss.str());

		ss.str("");

		//// Make an adjacency matrix from the list
		ss << environment.outDir << slash << "adjacencyMatrix.miic.txt";
		saveAdjMatrix(environment, ss.str());
		//adj.mat = saveEdgesListAsAdjMat( environment, myAllProp = colnames(gV$data), myOutputFilePath = outputFilePath )
	}	

	return 0;
}
