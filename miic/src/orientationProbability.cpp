#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <algorithm>
#include <ctime>
#include <unistd.h>
#include <string.h>
#include <cmath>

#include "structure.h"
#include "utilities.h"
#include "computeEnsInformation.h"


#include "probaOrientation_interface.h"


using namespace std;

int sign(double val){
	if(val < 0) 
		return -1;
	else if(val > 0) 
		return 1;
	else 
		return 0;
}

bool getSStructure(Environment& environment, const int posX, const int posY, const int posZ, bool isVerbose, vector< vector<int> >& allTpl, vector<double>& allI3){
	bool structFound = false;	
	//// To check if xk belongs to the {ui} of the base
	bool isXkInUi = false;

	vector<int> u(environment.edges[posX][posZ].edgeStructure->ui_vect_idx);
	if(environment.edges[posX][posZ].edgeStructure->ui_vect_idx.size() > 0){
		//// Remove xk from the {ui} if needed
		for(int i = 0; i < environment.edges[posX][posZ].edgeStructure->ui_vect_idx.size(); i++){
			if(u[i] == posY){
				u.erase(u.begin() + i);
				isXkInUi = true;
				break;
			}
		}
	}

	int* ui;
	int* zi = NULL;

	if(u.empty())
		ui = NULL;
	else
		ui = &u[0];

	vector<int> z;
	z.clear();
	z.push_back(posY);

	zi = &z[0];

	double* res = NULL;
	double Is = -1;
	double Cs = -1;
	res = computeEnsInformationNew(environment, ui, u.size(), zi, z.size(), -1, posX, posZ, environment.cplx);

	Is = res[7];
	Cs = res[8];

	if(environment.isK23){

		if(environment.isDegeneracy){
				Cs += log(3);
		}

		Is = Is + Cs;
	} 

	delete(res);
	
	

	vector<int> v;
	v.push_back(posX+1);
	v.push_back(posY+1);
	v.push_back(posZ+1);
	allTpl.push_back(v);
	allI3.push_back(Is);
}


bool orientationProbability(Environment environment, string slash, bool isVerbose){
	cout << "# -> START miic PROBA Orient...\n";
	std::stringstream ss;
	ss.str("");
	ss << environment.outDir << slash << "log.miic.orientProba.txt";

	ExecutionTime execTime;

	vector< vector<int> > allTpl;
	vector<double> allI3;
	double* ptrRetProbValues;

	//// Define a log file
	char *addr = new char[ss.str().length()+1];
	strcpy(addr,ss.str().c_str());
	Log* pLog = new Log(addr);

	delete [] addr;

	// print the environment
	printEnvironment(environment, pLog);

	////// GET ALL TPL that could be V/NON-V-STRUCTURES #######
	if(isVerbose){ cout << "\n# ---- Get all the unshielded triplets ----\n";}
	clock_t startTime_algo = std::clock();

	for(int pos = 0; pos < environment.noMoreAddress.size(); pos++){
		int posX = environment.noMoreAddress[pos]->i;
		int posY = environment.noMoreAddress[pos]->j;

		if(isVerbose){ cout << "\n# (" << pos << ") ---- Edge" << environment.nodes[posX].name << "--" << environment.nodes[posY].name << "\n" ; }

		//// Prepare a list that will contain the neighbors of "x" and the neighbors of "y"
		vector<int> neighboursX;
		vector<int> neighboursY;

		for(int i = pos + 1; i < environment.noMoreAddress.size(); i++){
			int posX1 = environment.noMoreAddress[i]->i;
			int posY1 = environment.noMoreAddress[i]->j;

			if(posY1 == posX && !environment.edges[posY][posX1].isConnected)
				neighboursX.push_back(posX1);
			else if(posY1 == posY && !environment.edges[posX][posX1].isConnected)
				neighboursY.push_back(posX1);
		}

		for(int i = pos + 1; i < environment.noMoreAddress.size(); i++){
			int posX1 = environment.noMoreAddress[i]->i;
			int posY1 = environment.noMoreAddress[i]->j;
			
			if(posX1 == posX && !environment.edges[posY][posY1].isConnected)
				neighboursX.push_back(posY1);
			else if(posX1 == posY && !environment.edges[posX][posY1].isConnected)
				neighboursY.push_back(posY1);
		}


		int sizeX = neighboursX.size();
		int sizeY = neighboursY.size();

		if(sizeX == 0 && sizeY == 0)
			continue;


		for(int i = 0; i < sizeX; i++){
			if(isVerbose){ cout << "\n# (" << pos << ") -----| Test tpl (" << environment.nodes[posY].name << ", " << environment.nodes[posX].name << ", " << environment.nodes[neighboursX[i]].name << ")\n" ; } 
			//// Get the structure if any
		   	getSStructure(environment, posY, posX, neighboursX[i], isVerbose, allTpl, allI3);
		}

		// iterate on neighbours of y
		for(int i = 0; i < sizeY; i++){
			if(isVerbose){ cout << "\n# (" << pos << ") -----| Test tpl (" << environment.nodes[posX].name << ", " << environment.nodes[posY].name << ", " << environment.nodes[neighboursY[i]].name << ")\n" ; } 
			//// Get the structure if any
			getSStructure(environment, posX, posY, neighboursY[i], isVerbose,  allTpl, allI3);
		}

	}


	//create the one line matrix
	int* oneLineMatrixallTpl = new int[allTpl.size()*3];
	for(int i = 0; i < allTpl.size();i++){
		for(int j = 0; j < 3;j++){
			// cout << j * environment.numSamples + i << " ";
			oneLineMatrixallTpl[j * allTpl.size() + i] = allTpl[i][j];
			allTpl[i][j]--;
		}
	}

	//// Compute the arrowhead probability of each edge endpoint
	int myNbrTpl = allTpl.size();
	if(myNbrTpl > 0){
		int propag = 0;
		if(environment.isPropagation)
			propag = 1;
		int degeneracy = 0;
		if(environment.isDegeneracy)
			degeneracy = 1;
		int latent = 0;
		if(environment.isLatent)
			latent = 1;

		ptrRetProbValues = getOrientTplLVDegPropag(myNbrTpl, oneLineMatrixallTpl, &allI3[0], latent, degeneracy, propag, environment.halfVStructures);

		// write results to file 
		stringstream ss;
		ss.str("");
	    ss << environment.outDir << slash << "orientations.probaArrowhead.txt";
		ofstream output;
		output.open(ss.str().c_str());

		output << "source1" << "\t" << "p1" << "\t" << "p2" << "\t" << "target" << "\t" << "p3" << "\t" << "p4" << "\t" << "source2" << "\t" << "NI3" << "\t" <<
		"Error" << "\n";

		for (int i=0; i < allTpl.size(); i++){
			int error = 0;
			int info_sign = allI3[i];
			int proba_sign = 0;

			if(ptrRetProbValues[i + (1 * allTpl.size())] > 0.5 && ptrRetProbValues[i + (2 * allTpl.size())] > 0.5)
				proba_sign = -1;
			else
				proba_sign = 1;

			if((sign(info_sign) != proba_sign && info_sign != 0) || (info_sign == 0 && proba_sign == -1))
				error = 1;

			output << environment.nodes[allTpl[i][0]].name << "\t" << ptrRetProbValues[i + (0 * allTpl.size())] << "\t" << ptrRetProbValues[i + (1 * allTpl.size())] << "\t" <<
			environment.nodes[allTpl[i][1]].name << "\t" << ptrRetProbValues[i + (2 * allTpl.size())] << "\t" << ptrRetProbValues[i + (3 * allTpl.size())] << "\t" <<
			environment.nodes[allTpl[i][2]].name << "\t" << allI3[i] << "\t" << error << "\n";
		}
	}

	long double spentTime = std::clock() - startTime_algo;
	if(isVerbose){ cout << "\n# ----> Propagate Orientations elapsed time:", spentTime, "sec.\n\n" ; }

	// ALGO STOP TIME
	long double spentTime_algo = (std::clock() - startTime_algo) / (double)(CLOCKS_PER_SEC / 1000) /1000;

	//Save the time spent on the algo (if the time for the skeleton already exists, copy the file, open and add before saving)
	ss.str("");
	//// Save the execTime
	ss << environment.outDir << slash << "execTime.miic.skeleton.txt";

	if(existsTest(ss.str())){
		readTime(ss.str(), execTime);
	     execTime.initIter += spentTime_algo;
	     execTime.initIterSave += spentTime_algo;
	 } else{ 
	    execTime.initIter = spentTime_algo;
	}
	
	ss.str("");
	//// Save the execTime
	ss << environment.outDir << slash << "execTime.miic.orient.skeleton.txt";
	saveExecTime(environment, execTime, ss.str());
	vector<double> myProba;
	//// UPDATE ADJ MATRIX #######
	if( myNbrTpl > 0 ){
		for(int i = 0; i < allTpl.size(); i++){
			myProba.clear();;
			myProba.push_back(ptrRetProbValues[i + (0 * allTpl.size())]);
			myProba.push_back(ptrRetProbValues[i + (1 * allTpl.size())]);

			if( (*max_element(myProba.begin(), myProba.end() ) - 0.5) > 0 ) // If there is AT LEAST ONE arrowhead
			{
	            if( (*min_element( myProba.begin(), myProba.end()  ) - 0.5) <= 0 ) // If there is ONLY ONE arrowhead
	            {

	            	double a = myProba[0] - *max_element(myProba.begin(), myProba.end() );
	                if( a == 0 )  // If p1 is max: n1 <-- n2
	                {
	                	environment.edges[allTpl[i][0]][allTpl[i][1]].isConnected = -2;
	                	environment.edges[allTpl[i][1]][allTpl[i][0]].isConnected = 2;
	                } else {
	                	environment.edges[allTpl[i][0]][allTpl[i][1]].isConnected = 2;
	                	environment.edges[allTpl[i][1]][allTpl[i][0]].isConnected = -2;
	                }
				} else {
					environment.edges[allTpl[i][0]][allTpl[i][1]].isConnected = 6;
                	environment.edges[allTpl[i][1]][allTpl[i][0]].isConnected = 6;
				}
			}

			myProba.clear();;
			myProba.push_back(ptrRetProbValues[i + (2 * allTpl.size())]);
			myProba.push_back(ptrRetProbValues[i + (3 * allTpl.size())]);

			if( (*max_element( myProba.begin(), myProba.end()  ) - 0.5) > 0 ) // If there is AT LEAST ONE arrowhead
			{
	            if( (*min_element( myProba.begin(), myProba.end()  ) - 0.5) <= 0 ) // If there is ONLY ONE arrowhead
	            {
	            	double a = myProba[0] - *max_element(myProba.begin(), myProba.end() );
	            	if( a == 0 )  // If p1 is max: n1 <-- n2
	                {
	                	environment.edges[allTpl[i][1]][allTpl[i][2]].isConnected = -2;
	                	environment.edges[allTpl[i][2]][allTpl[i][1]].isConnected = 2;
	                } else {
	                	environment.edges[allTpl[i][1]][allTpl[i][2]].isConnected = 2;
	                	environment.edges[allTpl[i][2]][allTpl[i][1]].isConnected = -2;
	                }
                } else {
					environment.edges[allTpl[i][1]][allTpl[i][2]].isConnected = 6;
	            	environment.edges[allTpl[i][2]][allTpl[i][1]].isConnected = 6;
				}
			}
		}
	}

	//delete
	delete [] oneLineMatrixallTpl;

	//// Save the adjacency matrix
	ss.str("");
	    
	//// Make an adjacency matrix from the list
	ss << environment.outDir << slash << "adjacencyMatrix.miic.orientProba.txt";
	saveAdjMatrix(environment, ss.str());

	// Save the mutualInformation_final...txt with a new name

	spentTime = (std::clock() - startTime_algo) / (double)(CLOCKS_PER_SEC / 1000) /1000;
	if(isVerbose)
		cout << "\n# ----> Update Adjacency Matrix elapsed time:" << spentTime << "sec\n\n";
}
