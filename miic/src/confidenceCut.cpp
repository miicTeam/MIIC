#include <math.h>
#include <algorithm>
#include <iostream>
#include <sstream>

#include "structure.h"
#include "computeEnsInformation.h"
#include "utilities.h"
using namespace std;

void shuffle_lookup(int *array, size_t n)
{
    if (n > 1) 
    {
        size_t i;
        for (i = 0; i < n - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}


bool saveConfidenceVector(Environment& environment, int** inferredEdges_tab, double* confVect,  string outDir, string slash){
	stringstream ss;
	ss.str("");
	ss << outDir << slash << "confRatios.txt";
	string filename = ss.str();

	if(environment.isVerbose)
		cout << "Saving confidence matrix\n";
	ofstream output;
	output.open(filename.c_str());
	output << "x" << "\t" << "y" << "\t" << "confidence_ratio" << endl;
	for (int i=0; i < environment.numNoMore; i++){
		int X = inferredEdges_tab[i][0] ;
 		int Y = inferredEdges_tab[i][1] ;

		output << environment.nodes[X].name << "\t" << environment.nodes[Y].name << "\t" << confVect[i] << endl;
	}
	
	output.close();
}

bool SortFunctionNoMore2(const XJAddress* a, const XJAddress* b, const Environment& environment) {
	 return environment.edges[a->i][a->j].edgeStructure->Ixy_ui > environment.edges[b->i][b->j].edgeStructure->Ixy_ui;
}

class sorterNoMore2 {
	  Environment environment;
		public:
	  sorterNoMore2(Environment env) : environment(env) {}
	  bool operator()(XJAddress const* o1, XJAddress const* o2) const {
			return SortFunctionNoMore2(o1, o2, environment );
	  }
};

bool confidenceCut(Environment& environment, string slash){

	int** safe_state;

	int* lookup = new int[environment.numSamples];
	for(int i=0; i<environment.numSamples; i++)
		lookup[i] = i;

	double noMore = environment.numNoMore;

	// Allocate the true edges table
	int** inferredEdges_tab;
	inferredEdges_tab = new int*[environment.numNoMore];
	for(int i=0; i < environment.numNoMore; i++)
		inferredEdges_tab[i] = new int[2];

	int pos = 0;

	for(int i = 0; i < environment.numNodes -1; i++){
		for(int j = i+1; j < environment.numNodes; j++){
			if(environment.edges[i][j].isConnected){
				inferredEdges_tab[pos][0] = i;
				inferredEdges_tab[pos][1] = j;
				pos++;
			}
		}
	}


	// Create a back up of the data, for later randomization
	safe_state = new int*[environment.numSamples];
	for(int i=0; i<environment.numSamples; i++)
		safe_state[i] = new int[environment.numNodes];

	//copy to safe state
	for(int i=0; i<environment.numSamples; i++)
	{
		for(int j=0; j<environment.numNodes;j++){
			safe_state[i][j]=environment.dataNumeric[i][j];
		}
	}

	int* nodes_toShf = new int[environment.numNodes];
  	for(int i=0; i<environment.numNodes; i++)
  		nodes_toShf[i]=0;
  	
  	//indexes of nodes to shuffle
  	for(int i=0; i<environment.numNoMore; i++)
	{
		nodes_toShf[inferredEdges_tab[i][0]] = 1;
	}

	double* confVect = new double[environment.numNoMore];
	for(int nb=0;nb<environment.numNoMore;nb++)
		confVect[nb] = 0;


	int* ptrVarIdx = new int[2];

	// loop on the number of shuffling
	for(int nb=1; nb<=environment.numberShuffles; nb++)
	{
		//Shuffle the dataset only for the variables present in nodes_toShf
		for(int col=0; col<environment.numNodes; col++)
		{
			if( nodes_toShf[col] == 1)
			{
				shuffle_lookup(lookup,environment.numSamples);
				for(int row=0; row<environment.numSamples; row++){
					environment.dataNumeric[row][col] = safe_state[lookup[row]][col];
				}
			}
		}

		for(int i = 0; i < environment.numSamples;i++){
			for(int j = 0; j < environment.numNodes;j++){
				environment.oneLineMatrix[j * environment.numSamples + i] = environment.dataNumeric[i][j];
			}
		}

		int X,Y;
		int N_xy_ui;
		double NIxy_ui,k_xy_ui;
		// evaluate the mutual information for every edge
		for(int i=0; i<environment.numNoMore; i++)
		{
		 	X = inferredEdges_tab[i][0] ;
	 		Y = inferredEdges_tab[i][1] ;

			double* res; 

			res = computeEnsInformationNew(environment, NULL, 0, NULL, 0, -1, X, Y, environment.cplx);
			N_xy_ui = res[0]; 
			NIxy_ui = res[1];
			k_xy_ui = res[2];
			free(res);

			double ni = NIxy_ui-k_xy_ui;
			if(ni <= 0){ ni = 0 ;}
			confVect[i] += exp(-ni);
	    }
	}

	//evaluate the average confidence
	for(int nb=0;nb<environment.numNoMore;nb++){
		confVect[nb] /= environment.numberShuffles;	
	}

	//put values > 1 to 1
	for(int nb=0;nb<environment.numNoMore;nb++)
		if(confVect[nb] > 1) 
			confVect[nb] = 1;	

	//remove edges based on confidence cut
	double confidence;
	vector<int> toDelete ;

	for(int i = 0; i < environment.numNoMore; i++){
		int X = inferredEdges_tab[i][0] ;
 		int Y = inferredEdges_tab[i][1] ;
 		confidence = exp (- (environment.edges[X][Y].edgeStructure->Ixy_ui - environment.edges[X][Y].edgeStructure->cplx));
 		confVect[i] = confidence / confVect[i];
 		if(confVect[i] > environment.confidenceThreshold){
 			if(X > Y)
 				environment.edges[X][Y].edgeStructure->status = 1; 
 			else
 				environment.edges[Y][X].edgeStructure->status = 1; 
			environment.edges[X][Y].isConnected = 0;
			environment.edges[Y][Y].isConnected = 0;


			toDelete.push_back(i);
 		}
	}

	cout << "\n# -- number of edges cut: " << toDelete.size() << "\n";

	//delete from vector
	for(int i = 0 ; i < environment.numNoMore; i++)
		delete environment.noMoreAddress[i];
	environment.noMoreAddress.clear();

	
	for(int i = 0; i < environment.numNoMore; i++){
		if(!(std::find(toDelete.begin(), toDelete.end(), i) != toDelete.end())) {
			int X = inferredEdges_tab[i][0] ;
 			int Y = inferredEdges_tab[i][1] ;
			XJAddress* s = new XJAddress();
			s->i=X;
			s->j=Y;
			environment.noMoreAddress.push_back(s);
		}
	}

	//put back values without orientation	
	for(int X = 0; X < environment.numNodes -1; X++){
		for(int Y = X+1; Y < environment.numNodes; Y++){
			if(environment.edges[X][Y].isConnected == -2 || environment.edges[X][Y].isConnected == 2 || environment.edges[X][Y].isConnected == 6){
	 			environment.edges[X][Y].isConnected = 1;
	 			environment.edges[Y][X].isConnected = 1;
	 		}
	 	}
	 }



	//////////////////////////////////////////////// copy data back
	for(int i=0; i<environment.numSamples; i++)
	{
		for(int j=0; j<environment.numNodes;j++){
			environment.dataNumeric[i][j]=safe_state[i][j];
		}
	}

	for(int i = 0; i < environment.numSamples;i++){
		for(int j = 0; j < environment.numNodes;j++){
			environment.oneLineMatrix[j * environment.numSamples + i] = environment.dataNumeric[i][j];
		}
	}

	saveConfidenceVector(environment, inferredEdges_tab, confVect, environment.outDir, slash);


	std::sort(environment.noMoreAddress.begin(), environment.noMoreAddress.end(), sorterNoMore2(environment));
	environment.numNoMore = environment.noMoreAddress.size();

		delete [] ptrVarIdx;

	for(int i=0; i<environment.numSamples; i++)
		delete [] safe_state[i];
	delete [] safe_state;

	delete [] lookup;

	delete [] nodes_toShf;

	for(int i=0; i < noMore; i++)
		delete [] inferredEdges_tab[i];
	delete [] inferredEdges_tab;
	delete [] confVect; 
	
	return true;
}