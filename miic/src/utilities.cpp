#include <fstream>
#include <sstream>
#include <string>
#include <getopt.h>
#include <cstdlib>
#include <iostream>
#include <map>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <limits>

#include "utilities.h"

using namespace std;


class sort_indices
{
   private:
     int* mparr;
   public:
     sort_indices(int* parr) : mparr(parr) {}
     bool operator()(int i, int j) const { return mparr[i]<mparr[j]; }
};


void sort2arrays(int len, int a[], int brr[], int bridge[]){
    int i;

    int *pArray = &a[1];
    int *pArray2 = &brr[1];

    std::sort(pArray2, pArray2+len, sort_indices(pArray));

    for(i=1; i < len+1; i++){
    	bridge[i] = pArray2[i-1];
    }

    brr = bridge;
    brr[0] = 0;
}


bool createMemorySpace(Environment& environment, MemorySpace& m){

	int maxLevel = 0;
	for(int i =0; i<environment.numNodes; i++)
		if(environment.allLevels[i] > maxLevel)
			maxLevel = environment.allLevels[i];
	m.maxlevel = maxLevel;

	int nrow=environment.numSamples+1;
	int sampleSize = environment.numSamples;
	int ncol=7;
	int bin_max=maxLevel;
	int iii;

	m.sample = (int **)calloc(nrow, sizeof(int*));
	for(iii = 0; iii < nrow; iii++) 
		m.sample[iii] = (int *)calloc(ncol, sizeof(int));

	m.sortedSample = (int **)calloc(nrow, sizeof(int*));
	for(iii = 0; iii < nrow; iii++) 
	 	m.sortedSample[iii] = (int *)calloc(ncol, sizeof(int));

	m.Opt_sortedSample = (int **)calloc(nrow, sizeof(int*));
	for(iii = 0; iii < nrow; iii++) 
	 	m.Opt_sortedSample[iii] = (int *)calloc(ncol, sizeof(int));

	m.Nxuiz = (int **)calloc(bin_max+1, sizeof(int*));
	for(iii = 0; iii < bin_max+1; iii++) 
		m.Nxuiz[iii] = (int *)calloc(bin_max+1, sizeof(int));


	m.orderSample = (int *)calloc((sampleSize+2), sizeof(int));
	m.sampleKey = (int *)calloc((sampleSize+2), sizeof(int));

	m.Nxyuiz = (int *)calloc((bin_max+1), sizeof(int));
	m.Nyuiz = (int *)calloc((bin_max+1), sizeof(int));
	m.Nuiz = (int *)calloc((bin_max+1), sizeof(int));
	m.Nz = (int *)calloc((bin_max+1), sizeof(int));

	m.Ny = (int *)calloc((bin_max+1), sizeof(int));
	m.Nxui = (int *)calloc((bin_max+1), sizeof(int));
	m.Nx = (int *)calloc((bin_max+1), sizeof(int));		
	m.bridge = (int *)calloc(sampleSize+2, sizeof(int));
}


bool createMemorySpaceThreads(Environment& environment, ContainerMemory& m){

	int maxLevel = 0;
	for(int i =0; i<environment.numNodes; i++)
		if(environment.allLevels[i] > maxLevel)
			maxLevel = environment.allLevels[i];

	int nrow=environment.numSamples+1;
	int sampleSize = environment.numSamples;
	int ncol=7;
	int bin_max=maxLevel;
	int iii;

	m.sortedSample = (int **)calloc(nrow, sizeof(int*));
	for(iii = 0; iii < nrow; iii++) 
	 	m.sortedSample[iii] = (int *)calloc(ncol, sizeof(int));


	m.Nxuiz = (int **)calloc(bin_max+1, sizeof(int*));
	for(iii = 0; iii < bin_max+1; iii++) 
		m.Nxuiz[iii] = (int *)calloc(bin_max+1, sizeof(int));

	m.Nxyuiz = (int *)calloc((bin_max+1), sizeof(int));
	m.Nyuiz = (int *)calloc((bin_max+1), sizeof(int));
	m.Nuiz = (int *)calloc((bin_max+1), sizeof(int));
	m.Nz = (int *)calloc((bin_max+1), sizeof(int));

	m.Ny = (int *)calloc((bin_max+1), sizeof(int));
	m.Nxui = (int *)calloc((bin_max+1), sizeof(int));
	m.Nx = (int *)calloc((bin_max+1), sizeof(int));	
	m.bridge = (int *)calloc(sampleSize+2, sizeof(int));
}


bool deleteMemorySpace(Environment& environment, MemorySpace& m){
	int maxLevel = 0;
	for(int i =0; i<environment.numNodes; i++)
		if(environment.allLevels[i] > maxLevel)
			maxLevel = environment.allLevels[i];
		
	int nrow=environment.numSamples+1;
	int sampleSize = environment.numSamples;
	int ncol=7;
	int bin_max=maxLevel;
	int i;

	for(i=0; i<nrow;i++)
		free(m.sample[i]);
	free(m.sample);

	for(i=0; i<nrow;i++)
		free(m.sortedSample[i]);
	free(m.sortedSample);

	for(i=0; i<bin_max+1;i++)
		free(m.Nxuiz[i]);
	free(m.Nxuiz);


	free(m.orderSample);
	free(m.sampleKey);

	free(m.Nxyuiz);
	free(m.Nyuiz);
	free(m.Nuiz); 
	free(m.Nz);

	free(m.Ny);
	free(m.Nxui);
	free(m.Nx);
	free(m.bridge);
}


bool isOnlyDouble(const char* str) {
    char* endptr = 0;
    strtod(str, &endptr);

    if(*endptr != '\0' || endptr == str)
        return false;
    return true;
}

bool comparatorPairs ( const pair<double,int>& l, const pair<double,int>& r)
   { return l.first < r.first; }

bool SortFunctionNoMore1(const XJAddress* a, const XJAddress* b, const Environment environment) {
	 return environment.edges[a->i][a->j].edgeStructure->Ixy_ui > environment.edges[b->i][b->j].edgeStructure->Ixy_ui;
}

class sorterNoMore {
	  Environment environment;
		public:
	  sorterNoMore(Environment env) : environment(env) {}
	  bool operator()(XJAddress const* o1, XJAddress const* o2) const {
			return SortFunctionNoMore1(o1, o2, environment );
	  }
};

bool SortFunction(const XJAddress* a, const XJAddress* b, const Environment environment) {

	if(environment.edges[a->i][a->j].edgeStructure->status > environment.edges[b->i][b->j].edgeStructure->status)
		return true;
	else if(environment.edges[a->i][a->j].edgeStructure->status < environment.edges[b->i][b->j].edgeStructure->status)
		return false;

	if(environment.edges[a->i][a->j].edgeStructure->status == 1 && environment.edges[b->i][b->j].edgeStructure->status == 1){
		if(environment.edges[a->i][a->j].edgeStructure->Rxyz_ui == 0 && environment.edges[b->i][b->j].edgeStructure->Rxyz_ui != 0)
			return true;
		else if(environment.edges[a->i][a->j].edgeStructure->Rxyz_ui != 0 && environment.edges[b->i][b->j].edgeStructure->Rxyz_ui == 0)
			return false;

		if(environment.edges[a->i][a->j].edgeStructure->Rxyz_ui > environment.edges[b->i][b->j].edgeStructure->Rxyz_ui)
			return true;
		else if(environment.edges[a->i][a->j].edgeStructure->Rxyz_ui < environment.edges[b->i][b->j].edgeStructure->Rxyz_ui)
			return false;
	}

	if(environment.edges[a->i][a->j].edgeStructure->status == 3 && environment.edges[b->i][b->j].edgeStructure->status == 3){
		if(environment.edges[a->i][a->j].edgeStructure->Ixy_ui > environment.edges[b->i][b->j].edgeStructure->Ixy_ui)
			return true;
		else if(environment.edges[a->i][a->j].edgeStructure->Ixy_ui < environment.edges[b->i][b->j].edgeStructure->Ixy_ui)
			return false;
	}
	return false;
}

class sorter {
	  Environment environment;
		public:
	  sorter(Environment env) : environment(env) {}
	  bool operator()(XJAddress const* o1, XJAddress const* o2) const {
			return SortFunction(o1, o2, environment );
	  }
};


bool readTime(string name, ExecutionTime& execTime){
	const char * c = name.c_str();
	ifstream input (c);
	string lineData;
	string s;
	int row = 0;
	int col = 0;
	while(getline(input, lineData))
	{
		if(row == 1){
		istringstream f(lineData);
			while (getline(f, s, '\t')) {
				if(col == 0)
					execTime.init = atof(s.c_str());
				else if(col == 1)
					execTime.iter = atof(s.c_str());
				else if(col == 2)
					execTime.initIter = atof(s.c_str());
				else if(col == 3)
					execTime.initIterSave = atof(s.c_str());
				col++;
			}
		}
		row++;
	}
}


/*
 * Save the adjacency matrix
 */
bool saveAdjMatrix(const Environment environment, const string filename){
	if(environment.isVerbose)
		cout << "Saving adjacency matrix\n";
	ofstream output;
	output.open(filename.c_str());
	for (int i=0; i < environment.numNodes; i++){
		output << environment.nodes[i].name;
		if(i + 1 < environment.numNodes)
			output << "\t";
	}
	output << endl;

	for (int i=0; i < environment.numNodes; i++)
	{	
		output << environment.nodes[i].name << "\t";
		for (int j=0; j < environment.numNodes; j++)
		{
			output << environment.edges[i][j].isConnected;
			if(j + 1 < environment.numNodes)
				output << "\t";
		}
		output << endl;
	} 
	output.close();
}

bool saveAdjMatrixState(const Environment environment, const string filename){
	if(environment.isVerbose)
		cout << "Saving adjacency matrix\n";
	ofstream output;
	output.open(filename.c_str());
	output << "\t";
	for (int i=0; i < environment.numNodes; i++){
		output << environment.nodes[i].name;
		if(i + 1 < environment.numNodes)
			output << "\t";
	}
	output << endl;

	for (int i=0; i < environment.numNodes; i++)
	{	
		output << environment.nodes[i].name << "\t";
		for (int j=0; j < environment.numNodes; j++)
		{
			if(j > i){
				if(environment.edges[i][j].edgeStructure->status == 3)
					output << "1";
			   	else
			   		output << "0";
				if(j + 1 < environment.numNodes)
					output << "\t";
			} else if( i > j){
				if(environment.edges[j][i].edgeStructure->status == 3)
					output << "1";
			   	else
			   		output << "0";
				if(j + 1 < environment.numNodes)
					output << "\t";
			} else {
				output << "0";
				if(j + 1 < environment.numNodes)
					output << "\t";
			}
		}
		output << endl;
	} 
}

/*
 * Transform a vector to a string
 */
string vectorToStringNodeName(Environment environment, const vector<int> vec){
	stringstream ss;
	int length = vec.size();
	if(length > 0){
	  	for (int temp = 0; temp < length; temp++){
	  		if(vec[temp] != -1)
				ss << environment.nodes[vec[temp]].name;
			if(temp+1 < length)
				ss << ",";
		}
	} else {
		ss << "NA";
	}
  	return ss.str();
}

string vectorToString(const vector<int> vec){
	stringstream ss;
	int length = vec.size();
	if(length > 0){
	  	for (int temp = 0; temp < length; temp++){
			ss << vec[temp];
			if(temp+1 < length)
				ss << ",";
		}
	}
  	return ss.str();
}

string arrayToString1(const double* int_array, const int length){
	stringstream ss;
	if(length > 0){
	  	for (int temp = 0; temp < length; temp++){
	  		if(int_array[temp] != -1)
				ss << int_array[temp] << ", ";
		}
	} else {
		ss << "NA";
	}
  	return ss.str();
}

string zNameToString(Environment environment, vector<int> vec, int pos){
	stringstream ss;
	if(pos != -1)
		ss << environment.nodes[vec[pos]].name;
	else
		ss << "NA";
	return ss.str();
}

/*
 * Save the edges list with their properties (Info, ui, zi..)
 */
bool saveEdgesListAsTable(const Environment environment, const string filename){
	
	bool sortedPrint = true;

	if(environment.isVerbose)
		cout << "Saving edges file\n";
	ofstream output;
	output.open(filename.c_str());

	if(sortedPrint){
		vector<XJAddress*> allEdges;

		for(int i = 0; i < environment.numNodes -1; i++){
		 	for(int j = i + 1; j < environment.numNodes; j++){
		 		XJAddress* s = new XJAddress();
				s->i=i;
				s->j=j;
				allEdges.push_back(s);
		 	}
		}

		std::sort(allEdges.begin(), allEdges.end(), sorter(environment));


		output << "x" << "\t" << "y" << "\t" << "z.name" << "\t" << "ui.vect" << "\t" << "zi.vect" << "\t" 
				<< "Ixy_ui" << "\t" << "cplx" << "\t" << "Rxyz_ui" << "\t" << "category" << "\t" << "Nxy_ui\n";

		for(int i = 0; i < allEdges.size();i++){

			for(int j = 0; j < environment.numNodes;j++){
				if(j == allEdges[i]->i || j == allEdges[i]->j)
					output << "1";
				else
					output << "0";
			}
			output << "\t";

			output << environment.nodes[allEdges[i]->i].name << "\t" << environment.nodes[allEdges[i]->j].name << "\t" << 
		 		zNameToString(environment, environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->zi_vect_idx, environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->z_name_idx) << "\t"
	 			 << vectorToStringNodeName(environment, environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->ui_vect_idx) << "\t"
	 			 << vectorToStringNodeName(environment, environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->zi_vect_idx) << "\t"
	 			 << environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->Ixy_ui << "\t" << environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->cplx << "\t"
	 			 << environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->Rxyz_ui << "\t" << environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->status << "\t"
	 			 << environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->Nxy_ui;
		 	output << "\n";	

		}
	} else{
		output << "x" << "\t" << "y" << "\t" << "z.name" << "\t" << "ui.vect" << "\t" << "zi.vect" << "\t" 
				<< "Ixy_ui" << "\t" << "cplx" << "\t" << "Rxyz_ui" << "\t" << "category" << "\t" << "Nxy_ui\n";

		for(int i = 0; i < environment.numNodes -1; i++){
			for(int j = i + 1; j < environment.numNodes; j++){
				output << environment.nodes[i].name << "\t" << environment.nodes[j].name << "\t" << 
				zNameToString(environment, environment.edges[i][j].edgeStructure->zi_vect_idx, environment.edges[i][j].edgeStructure->z_name_idx) << "\t"
					 << vectorToStringNodeName(environment, environment.edges[i][j].edgeStructure->ui_vect_idx) << "\t"
					 << vectorToStringNodeName(environment, environment.edges[i][j].edgeStructure->zi_vect_idx) << "\t"
					 << environment.edges[i][j].edgeStructure->Ixy_ui << "\t" << environment.edges[i][j].edgeStructure->cplx << "\t"
					 << environment.edges[i][j].edgeStructure->Rxyz_ui << "\t" << environment.edges[i][j].edgeStructure->status << "\t"
					 << environment.edges[i][j].edgeStructure->Nxy_ui;
 				output << "\n";

			}
		}
	}
	output.close();
}

/*
 * Save the runtime file
 */
bool saveExecTime(const Environment environment, ExecutionTime execTime, const string filename){
	if(environment.isVerbose)
		cout << "Saving execution time\n";
	ofstream output;
	output.open(filename.c_str());

	output << "init" <<	"\t" << "iter" << "\t" << "initIter" << "\t" << "initIterSave" << "\n";
	output << execTime.init <<	"\t" << execTime.iter << "\t" << execTime.initIter << "\t" << execTime.initIterSave;
}

bool existsTest(const string& name) {
	  ifstream f(name.c_str());
	if (f.good()) {
		f.close();
		return true;
	} else {
		f.close();
		return false;
	}   
}

int** copyMatrix(int**oldmatrix, int numRows, int numColumns){
	int** newMatrix;
	newMatrix = new int*[numRows];
	for(int i = 0; i < numRows; i++)
		newMatrix[i] = new int[numColumns];

	for(int i = 0; i < numRows; i++){
		for(int j = 0; j < numColumns; j++){
			newMatrix[i][j] = oldmatrix[i][j];
		}
	}
	return newMatrix;
}

bool checkNA(int** data, int numRows, int numColumns){
	for(int i = 0; i < numRows; i++)
		for(int j = 0; j < numColumns; j++)
			if(data[i][j] == -1)
				return true;

	return false;
}

string printNodesName(Environment environment){
	string s = "";
	for(int i = 0; i < environment.numNodes; i++){
		cout << environment.nodes[i].name;
		if(i + 1 < environment.numNodes)
			cout << " ";
	}
	cout << "\n";
	return s;
}

/*
 * Print input data, as strings or factors
 */
void printMatrix(Environment environment, string type){
	if(type.compare("string") == 0){
		cout << "Data matrix of strings\n";
		printNodesName(environment);
		for(int i = 0; i < environment.numSamples; i++){
			for(int j = 0; j < environment.numNodes; j++){
				cout << environment.data[i][j] << " ";
			}
			cout << endl;
		}
	} else if(type.compare("factors") == 0){
		cout << "Data matrix of factors\n";
		printNodesName(environment);
		for(int i = 0; i < environment.numSamples; i++){
			for(int j = 0; j < environment.numNodes; j++){
				cout << environment.dataNumeric[i][j] << " ";
			}
			cout << endl;
		}
	}
}

/*
 * Initialize all the elements of the array to the given value
 */
bool setArrayValuesInt(int* array, int length, int value){
	for(int i = 0; i < length; i++){
			array[i] = value;
		}
		return true;
}

/*
 * Check if a value is of type integer
 */
bool isInteger(const string &s)
{
   if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false ;

   char * p ;
   strtol(s.c_str(), &p, 10) ;

   return (*p == 0) ;
}

/* 
 * Read data matrix from file
 */
 bool readData(Environment& environment, bool& isNA){

 	struct stat buffer; 
	if((!stat (environment.inData.c_str(), &buffer) == 0)){
		cout << "Input file does not exist. STOP!";
		exit(1);
	}
		

 	const char * c = environment.inData.c_str();
	ifstream input (c);
	string lineData;

	environment.numSamples = 0;
	environment.numNodes = 0;

	while(getline(input, lineData))
	{
		environment.numSamples++;
		string d;
		stringstream lineStream(lineData);
		if(environment.numNodes == 0){
			while(std::getline(lineStream, d, '\t')) {
				environment.numNodes++;
			}
		}
	}

	ifstream input1 (c);
	//cout << environment.numSamples << " " << environment.numNodes << endl;
	
	//initialize matrix to void string, skip the header
	environment.data = new string*[environment.numSamples - 1];
	for(int i = 0; i < environment.numSamples -1; i++){
		environment.data[i] = new string[environment.numNodes];
		for(int j = 0; j < environment.numNodes; j++){
			environment.data[i][j] = "";
		}
	}

	environment.numSamples--;
	
	//load the input as string
	int i = -1;
	while(getline(input1, lineData))
	{
		int j = 0;
		istringstream lineStream(lineData);
		string d;

		// save the name of the nodes
		if(i == -1){
			environment.nodes = new Node[environment.numNodes];
			while(std::getline(lineStream, d, '\t')) {

				environment.nodes[j].name = d;
				j++;
			}
		} else {

			while(std::getline(lineStream, d, '\t')) {
				if((d.compare("NA")  == 0) || (d.compare("") == 0)){
					isNA = true;
				}

				environment.data[i][j] = d;

				j++;
			}
		}
		i++;
	}

	// set effN if not set before to the number of rows from the input environment.data
	if(environment.effN == -1)
		environment.effN = environment.numSamples;

	//set maxbin coarse
	environment.maxbins=50;

		
	return true;
}

/* 
 * Remove all the lines that contain only NA
 */ 
bool removeRowsAllNA(Environment& environment){
	if(environment.isVerbose)
		cout << "# Removing NA rows\n";

	int* indexNA = new int[environment.numSamples];
	setArrayValuesInt(indexNA, environment.numSamples, -1);

	int pos = 0;
	for(int i = 0; i < environment.numSamples; i++){
		bool isNA = true;
		for(int j = 0; j < environment.numNodes && isNA; j++){
			if((environment.data[i][j].compare("NA")  != 0) && (environment.data[i][j].compare("") != 0)){
				isNA = false;
			}
		}
		if(!isNA){
			indexNA[pos] = i;
			pos++;
		}
	}


	if(environment.isVerbose)
		cout << "-------->Number of NA rows: " << pos << "\n";
	cout << "environment.numSamples :" << environment.numSamples << endl;

	// if there are rows of NA value
	if(pos != 0){		
		//correct variable numSamples

		// save the values
		int pos = 0;
		for(int i = 0; i < environment.numSamples; i++){
			if(indexNA[i] != -1){
				for(int j = 0; j < environment.numNodes; j++){
					environment.data[pos][j] = environment.data[indexNA[i]][j];
				}
				pos++;
			}
		}
		environment.numSamples = pos;
		if(environment.effN > pos){
			environment.effN = pos;
		}
	}
	cout << "environment.numSamples :" << environment.numSamples << endl;
	return true;
}

/* 
 * Transforms the string into factors
 */ 
bool transformToFactors(Environment& environment, int i){
	if(environment.isVerbose)
		cout << "# Transforming matrix to factors\n";

	 // create a dictionary to store the factors of the strings
 	map<string,int> myMap;

	//clean the dictionary since it is used column by column
	myMap.clear();
	myMap["NA"] = -1;
	myMap[""] = -1;
	int factor = 0;

 	for(int j = 0; j < environment.numSamples; j++){
		map<string,int>::iterator it = myMap.find(environment.data[j][i]);
		if ( it != myMap.end() ){
			environment.dataNumeric[j][i] = it->second;
		}
		else {
			myMap[environment.data[j][i]] = factor;
			environment.dataNumeric[j][i] = factor;
			factor++;

		}
	}

}

bool copyValue(Environment& environment, int i){

 	for(int j = 0; j < environment.numSamples; j++){
		environment.dataNumeric[j][i] = atof(environment.data[j][i].c_str());
	}
}

/* 
 * Set the number of levels for each node (the maximum level of each column)
 */ 
bool setNumberLevels(Environment& environment){
	int max;
 	environment.allLevels = new int[environment.numNodes];

	for(int i = 0; i < environment.numNodes;i++){
		max = 0;
	 	for(int j = 0; j < environment.numSamples; j++){
	 		if(environment.dataNumeric[j][i] > max)
	 			max = environment.dataNumeric[j][i];
	 	}


	 	environment.allLevels[i] = max+1;
	 	// cout << environment.nodes[i].name << " has "  << environment.allLevels[i] << " levels\n";
	 }
}

/*
 * Parse the command line inserted by the user
 */
bool parseCommandLine(Environment& environment, int argc, char** argv) {
	int c;
	string s;
	// set the default values if not the test case
	if(!environment.myTest){
		environment.inData = "";
		environment.outDir = "";
		environment.blackbox_name = "";
		environment.edgeFile = "";
		environment.sampleWeightsFile="";
		environment.effN = -1;
		environment.cplx = 1;
		environment.isVerbose = false;
		environment.numberShuffles = 0;
		environment.isLatent = false;
		environment.isTplReuse = true;
		environment.isK23 = true;
		environment.isDegeneracy = false;
		environment.isNoInitEta = false;
		environment.isPropagation = true;
		environment.halfVStructures = 0;
		environment.nThreads = 0;
		cout << 	environment.sampleWeightsFile << endl;

	}
	
	// parse the command line

	while ((c = getopt (argc, argv, "i:o:b:d:c:e:s:r:q:k:n:p:a:h:m:t:u:z:x:glfv?")) != -1){
		switch (c){
			case 'i':{
				environment.inData.append(optarg);
				break;
			}
			case 'o':{
				environment.outDir.append(optarg);		
				break;
			}
			case 'x':{
				environment.seed = 0;
				s = optarg;
				if(!isInteger(s)){
					cout << "[ERR] Seed should be an integer!\n";
				}
				else
					environment.seed = atoi(optarg);
				break;
			}
			case 'b':{
				environment.blackbox_name.append(optarg);		
				break;
			}
			case 'd': {
				s = optarg;
				stringstream ss(s); // Turn the string into a stream.
				string tok;
				char delimiter = ',';
				while(getline(ss, tok, delimiter)) {
					int ival = atoi(tok.c_str());
					environment.steps.push_back(ival);
				}
				break;
			}
			case 'n':{
				s = optarg;
				if(!isInteger(s) && atoi(optarg) > 1) {
					cout << "[ERR] EffeN should be a positive integer!\n";
				}
				else
					environment.effN = atoi(optarg);
				break;
			}
			case 'u':{
				s = optarg;
				environment.dataTypeFile.append(optarg);
				break;
			}
			case 'q':{
				environment.sampleWeightsFile.append(optarg);	
				break;
			}
			case 'm':{
				environment.edgeFile.append(optarg);
				break;		
			}
			case 'z':{
				environment.nThreads = atoi(optarg);
				break;		
			}
			case 'c':{
				environment.cplxType.append(optarg);

				if(environment.cplxType.compare("mdl")!=0 & environment.cplxType.compare("nml")!=0){
					cout << "[ERR] Wrong complexity check option!\n";
					exit(1);
				} else if(environment.cplxType.compare("mdl") == 0){
					environment.cplx = 0;
				}
				break;
			}
			case 'e':{
				s = optarg;
				if(!isOnlyDouble(optarg)){
					cout << "[ERR] Confidence cut should be a double!\n";
					exit(1);
				} 

				environment.confidenceThreshold = atof(optarg);
				break;
			}
			case 's':{
				s = optarg;
				if(!isInteger(s)){
					cout << "[ERR] Shuffle should be an integer!\n";
					exit(1);
				} 
				else
					environment.numberShuffles = atoi(optarg);
				break;
			}
			case 'r':{
				s = optarg;
				if(s.compare("1") != 0 && s.compare("0") != 0){
					cout << "[ERR] Wrong reuse/not reuse tpl argument!\n";
					return false;
				} else if(s.compare("0") == 0)
					environment.isTplReuse = false;
				break;
			}
			case 'k': {
				s = optarg;
				if(s.compare("1") != 0 && s.compare("0") != 0){
					cout << "[ERR] Case k: Wrong k23 argument!\n";
					exit(1);
				} else if(s.compare("0") == 0)
					environment.isK23 = false;
				break;
			}
			case 'p' : {
				s = optarg;
				if(s.compare("0") != 0){
					cout << "[ERR] Case p: Wrong propagation argument!*" << s <<"*\n";
					exit(1);
				} else if(s.compare("0") == 0)
					environment.isPropagation = false;
				break;
			}
			case 'a' : {
				s = optarg;
				if(s.compare("0") != 0 && s.compare("1") != 0){
					cout << "[ERR] Case a: Wrong half V-structures argument!*" << s <<"*\n";
					exit(1);
				} else
					environment.halfVStructures = atoi(optarg);
				break;
			}
			case 'l':{
				environment.isLatent = true;
				break;
			}
			case 'g': {
				environment.isDegeneracy = true;
				break;
			}
			case 'f': {
				environment.isNoInitEta = true;
				break;
			}
			case 'v':{
				environment.isVerbose = true;
				break;
			}


			case '?':{
				if (optopt == 'c')
					fprintf (stderr, "Option -%c requires an argument.\n", optopt);
				else if (isprint (optopt))
					fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				else
				fprintf (stderr,"Unknown option character `\\x%x'.\n", optopt);
				exit(1);
			}
		}
	}

	if(!environment.inData.compare("")){
		cout << "The input data file is required (-i)\n";
		exit(1);
	}

	if(!existsTest(environment.inData)){
		cout << "The input file does not exist\n";
		exit(1);
	}

	if(!environment.outDir.compare("")){
		cout << "The output dir path is required (-o)\n";
		exit(1);
	}

	if(environment.steps.empty()){
		environment.steps.push_back(1);
		environment.steps.push_back(2);
	}



	return true;
}

/* 
 * Print the most important variables of the environment
 */ 
void printEnvironment(Environment environment, Log* pLog){
	//// Recall the main parameters
	stringstream s;
	s << "# --------\n# Inputs:\n# ----\n"
		<< "# Input data file --> " << environment.inData << "\n"
		<< "# Output directory --> " << environment.outDir << "\n"
		<< "# All properties --> " << printNodesName(environment) << "\n"
		<< "# Eff. N --> " << environment.effN << "\n"
		<< "# Thres. Pc --> " << environment.thresPc << "\n"
		<< "# N min --> " << environment.minN << "\n"
		<< "# Clpx type --> " << environment.cplxType << "\n"
		<< "# Clpx check --> " << environment.cplx << "\n"
		<< "# Latent --> " << environment.isLatent << "\n"
		<< "# Reuse --> " << environment.isTplReuse << "\n"
		<< "# K23 --> " << environment.isK23 << "\n"
		<< "# propagation --> " << environment.isPropagation << "\n"
		<< "# half V structures --> " << environment.halfVStructures << "\n"
		<< "# Degeneracy --> " << environment.isDegeneracy << "\n"
		<< "# No Init Eta --> " << environment.isNoInitEta << "\n"
		<< "# VERSION --> " << environment.myVersion << "\n"
		<< "# --------\n";

		cout << s.str();

		char* msg =  const_cast<char*>  (s.str().c_str());
		pLog->write(msg);
}


/* 
 * Set the variables in the environment structure
 */ 
bool setEnvironment(Environment& environment){
	// Load the data
	// ----
	environment.noMoreAddress.clear();
	environment.numNoMore = 0;
	environment.searchMoreAddress.clear();
	environment.numSearchMore = 0;


	if( environment.isVerbose == true )
	{ cout << "# --------\n# Load the input data...\n"; }
	
	bool isNA = false;

	readData(environment, isNA);

	if(isNA){
		//// Remove the lines that are all 'NA'
		removeRowsAllNA(environment);
	}
	

	//create the data matrix for factors
	environment.dataNumeric = new int*[environment.numSamples];
	for(int i = 0; i < environment.numSamples; i++){
		environment.dataNumeric[i] = new int[environment.numNodes];
	}


	//transform to factors
	for(int i = 0; i < environment.numNodes; i++){
		transformToFactors(environment, i);
	}

	if(environment.effN > environment.numSamples)
		environment.effN = environment.numSamples;


	//// Set the effN if not already done
	if(environment.effN == -1 )
		environment.effN = environment.numSamples;

	//// Set a variables with all properties name and levels
	setNumberLevels(environment);

	// create the 1000 entries to store c2 values
	environment.c2terms = new double[environment.numSamples+1];
	for(int i = 0; i < environment.numSamples+1; i++){
		environment.c2terms[i] = -1;
	}

	//// Set the number of digits for the precision while using round( ..., digits = ... )
	//// Make sure the min levels for the data is 0
	environment.minN = 1;

	//// Set the probability threshold for the rank
	environment.thresPc = 0;	// if the contribution probability is the min value
	environment.l = (environment.numNodes*(environment.numNodes-1)/2);


	// create the edge structure and keep track of how many searchMore we have
	environment.edges = new Edge*[environment.numNodes];

	for(int i = 0; i < environment.numNodes; i++)
		environment.edges[i] = new Edge[environment.numNodes];

	for(int i = 0; i < environment.numNodes; i++){
		for(int j = 0; j < environment.numNodes; j++){
			environment.edges[i][j].isConnected = 1;
		}
	}
}

bool readFilesAndFillStructures(Environment& environment, string slash){
	//fill nodes
	setEnvironment(environment);

	//create the one line matrix
	environment.oneLineMatrix = new int[environment.numSamples*environment.numNodes];
	for(int i = 0; i < environment.numSamples;i++){
		for(int j = 0; j < environment.numNodes;j++){
			// cout << j * environment.numSamples + i << " ";
			environment.oneLineMatrix[j * environment.numSamples + i] = environment.dataNumeric[i][j];
		}
	}

	//create edges
	environment.edges = new Edge*[environment.numNodes];

	for(int i = 0; i < environment.numNodes; i++)
		environment.edges[i] = new Edge[environment.numNodes];

	for(int i = 0; i < environment.numNodes; i++)
		environment.edges[i][i].isConnected = 0;

	for(int i = 0; i < environment.numNodes - 1; i++){
		for(int j = i + 1; j < environment.numNodes; j++){
			// create a structure for the nodes that need to store information about them
			environment.edges[i][j].edgeStructure = new EdgeStructure();
			environment.edges[j][i].edgeStructure = environment.edges[i][j].edgeStructure ;
			// initialize the structure
			environment.edges[j][i].edgeStructure->z_name_idx = -1;
			environment.edges[j][i].edgeStructure->status = -1;
		}
	}
	std::stringstream ss;
	ss.str("");
	if(environment.edgeFile.compare("") == 0)
		ss << environment.outDir << slash << "edgesList.miic.txt";
	else
		ss << environment.edgeFile;
	string filename = ss.str();

	//fill edges reading the edge table
	struct stat buffer; 
	if((!stat (filename.c_str(), &buffer) == 0)){
		cout << "Input file does not exist. STOP!";
		exit(1);
	}
		

 	const char * c = filename.c_str();
	ifstream input (c);
	string lineData;
	string s;
	int row = 0;
	int col = 0;
	int posX = -1;
	int posY = -1;

	while(getline(input, lineData))
	{
		if(row != 0){
			col = 0;
			istringstream f(lineData);
			while (getline(f, s, '\t')) {
				if(col == 0){			}
				if(col == 1){
					for(int i = 0; i < environment.numNodes;i++)
						if(environment.nodes[i].name.compare(s) == 0)
							posX=i;

				}
				else if(col == 2){
					for(int i = 0; i < environment.numNodes;i++)
						if(environment.nodes[i].name.compare(s) == 0)
							posY=i;
				}
				else if(col == 3){
				}
				else if(col == 4){
					if(s.compare("NA") != 0){
						stringstream ss(s); // Turn the string into a stream.
						string tok;
						char delimiter = ',';
						while(getline(ss, tok, delimiter)) {
							int ival;
							for(int i = 0; i < environment.numNodes;i++){
								if(environment.nodes[i].name.compare(tok) == 0){
									ival=i;
									break;
								}
							}
							environment.edges[posX][posY].edgeStructure->ui_vect_idx.push_back(ival);
						}
					}
				} else if(col == 5){
					if(s.compare("NA") != 0){
						stringstream ss(s); // Turn the string into a stream.
						string tok;
						char delimiter = ',';
						while(getline(ss, tok, delimiter)) {
							int ival;
							for(int i = 0; i < environment.numNodes;i++){
								if(environment.nodes[i].name.compare(tok) == 0){
									ival=i;
									break;
								}
							}
							environment.edges[posX][posY].edgeStructure->zi_vect_idx.push_back(ival);
						}
					}
				} else if(col == 6){
					environment.edges[posX][posY].edgeStructure->Ixy_ui = atof(s.c_str());
				} else if(col == 7){
					environment.edges[posX][posY].edgeStructure->cplx = atof(s.c_str());
				} else if(col == 8){
					environment.edges[posX][posY].edgeStructure->Rxyz_ui = atof(s.c_str());
				} else if(col == 9){
					int state = atoi(s.c_str());
					environment.edges[posX][posY].edgeStructure->status = state;
					if(state == 3){
						environment.edges[posX][posY].isConnected = 1;
						environment.edges[posY][posX].isConnected = 1;
						// add the edge to Nomore
						XJAddress* ij = new XJAddress();
						ij->i = posX;
						ij->j = posY;  
						environment.noMoreAddress.push_back(ij);	
					}
					else{
						environment.edges[posX][posY].isConnected = 0;
						environment.edges[posY][posX].isConnected = 0;
					}
				} else if(col == 10){
					environment.edges[posX][posY].edgeStructure->Nxy_ui = atof(s.c_str());
				}
				col++;
			}
		}
		row++;
	}

	environment.numNoMore = environment.noMoreAddress.size();

	std::sort(environment.noMoreAddress.begin(), environment.noMoreAddress.end(), sorterNoMore(environment));
}

// read the list of the edges to remove
bool readBlackbox(Environment& environment, string slash){
	std::stringstream ss;
	ss.str("");
	ss << environment.blackbox_name;
	string filename = ss.str();

	const char * c = filename.c_str();
	ifstream input (c);
	string lineData;
	string s;
	string s1;
	string s2;
	int col;
	int posX;
	int posY;

	while(getline(input, lineData)){
		col = 0;
		posX = -1;
		posY = -1;
		istringstream f(lineData);
		while (getline(f, s, '\t')) {
			if(col == 0){
				s1 = s;
				for(int i = 0; i < environment.numNodes;i++){
					if(environment.nodes[i].name.compare(s) == 0)
						posX=i;
				}

			}
			else if(col == 1){
				s2 = s;
				for(int i = 0; i < environment.numNodes;i++){
					if(environment.nodes[i].name.compare(s) == 0)
						posY=i;
				}
			}
			col++;
		}
		if(posX != -1 && posY != -1){
			environment.edges[posX][posY].isConnected = 0;
			environment.edges[posY][posX].isConnected = 0;
		}
		else
			cout << "[WARNING] Edge " << s1 << "-" << s2 << " is not an edge of the network\n";
	}

	return true;
}
