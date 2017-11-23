#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include "structure.h"
#include "log.h"
bool createMemorySpace(Environment&, MemorySpace&);
bool createMemorySpaceThreads(Environment&, ContainerMemory&);
bool deleteMemorySpace(Environment&, MemorySpace&);
std::string printNodesName(Environment);
void printMatrix(Environment, std::string);
void readData(Environment&);
bool parseCommandLine(Environment&, int, char**);
void printEnvironment(Environment, Log*);
bool setEnvironment(Environment&);
int** copyMatrix(int**, int, int);
bool setNumberLevels(Environment&);
bool existsTest(const std::string&);
bool checkNA(int**, int, int);
bool saveExecTime(const Environment, const std::string);
bool saveAdjMatrix(const Environment, const std::string);
double findAvg(const Environment, double**, int);
bool saveEdgesListAsTable(const Environment, const std::string);
std::string arrayToString(Environment, const int* , const int );
std::string arrayToString1(const double*, const int);
std::string vectorToString(const std::vector<int> vec);
bool readTime(Environment&, std::string);
bool readFilesAndFillStructures(Environment&, std::string);
int sign(double val);
bool readBlackbox(Environment&, std::string);
bool removeBlackboxEdges(Environment& environment);
void sort2arrays(int len, int a[], int brr[], int bridge[]);
vector< vector <string> > saveEdgesListAsTable1(Environment& environment);
void deleteMemorySpaceThreads(Environment& environment, ContainerMemory& m);
#endif
