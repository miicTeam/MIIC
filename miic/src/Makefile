CXX = g++ 
CFLAGS = -c -O3 -static -pthread
#CFLAGS = -g -O0 -c -static -pthread

LDFLAGS = -static -pthread -O3
#LDFLAGS = -g -O0 -static -pthread

CXXSOURCES = main.cpp skeleton.cpp utilities.cpp log.cpp computeEnsInformation.cpp skeletonInitialization.cpp skeletonIteration.cpp orientationProbability.cpp probaOrientation.cpp probaOrientation_interface.cpp confidenceCut.cpp computeInfo.cpp computeInfo_interface.cpp
OBJECTS = $(CXXSOURCES:.cpp=.o) $(CSOURCES:.c=.o)
EXECUTABLE = miic

all: $(SOURCES) $(EXECUTABLE)
	cp miic ../miic
	
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $(EXECUTABLE)

.cpp.o:
	$(CXX) $(CFLAGS) $< -o $@

clean:
	rm -rf $(EXECUTABLE) *.o *~

.PHONY: all clean
