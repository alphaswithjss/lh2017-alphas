g++    -c -o qcd-constants.o qcd-constants.cc
g++    -c -o building-blocks.o building-blocks.cc
g++    -c -o resum-SD.o resum-SD.cc
g++    -c -o resum-mMDT.o resum-mMDT.cc
g++ -lgsl -lgslcblas -lm -Wall `root-config --cflags --libs` qcd-constants.o building-blocks.o resum-SD.o resum-mMDT.o NLL.C -o NLL
