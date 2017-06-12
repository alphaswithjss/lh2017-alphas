CXX=g++
ifeq ($(shell uname), Darwin)
  CXX = clang++
endif

CXXFLAGS = -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -std=c++14 

HEPMC2INC  =-I$(shell rivet-config --includedir)
PYTHIA8INC = $(shell pythia8-config --cxxflags)
#LHAPDFINC  = $(shell lhapdf-config --cxxflags)

HEPMC2LIBDIR = $(shell rivet-config --libdir)
HEPMC2LIB  = -L$(HEPMC2LIBDIR) -Wl,-rpath,$(HEPMC2LIBDIR) -lHepMC
PYTHIA8LIB = $(shell pythia8-config --libs)
#LHAPDFLIB  = $(shell lhapdf-config --libs)

LIBRARIES += $(PYTHIA8LIB) $(HEPMC2LIB) -lz -ldl 
INCLUDE   += $(PYTHIA8INC) $(HEPMC2INC)

all: 
	make Zj-pythia

Zj-pythia: Zj-pythia.o CmdLine.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBRARIES)

Zj-pythia.o: Zj-pythia.cc
	$(CXX) $(INCLUDE) -o $@ $< -c $(CXXFLAGS)

CmdLine.o: CmdLine.cc
	$(CXX) -o $@ $< -c $(CXXFLAGS)

clean :
	rm -fv Zj-pythia *.o
