# --------------------------------------------- #
# Makefile for EventGen code                    #   
# Pascal Nef, March 6th 2014                    #
# Appended by Ben Nachman an F. Rubbo           #
# Note: source setup.sh before make             #
# --------------------------------------------- #

CXXFLAGS =   -O2 -Wall 

.PHONY: clean debug all

all: dijets

dijets:  dijets.so PUPPI.so CmdLine.so
	$(CXX) dijets.so PUPPI.so CmdLine.so -o $@.exe \
	$(CXXFLAGS) -Wno-shadow  \
	`root-config --glibs`  \
	-L$(FASTJETLOCATION)/lib `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins ` -lVariableR -lSoftKiller -lRecursiveTools \
	-L$(PYTHIA8LOCATION)/lib -lpythia8 -llhapdfdummy \
	-L$(BOOSTLIBLOCATION) -lboost_program_options 

CmdLine.so: CmdLine.cc
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` \
	-I$(PYTHIA8LOCATION)/include \
	-I $(BOOSTINCDIR) \
	`root-config --cflags` 

PUPPI.so: PUPPI.cc PUPPI.hh
	 $(CXX) -o $@ -c $<  \
	 $(CXXFLAGS) -Wno-shadow -fPIC -shared \
        `$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins`  \
        -I$(PYTHIA8LOCATION)/include \
        `root-config --cflags`   

dijets.so : dijets.cc PUPPI.so CmdLine.so
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins`  \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags`	

clean:
	rm -rf *.exe
	rm -rf *.o
	rm -rf *.so
	rm -f *~
