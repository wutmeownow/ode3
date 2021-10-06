ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)
ROOTFLAGS   = $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS) 
CXXFLAGS  += $(ROOTCFLAGS) -I$(ODELIB) -Wall -O3
LDFLAGS    = $(ROOTLIBS) $(ROOTGLIBS) -Wl,-rpath,$(ODELIB) -L$(ODELIB) -lode
GXX	   = g++ $(CXXFLAGS)


EXES=baseball1 baseball2 RKnDemo RKnStep RKnTest


all: $(EXES)

RKnDemo: RKnDemo.cpp $(ODELIB)/RKn.hpp 
	$(GXX) $(CXXFLAGS) -o$@ $@.cpp $(LDFLAGS)

RKnStep: RKnStep.cpp $(ODELIB)/RKn.hpp 
	$(GXX) $(CXXFLAGS) -o$@ $@.cpp $(LDFLAGS)

RKnTest: RKnTest.cpp $(ODELIB)/RKn.hpp 
	$(GXX) $(CXXFLAGS) -o$@ $@.cpp $(LDFLAGS)

baseball1: baseball1.cpp $(ODELIB)/RKn.hpp 
	$(GXX) $(CXXFLAGS) -o$@ $@.cpp $(LDFLAGS) 

baseball2: baseball2.cpp $(ODELIB)/RKn.hpp 
	$(GXX) $(CXXFLAGS) -o$@ $@.cpp $(LDFLAGS)




clean:
	rm -f *.o *.so *.d *~ $(EXES)


