
CXX = g++

CPPFLAGS = -g -fno-inline -O0

ROOTFLAGS = `root-config --cflags --glibs`
ROOTINC = `root-config --cflags`
ROOTLIB = `root-config --glibs`

INCLUDE += -I$(GENIE)/src -I$(LOG4CPP_INC)

# If 32-bit, you may need to update the libxml flag.
LDLIBS += -L$(LOG4CPP_LIB) -llog4cpp
LDLIBS += -L/usr/lib64 -lxml2
LDLIBS += -L$(PYTHIA6) -lPythia6
# ROOT Flags are incomplete.
LDLIBS += -L$(ROOTSYS)/lib -lGeom -lEGPythia6
# If you add a new library, you may need to update this list.
LDLIBS += -L$(GENIE)/lib \
					-lGAlgorithm \
					-lGBaryonResonance \
					-lGBase \
					-lGBodekYang \
					-lGCharm \
					-lGCoh \
					-lGCrossSections \
					-lGDecay \
					-lGDfrc \
					-lGDIS \
					-lGElas \
					-lGElFF \
					-lGEVGCore \
					-lGEVGDrivers \
					-lGEVGModules \
					-lGFluxDrivers \
					-lGFragmentation \
					-lGGeo \
					-lGGiBUU \
					-lGHadronTransp \
					-lGHEP \
					-lGInteraction \
					-lGLlewellynSmith \
					-lGMEC \
					-lGMessenger \
					-lGMuELoss \
					-lGNtuple \
					-lGNuclear \
					-lGNuE \
					-lGNuGamma \
					-lGNumerical \
					-lGPDF \
					-lGPDG \
					-lGQEL \
					-lGQPM \
					-lGRegistry \
					-lGReinSeghal \
					-lGRES \
					-lGUtils 

all: transStudy diffXSec 

transStudy: transStudy.o 
	@if test -d ../bin; then echo BinOk >& /dev/null ; else mkdir -p ../bin; fi
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTFLAGS) -o transStudy transStudy.cxx $(LDLIBS) 
	mv transStudy ../bin

transStudy.o: transStudy.cxx
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTINC) -c transStudy.cxx -o transStudy.o 

diffXSec: diffXSec.o 
	@if test -d ../bin; then echo BinOk >& /dev/null ; else mkdir -p ../bin; fi
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTFLAGS) -o diffXSec diffXSec.cxx $(LDLIBS) 
	mv diffXSec ../bin

diffXSec.o: diffXSec.cxx
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTINC) -c diffXSec.cxx -o diffXSec.o 


clean:
	-rm -f ../bin/*
	-rm -f *.o

remake:
	make clean
	make all
