
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

all: simpleAna

simpleAna: simpleAna.o 
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTFLAGS) -o simpleAna simpleAna.cxx $(LDLIBS) 

simpleAna.o: simpleAna.cxx
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTINC) -c simpleAna.cxx -o simpleAna.o 

clean:
	-rm -f simpleAna
	-rm -f simpleAna.o

remake:
	make clean
	make all
