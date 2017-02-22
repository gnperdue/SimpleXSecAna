
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
					-lGReinSehgal \
					-lGRES \
					-lGUtils 

all: qetotalxsec qelikediff

qetotalxsec: qetotalxsec.o 
	@if test -d ../bin; then echo BinOk >& /dev/null ; else mkdir -p ../bin; fi
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTFLAGS) -o qetotalxsec qetotalxsec.cxx $(LDLIBS) 
	mv qetotalxsec ../bin

qetotalxsec.o: qetotalxsec.cxx
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTINC) -c qetotalxsec.cxx -o qetotalxsec.o 

qelikediff: qelikediff.o 
	@if test -d ../bin; then echo BinOk >& /dev/null ; else mkdir -p ../bin; fi
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTFLAGS) -o qelikediff qelikediff.cxx $(LDLIBS) 
	mv qelikediff ../bin

qelikediff.o: qelikediff.cxx
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTINC) -c qelikediff.cxx -o qelikediff.o 

clean:
	-rm -f *.o
	-rm -f ../bin/qetotalxsec
	-rm -f ../bin/qelikediff

remake:
	make clean
	make all
