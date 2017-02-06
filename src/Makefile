
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

all: diffractiveStudy countPi0 diffXSec gstana_r grooana_r qetotalxsec

diffractiveStudy: diffractiveStudy.o 
	@if test -d ../bin; then echo BinOk >& /dev/null ; else mkdir -p ../bin; fi
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTFLAGS) -o diffractiveStudy diffractiveStudy.cxx $(LDLIBS) 
	mv diffractiveStudy ../bin

diffractiveStudy.o: diffractiveStudy.cxx
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTINC) -c diffractiveStudy.cxx -o diffractiveStudy.o 

countPi0: countPi0.o 
	@if test -d ../bin; then echo BinOk >& /dev/null ; else mkdir -p ../bin; fi
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTFLAGS) -o countPi0 countPi0.cxx $(LDLIBS) 
	mv countPi0 ../bin

countPi0.o: countPi0.cxx
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTINC) -c countPi0.cxx -o countPi0.o 

diffXSec: diffXSec.o 
	@if test -d ../bin; then echo BinOk >& /dev/null ; else mkdir -p ../bin; fi
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTFLAGS) -o diffXSec diffXSec.cxx $(LDLIBS) 
	mv diffXSec ../bin

diffXSec.o: diffXSec.cxx
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTINC) -c diffXSec.cxx -o diffXSec.o 

gstana_r: Rungst.cxx gstana.cxx gstana.h
	@if test -d ../bin; then echo BinOk >& /dev/null ; else mkdir -p ../bin; fi
	$(CXX) -o gstana_r Rungst.cxx gstana.cxx $(ROOTFLAGS)
	mv gstana_r ../bin

grooana_r: Rungroo.cxx grooana.cxx grooana.h
	@if test -d ../bin; then echo BinOk >& /dev/null ; else mkdir -p ../bin; fi
	$(CXX) -o grooana_r Rungroo.cxx grooana.cxx $(ROOTFLAGS)
	mv grooana_r ../bin

qetotalxsec: qetotalxsec.o 
	@if test -d ../bin; then echo BinOk >& /dev/null ; else mkdir -p ../bin; fi
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTFLAGS) -o qetotalxsec qetotalxsec.cxx $(LDLIBS) 
	mv qetotalxsec ../bin

qetotalxsec.o: qetotalxsec.cxx
	$(CXX) $(CPPFLAGS) $(INCLUDE) $(ROOTINC) -c qetotalxsec.cxx -o qetotalxsec.o 


clean:
	-rm -f ../bin/diffractiveStudy
	-rm -f ../bin/countPi0
	-rm -f ../bin/diffXSec
	-rm -f *.o
	-rm -f ../bin/gstana_r
	-rm -f ../bin/grooana_r

remake:
	make clean
	make all
