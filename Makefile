## Make file for CyMiniAna c++ scripts
ROOTCFLAGS = $(shell root-config --cflags)

SUFFIXES += .d

MYCF=""

ifneq ($(findstring m32,$(ROOTCFLAGS)),)
	MYCF= CXXFLAGS=-m32 CFLAGS=-m32
endif


CXXFLAGS  = -g -I. -I../lwtnn/include/ -I../CondFormats/
CXXFLAGS += -I../CondFormats/JetMetObjects/interface/
CXXFLAGS += -Wno-long-long -fPIC
CXXFLAGS += $(shell root-config --cflags)
CXXFLAGS += -I$(ALRB_BOOST_ROOT)/include
CXXFLAGS += -I$(SFT_HOME_eigen)/include/eigen3/

ROOTLIB   = $(shell root-config --glibs) -lMinuit -lTreePlayer

LDFLAGS   = $(ROOTLIB)

# base classes
OBJS  = Root/Event.o
OBJS += Root/eventSelection.o
OBJS += Root/miniTree.o
OBJS += Root/histogrammer.o
OBJS += Root/histogrammerTruth.o
OBJS += Root/efficiency.o
OBJS += Root/efficiencyTruth.o
OBJS += Root/configuration.o
OBJS += Root/tools.o
OBJS += Root/AMWT.o
OBJS += Root/MassSolver.o
OBJS += Root/TopMassVariables.o
OBJS += lwtnn/src/LightweightNeuralNetwork.o
OBJS += lwtnn/src/parse_json.o
OBJS += lwtnn/src/Stack.o

# base class steering macro
OBJSBASE += ${OBJS}
OBJSBASE += util/run.o

# skimming macro
OBJSSKIM += ${OBJS}
OBJSSKIM += util/skim.o

# AC
OBJSAC += ${OBJS}
OBJSAC += util/runAC.o


all: run skim runAC
# steering macro goes in 'all'. Separate by spaces
# e.g., 
# all: run runCustom1 runCustom2


NODEPS:=clean
SOURCES:=$(shell find Root/ util/ -name "*.cxx")
DEPFILES:=$(patsubst %.cxx,%.d,$(SOURCES))

ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
	-include $(DEPFILES)
endif


%.o: %.cxx %.d
	g++ -c $(CXXFLAGS) -o $@ $<

%.d: %.cxx
	g++ $(CXXFLAGS) -MM -MT '$(patsubst %.cxx,%.o,$<)' $< -MF $@

run: $(OBJSBASE)
	g++ $(CXXFLAGS) -o run $(OBJSBASE) $(LDFLAGS)

skim: $(OBJSSKIM)
	g++ $(CXXFLAGS) -o skim $(OBJSSKIM) $(LDFLAGS)

runAC: $(OBJSAC)
	g++ $(CXXFLAGS) -o runAC $(OBJSAC) $(LDFLAGS)

clean:
	rm -f lwtnn/src/*.o lwtnn/src/*.d Root/*.o Root/*.d util/*.o util/*.d run skim runAC runAC_deltay
