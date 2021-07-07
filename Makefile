all: replay_gmn
#ReplayMCProd

#DEBUG = 1

ifdef DEBUG
  CXXFLAGS    = -g -O0 -Wall -Wextra
else
  CXXFLAGS    = -g -O2 -Wall -Wextra
endif

ifndef ANALYZER
  $(error $$ANALYZER environment variable not defined)
endif

INCDIRS  = $(wildcard $(addprefix $(ANALYZER)/, include src hana_decode))
#INCDIRS     = ../src ../hana_decode
INCDIRS  +=${LIBSBSDIG}/src
INCDIRS  +=${SBS}/include
INCLUDES    = $(addprefix -I, $(INCDIRS))

CXX         = $(shell root-config --cxx)
ROOTCFLAGS  = $(shell root-config --cflags)
LD          = $(shell root-config --ld)
LDFLAGS     = $(shell root-config --ldflags)

CXXFLAGS   += $(ROOTCFLAGS) $(INCLUDES)
#ALLINCLUDES  = -I$(shell root-config --incdir) $(INCLUDES)
LIBS        = $(shell root-config --libs) -L${ANALYZER}/lib64 -lHallA -ldc -lPodd 

# Add EVIO lib, needed by libdc
ifndef EVIO_LIBDIR
  EVIO_LIBDIR = ${ANALYZER}/lib64
endif
LIBS += -L$(EVIO_LIBDIR) -levio

ifndef SBS_LIBDIR
  SBS_LIBDIR = ${SBS}/lib64
endif
LIBS += -L$(SBS_LIBDIR) -lsbs

replay_gmn:	replay_gmn.o
		$(LD) $(LDFLAGS) $(LIBS) -o $@ $^

clean:
		rm -f replay_gmn replay_gmn.o

%.o:		%.cxx Makefile
		$(CXX) $(CXXFLAGS) -o $@ -c $<

.PHONY: all clean
