ARCH = linux

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
ROOTGLIBS := $(shell root-config --glibs)

ifeq ($(ARCH),linux)
CXX = g++
CXXFLAGS = -O -Wno-deprecated -fPIC -m32 -fno-inline -Wno-write-strings
LD = g++
LDFLAGS = -O
SOFLAGS = -shared
INCLUDE = $(ROOTSYS)/include/
endif

CXXFLAGS += $(ROOTCFLAGS)
LIBS = $(ROOTLIBS) $(SYSLIBS)
GLIBS = $(ROOTGLIBS) $(SYSLIBS)

#-----------------------------------------------

CXXFILES = RunData12.cxx
HEADERS = RunData12.h

#-----------------------------------------------

All: RunData12.so

RunData12.so: LinkDef.h $(HEADERS) $(CXXFILES)
	@echo "generating dictionary $@..."
	@rootcint -f RunData12Dict.cxx -c $(HEADERS) LinkDef.h
	g++ $(CXXFILES) RunData12Dict.cxx -shared -o RunData12.so $(CXXFLAGS) $(GLIBS)

clean:
	@rm RunData12.so
