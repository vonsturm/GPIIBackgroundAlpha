###################################################################
# This Makefile was created using the CreateProject.sh script
# for project BEGE_backgrounds.
# CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
# BAT can be downloaded from http://www.mppmu.mpg.de/bat
###################################################################
#
# Run 'make' to compile the program and 'make clean' to remove
# all compiled parts and 'clean' the directory.
#
# You might need to adjust the CFLAGS, LIBS, and GLIBS based on
# the BAT installation on your system. Consult the gmake manual
# for details.
#
###################################################################

# Root variables
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs) -lMinuit
ROOTGLIBS    := $(shell root-config --glibs)

#ROOTLIBS     += -lRooFitCore -lRooFit -lRooStats -lFoam -lMathMore

# compiler and flags
CXX          = g++ 
CXXFLAGS     = -g -Wall -fPIC -Wno-deprecated -O2
LD           = g++ 
LDFLAGS      = -g -O2
SOFLAGS      = -shared

# standard commands
RM           = rm -f
MV           = mv
ECHO         = echo
CINT         = rootcint

# add ROOT flags
CXXFLAGS    += $(ROOTCFLAGS) -I$(CUBA_BASE_DIR)/include

LIBS        += -L$(CUBA_BASE_DIR)/lib -lcuba
GLIBS       += -L$(CUBA_BASE_DIR)/lib -lcuba


# add GERDA-ADA flags
CXXFLAGS    += -I$(GERDA_BASE_DIR)/include/gerda-ada/

LIBS        += -L$(GERDA_BASE_DIR)/lib -lgerda-ada
GLIBS       += -L$(GERDA_BASE_DIR)/lib -lgerda-ada

# add GELATIO flags
CXXFLAGS    += -I$(GERDA_BASE_DIR)/include/gelatio

LIBS        += -L$(GERDA_BASE_DIR)/lib -lGELATIODecoders -lGELATIOManagement -lModules -lGELATIOUtilities
GLIBS       += -L$(GERDA_BASE_DIR)/lib -lGELATIODecoders -lGELATIOManagement -lModules -lGELATIOUtilities


# ----------------------------------------------------------------------
# The following definitions depend on the setup of the system where
# the project is being compiled. If BAT is installed in the standard
# system search path or the installation directory is defined in the
# BATINSTALLDIR environmental variable then the lines below are correct
# and the compilation will work
CXXFLAGS    += -I. -I./include -I$(shell bat-config --prefix)/include
LIBS        += $(ROOTLIBS)  -L$(shell bat-config --prefix)/lib -lBAT -lBATmodels
GLIBS       += $(ROOTGLIBS) -L$(shell bat-config --prefix)/lib -lBAT -lBATmodels

# List of all classes (models) used in the program
# Add classes to the end. Backslash indicates continuation
# on the next line
CXXSRCS      = \
        BEGE_backgrounds.cxx RunPhaseII.cxx DetectorPhaseII.cxx

# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#

CXXOBJS      = $(patsubst %.cxx,%.o,$(CXXSRCS))
EXEOBJS      =
MYPROGS     = \
        runBEGE_backgrounds test_classes

GARBAGE      = $(CXXOBJS) $(EXEOBJS) *.o *~ link.d $(MYPROGS)


# targets
all : project project_test

link.d : $(patsubst %.cxx,%.h,$(CXXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d;

-include link.d

%.o : %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	$(RM) $(GARBAGE)

project : runBEGE_backgrounds.cxx $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(LIBS) runBEGE_backgrounds.o $(CXXOBJS) -o runBEGE_backgrounds

project_test : test_classes.cxx $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(LIBS) test_classes.o $(CXXOBJS) -o test_classes

print :
   echo compiler  : $(CXX)
   echo c++ srcs  : $(CXXSRCS)
   echo c++ objs  : $(CXXOBJS)
   echo c++ flags : $(CXXFLAGS)
   echo libs      : $(LIBS)
   echo so flags  : $(SOFLAGS)

   echo rootlibs  : $(ROOTLIBS)
   echo rootglibs : $(ROOTGLIBS)

