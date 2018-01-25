###################################################################
# This Makefile was created using the CreateProject.sh script
# for project GPIIBackgroundAlpha.
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

# BAT
CXXFLAGS += $(shell bat-config --cflags) -lBATmvc
LIBS := $(shell bat-config --libs)

# CUBA
#CXXFLAGS += -I$(CUBA_BASE_DIR)/include
#LIBS += -L$(CUBA_BASE_DIR)/lib -lcuba

# ROOT
CXXFLAGS += $(shell root-config --cflags)
LIBS += $(shell root-config --libs) -lMinuit

# GERDA-ADA
CXXFLAGS += $(shell gerda-ada-config --cflags)/gerda-ada/
LIBS += $(shell gerda-ada-config --libs)

# MGDO
CXXFLAGS += $(shell mgdo-config --cflags)
LIBS += $(shell mgdo-config --libs)

# GELATIO
CXXFLAGS += $(shell gelatio-config --cflags)
LIBS += $(shell gelatio-config --libs)

# ProgressBar
CXXFLAGS += -I$(PROGRESS_BAR_DIR)
LIBS += -L$(PROGRESS_BAR_DIR) -lProgressBar

# jsoncpp
CXXFLAGS += -I$(JSONCPP_BASE_DIR)/include/
LIBS += -L$(JSONCPP_BASE_DIR)/lib -ljsoncpp

# List of all classes (models) used in the program
# Add classes to the end. Backslash indicates continuation
# on the next line
CXXSRCS      = \
    GPIIBackgroundAlpha.cxx

# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#

CXXOBJS      = $(patsubst %.cxx,%.o,$(CXXSRCS))
MYPROGS      = runGPIIBackgroundAlpha

GARBAGE      = $(CXXOBJS) *.o *~ link.d $(MYPROGS)


# targets
all : runGPIIBackgroundAlpha

link.d : $(patsubst %.cxx,%.h,$(CXXSRCS))
    $(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d;

-include link.d

%.o : %.cxx
    $(CXX) $(CXXFLAGS) -c $< -o $@

clean :
    $(RM) $(GARBAGE)

runGPIIBackgroundAlpha : runGPIIBackgroundAlpha.cxx $(CXXOBJS)
    $(CXX) $(CXXFLAGS) -c $<
    $(CXX) $(LDFLAGS) $(LIBS) runGPIIBackgroundAlpha.o $(CXXOBJS) -o runGPIIBackgroundAlpha

print :
   echo compiler  : $(CXX)
   echo c++ srcs  : $(CXXSRCS)
   echo c++ objs  : $(CXXOBJS)
   echo c++ flags : $(CXXFLAGS)
   echo libs      : $(LIBS)
   echo so flags  : $(SOFLAGS)

   echo rootlibs  : $(ROOTLIBS)
   echo rootglibs : $(ROOTGLIBS)
