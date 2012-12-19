#-----------------------------------------------------------------------
# This Makefile is for building MOZART2 Pre-processor
#------------------------------------------------------------------------

# Set up special characters
null  :=
space := $(null) $(null)

# Check for directory in which to put executable
ifeq ($(MODEL_EXEDIR),$(null))
MODEL_EXEDIR := ../bin
endif

# Check for name of executable
ifeq ($(EXENAME),$(null))
EXENAME := campp
endif

# Check for source list
ifeq ($(SRCLIST),$(null))
SRCLIST := Base_Srclist_f
endif

ifeq ($(OBJ_DIR),$(null))
OBJ_DIR := ./OBJ
endif

# Load dependency search path.
# Check for source directories
ifeq ($(SRCDIRS),$(null))
dirs := ./cam_chempp
else
dirs := ./ $(SRCDIRS)
endif

# Determine platform 
UNAMES := $(shell uname -s)
UNAMEM := $(shell uname -m)

# Set cpp search path, include netcdf
cpp_dirs := $(dirs)
cpp_path := $(foreach dir,$(cpp_dirs),-I$(dir)) # format for command line

# Expand any tildes in directory names. Change spaces to colons.
VPATH    := $(foreach dir,$(cpp_dirs),$(dir)) 
VPATH    := $(subst $(space),:,$(VPATH))               

# Get list of files and determine objects and dependency files
base_srclist_f := $(shell cat $(SRCLIST))
OBJS      := $(foreach file,$(base_srclist_f),$(OBJ_DIR)/$(file:.f=.o))

all: $(MODEL_EXEDIR)/$(EXENAME)

# Architecture-specific flags and rules

#------------------------------------------------------------------------
# Altix (columbia)
#------------------------------------------------------------------------
ifeq ($(UNAMEM),ia64)
ifeq ($(COMPILER),ifort)

# @echo $(COMPILER)

 objflag    = -o $@

  FC          := ifort
  FFLAGS      := -c -132 -ftz -g -FR  -I $(OBJ_DIR)

  ifeq ($(DEBUG),TRUE)
    FFLAGS    += -CB
  else
    # Check for override of default Fortran compiler optimizations
    ifeq ($(F_OPTIMIZATION_OVERRIDE),$(null))
      FORTRAN_OPTIMIZATION := -O2
    endif
    FFLAGS    += $(FORTRAN_OPTIMIZATION)
  endif

.SUFFIXES:
.SUFFIXES: .f .F .c .o

$(OBJ_DIR)/%.o : %.f
	$(FC) $(FFLAGS) $(objflag) $<

endif
endif

#------------------------------------------------------------------------
# Alpha ?
#------------------------------------------------------------------------
ifeq ($(UNAMEM),alpha)

objflag    = -o $@
FC          := f95
CPP         := /usr/bin/cpp
FFLAGS      := -O4 -c -tune host -arch host -free -module $(OBJ_DIR) -I $(OBJ_DIR)


.SUFFIXES:
.SUFFIXES: .f .F .c .o

$(OBJ_DIR)/%.o : %.f
	$(FC) $(FFLAGS) $(objflag) $<

endif

#------------------------------------------------------------------------
# AIX ?
#------------------------------------------------------------------------
ifeq ($(UNAMES),AIX)

objflag    = -o $@
FC          := xlf95
FFLAGS      := -g -c -qarch=auto -qnosave -qfree=f90 -qmoddir=$(OBJ_DIR) -I $(OBJ_DIR) -qstrict -O3


.SUFFIXES:
.SUFFIXES: .f .F .c .o

$(OBJ_DIR)/%.o : %.f
	$(FC) $(FFLAGS) $(objflag) $<

endif

#------------------------------------------------------------------------
# BG/L, BG/P
#------------------------------------------------------------------------
ifeq ($(UNAMEM),ppc64)

objflag    = -o $@
FC          := xlf95
FFLAGS      := -g -c -qarch=auto -qnosave -qfree=f90 -qmoddir=$(OBJ_DIR) -I $(OBJ_DIR) -qstrict -O3


.SUFFIXES:
.SUFFIXES: .f .F .c .o

$(OBJ_DIR)/%.o : %.f
	$(FC) $(FFLAGS) $(objflag) $<

endif

#------------------------------------------------------------------------
# Intel Pentium ? ; use Lahey/Fujistu F95 compiler
#------------------------------------------------------------------------
ifeq ($(UNAMEM),i686)

obj_flag  = -o $@ 

ifeq ($(COMPILER),Lahey)
FC        := lf95
ifneq ($(DEBUG),ON)
FFLAGS    := --nfix --nap --nchk --ng --npca --nsav --ntrace --tp4 -c --mod $(OBJ_DIR) -O
else
FFLAGS    := --nfix --nap --nchk --g --npca --nsav --trace --trap --tp4 -c --mod $(OBJ_DIR) -O
endif
else
FC          := pgf95
FFLAGS      := -O2 -c -g -Mfree -module $(OBJ_DIR)
endif


.SUFFIXES:
.SUFFIXES: .f .F .c .o

$(OBJ_DIR)/%.o : %.f
	$(FC) $(FFLAGS) -o $@ $<

endif

#------------------------------------------------------------------------
# Linux - lightning
#------------------------------------------------------------------------
ifeq ($(UNAMES),Linux)
ifeq ($(UNAMEM),x86_64)

obj_flag  = -o $@ 

ifeq ($(COMPILER),Lahey)
FC        := lf95
ifneq ($(DEBUG),ON)
FFLAGS    := --nfix --nap --nchk --ng --npca --nsav --ntrace --tp4 -c --mod $(OBJ_DIR) -O
else
FFLAGS    := --nfix --nap --nchk --g --npca --nsav --trace --trap --tp4 -c --mod $(OBJ_DIR) -O
endif
else
FC          := pgf95
FFLAGS      := -O2 -c -g -Mfree -module $(OBJ_DIR)
endif


.SUFFIXES:
.SUFFIXES: .f .F .c .o

$(OBJ_DIR)/%.o : %.f
	$(FC) $(FFLAGS) -o $@ $<

endif
endif

#------------------------------------------------------------------------
# SGI ?
#------------------------------------------------------------------------
ifeq ($(UNAMES),IRIX64)

objflag    = -o $@
FC          := f90
FFLAGS      := -c -freeform -I $(OBJ_DIR) -O2 -OPT:Olimit=8200

.SUFFIXES:
.SUFFIXES: .f .F .c .o

$(OBJ_DIR)/%.o : %.f
	$(FC) $(FFLAGS) $(objflag) $<

endif

#------------------------------------------------------------------------
# Default rules
#------------------------------------------------------------------------

$(MODEL_EXEDIR)/$(EXENAME):	$(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

RM := rm

clean:
	$(RM) -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(MODEL_EXEDIR)/$(EXENAME)

realclean:
	$(RM) -f $(OBJ_DIR)/*.o *.d $(MODEL_EXEDIR)/$(EXENAME)
