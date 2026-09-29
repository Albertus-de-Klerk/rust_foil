#*********************************************************
# XFOIL MinGW Build - Legacy Precision (Match Standard Binary)
#*********************************************************

SHELL = sh
BINDIR = .
PROGS = qfoil pplot pxplot

SRC = ../src
OSRC = ../osrc
VPATH = $(SRC):$(OSRC)

FC = gfortran
CC = gcc

#----------------------------------------------------
# FLAGS EXPLANATION FOR MATCHING STANDARD BINARY:
# 1. REMOVED -fdefault-real-8: Reverts to Single Precision (32-bit).
#    This matches the rounding errors of the original 90s binaries.
# 2. -O2: Optimization is required to match standard rounding behavior.
# 3. -fno-automatic: Critical for solver state retention.
# 4. -finit-local-zero: Initializes variables to 0.0 (fixes random glitches).
#----------------------------------------------------

FLAGS_LEGACY = -O2 -std=legacy -fno-automatic -finit-local-zero -fno-align-commons -fdefault-real-8 -mtune=generic -flto
CFLAGS  = -O2 -mtune=generic -flto -DUNDERSCORE
FFLAGS  = $(FLAGS_LEGACY)
FFLOPT  = $(FLAGS_LEGACY)

# Static Linking to avoid DLL errors
LDFLAGS = -static -s -flto

# Plot Libraries
PLTOBJ = ../plotlib/libPlt.a
PLTLIB = 
FTNLIB =

#----------------------------------------------------
# OBJECT LISTS
#----------------------------------------------------
XFOILOBJ = xfoil.o xpanel.o xoper.o xtcam.o xgdes.o xqdes.o xmdes.o \
xsolve.o xbl.o xblsys.o xpol.o xplots.o pntops.o xgeom.o xutils.o modify.o \
blplot.o polplt.o aread.o naca.o spline.o plutil.o iopol.o gui.o sort.o \
dplot.o profil.o

PPLOTOBJ = pplot.o polplt.o sort.o iopol.o
PXPLOTOBJ = pxplot.o plutil.o gui.o
XUTILOBJ  = userio.o
OSOBJ = frplot.o ntcalc.o osmap.o getosfile.o

#----------------------------------------------------
# BUILD TARGETS
#----------------------------------------------------

all: $(PROGS)

qfoil: $(XFOILOBJ) $(XUTILOBJ) $(OSOBJ)
	$(FC) -o qfoil $(XFOILOBJ) $(XUTILOBJ) $(OSOBJ) $(PLTOBJ) $(PLTLIB) $(LDFLAGS)

pxplot: $(PXPLOTOBJ) $(XUTILOBJ)
	$(FC) -o pxplot $(PXPLOTOBJ) $(XUTILOBJ) $(PLTOBJ) $(PLTLIB) $(LDFLAGS)

pplot: $(PPLOTOBJ) $(XUTILOBJ)
	$(FC) -o pplot $(PPLOTOBJ) $(XUTILOBJ) $(PLTOBJ) $(PLTLIB) $(LDFLAGS)

blu: blu.o profil.o
	$(FC) -o blu blu.o profil.o $(LDFLAGS)

#----------------------------------------------------
# GENERIC RULES
#----------------------------------------------------
%.o: %.f
	$(FC) -c $(FFLAGS) $<

%.o: %.c
	$(CC) -c $(CFLAGS) $<

.PHONY: clean

# Executable suffix on Windows
EXEEXT :=
ifeq ($(OS),Windows_NT)
  EXEEXT := .exe
endif

# Pick a remove command:
# - If rm exists (MSYS2/Git Bash/Cygwin), use it
# - Otherwise (plain cmd.exe), fall back to del
ifeq ($(OS),Windows_NT)
  RMTEST := $(strip $(shell where rm 2>nul))
  ifneq ($(RMTEST),)
    RM := rm -f
  else
    RM := cmd /C del /Q
  endif
else
  RM := rm -f
endif

clean:
ifeq ($(OS),Windows_NT)
	-$(RM) *.o 2>nul
	-$(RM) qfoil$(EXEEXT) pplot$(EXEEXT) pxplot$(EXEEXT) 2>nul
else
	$(RM) *.o qfoil$(EXEEXT) pplot$(EXEEXT) pxplot$(EXEEXT)
endif