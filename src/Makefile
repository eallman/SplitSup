# CONFIG

# Choose: mac, windows, or unix
ARCHITECTURE ?= unix

# Set MPI=yes for the parallel version
MPI ?= no

# Set DEBUG=yes for a debug version of mrbayes, otherwise set OPTFLAGS
# to the desired optimization level (e.g. -O2 -fomit-frame-pointer for
# less aggressive optimization than the default -O3)
#
# See the MrBayes manual for some optimization hints for your particular
# architecture or consult the manual of your compiler, eg. "info gcc".

DEBUG ?= yes
#OPTFLAGS = -O2 -march=pentium4 -mfpmath=sse -fomit-frame-pointers
OPTFLAGS ?= -O3

# set compiler for the non-MPI version (mpicc will be used for the MPI
# version regardless of this setting)
CC = gcc

# set to yes if you want to use the readline library (make sure you have it
# installed on your machine)
# defaults to yes for unix machines
ifeq ($(strip $(ARCHITECTURE)), unix)
   USEREADLINE ?= yes
else
   USEREADLINE ?= no
endif

#
# End of user configuration

ifeq ($(strip $(ARCHITECTURE)),mac)
        CFLAGS += -DMAC_VERSION
else
ifeq ($(strip $(ARCHITECTURE)), windows)
        CFLAGS += -DWIN_VERSION
else
ifeq ($(strip $(ARCHITECTURE)), unix)
        CFLAGS += -DUNIX_VERSION
else
        ARCHITECTURE = none
endif
endif
endif

ifeq  ($(strip $(USEREADLINE)),yes)
    CFLAGS += -DUSE_READLINE
        LIBS += -lncurses -lreadline
endif

ifeq  ($(strip $(MPI)),yes)
        CFLAGS += -DMPI_ENABLED
        CC = mpicc
endif

ifeq ($(strip $(DEBUG)), yes)
        CFLAGS += -ggdb
else
        CFLAGS += $(OPTFLAGS)
endif

CFLAGS   +=   -Wall

INCL_GSL = /opt/local/include
INCLUDE_FLAGS += -I$(INCL_GSL)

LIBS     +=   -lm -lgsl -lgslcblas

GSLLIBDIR = /opt/local/lib/
LBFLAGS  +=  -L$(GSLLIBDIR)

LDLIBS   +=  $(LBFLAGS) $(LIBS)
#LDLIBS   += $(LIBS)

OBJECTS  =   main.o flatBinary.o 

PROGS    = SplitSup

ifeq ($(ARCHITECTURE), none)
missing:
        @echo
        @echo Please set compilation target in this Makefile.
        @echo set ARCHITECTURE to one of \"mac\", \"windows\" or \"unix\"
        @echo set MPI to yes if you want to use the parallel version
        @echo set DEBUG to generate a debug version of MrBayes
        @echo
endif

all:            $(PROGS)

clean:
		rm -f *.o

SplitSup:	$(OBJECTS)
		$(CC) $(CFLAGS) $(OBJECTS) $(LDFLAGS) $(LDLIBS) -o SplitSup

main.o:         main.c
		$(CC) $(CFLAGS) $(INCLUDE_FLAGS) -c main.c

flatBinary.o:	flatBinary.c
		$(CC) $(CFLAGS) $(INCLUDE_FLAGS) -c flatBinary.c

