# Makefile for unilib
# subversion Id for THIS file : $Id: Makefile 35 2019-03-14 15:46:47Z dhconsultancy@skynet.be $
# $URL: https://www.mag-unilib.eu/svn/branches/VALIRENE/Makefile $
#-----------------------------------------------------------------------

NAME	= unilib

# install directory
# The following is overridden in the local version of UNILIB so that it installs to ./lib
#prefix = /usr/local
prefix = .
exec_prefix = $(prefix)
libdir = $(exec_prefix)/lib
includedir = $(prefix)/include

#-------------------
# expected directory structure:
# .
#|-- doc
#|-- include
#|-- python_examples
#|-- src
#`-- test

#-------------------
# select compiler and flags
#	-- gfortran
FC	= gfortran
FFLAGS = -ffree-line-length-0 -O2
#
#	-- ifort
#FC     = ifort
#FFLAGS = -O2
#
#	-- pgf
#FC	= pgf90
#FFLAGS = -O2

#-------------------
# select linker flags
LDFLAGS = -Bstatic

#-------------------
LD	= $(FC)
RM	= rm -f

SRCDIR	= ./src
SRC_fixed	=  \
    $(SRCDIR)/ext530.f \
    $(SRCDIR)/ext630.f \
    $(SRCDIR)/ext635.f \
    $(SRCDIR)/ua610.f \
    $(SRCDIR)/ua630.f \
    $(SRCDIR)/ud310.f \
    $(SRCDIR)/ud320.f \
    $(SRCDIR)/ud330.f \
    $(SRCDIR)/uf410.f \
    $(SRCDIR)/uf420.f \
    $(SRCDIR)/ul220.f \
    $(SRCDIR)/ul230.f \
    $(SRCDIR)/ul240.f \
    $(SRCDIR)/um510.f \
    $(SRCDIR)/um520.f \
    $(SRCDIR)/um530.f \
    $(SRCDIR)/ut540.f \
    $(SRCDIR)/ut550.f \
    $(SRCDIR)/ut980.f \
    $(SRCDIR)/ut990.f \
    $(SRCDIR)/uxidl.f \
    $(SRCDIR)/IRBEMAPI.f
SRC_free	= \
    $(SRCDIR)/version.f90
SRCS	= $(SRC_fixed) $(SRC_free)
SRCINCLUDE	= $(SRCDIR)/structure.h \
              $(SRCDIR)/IRBEMAPI.inc
OBJS	= \
    $(notdir $(SRC_fixed:.f=.o)) \
    $(notdir $(SRC_free:.f90=.o))
ifeq ($(OS), Windows_NT)
    LIBS = lib$(NAME).dll lib$(NAME).a
else
    LIBS = lib$(NAME).so lib$(NAME).a
endif

TESTDIR = ./test
TESTSRCS	= \
    $(TESTDIR)/ex1.f \
    $(TESTDIR)/ex2.f \
    $(TESTDIR)/ex3.f \
    $(TESTDIR)/ex4.f \
    $(TESTDIR)/ex5.f \
    $(TESTDIR)/testIRBEMAPI.f
TESTOBJS	= $(notdir $(TESTSRCS:.f=.o))
TESTS	= $(TESTOBJS:.o=.exe)
TESTNMLS	= ex5.nml
TESTS_OUT	= $(TESTS:.exe=.out)
TESTS_DIFF	= $(TESTS_OUT:.out=.diff)
ALL_TESTS	= $(TESTOBJS) $(TESTS) $(TESTNMLS) $(TESTS_OUT) $(TESTS_DIFF)

vpath %.f	$(SRCDIR):$(TESTDIR)
vpath %.f90	$(SRCDIR)

FFLAGS	+= -I $(SRCDIR)
LDFLAGS += -L .
LDLIBS	+= -l$(NAME)

LIBDIR    := ./lib
INCDIR    := ./include

#PYNEXE    := ${ANACONDA3_BIN}/python
PYNEXE    := python
PYNSRCDIR := ${SRCDIR}
PYNSRC    := $(PYNSRCDIR)/PyUNILIB.pyx
PYNHDR    := $(INCDIR)/PyUNILIB_Interface.h
PYNLIB    := $(LIBDIR)/PyUNILIB.*.so

.SUFFIXES:
.PHONY: lib test clean install reftest dist distclean

# create unilib library (static and dynamic versions)
lib: $(LIBS)

test: $(TESTNMLS) $(TESTS_OUT)

clean:
	$(RM) $(OBJS) $(LIBS) $(ALL_TESTS)

install: $(LIBS) $(SRCINCLUDE) $(PYNLIB)
	mkdir -p $(libdir) $(includedir)
	install -p -m 644 $(LIBS) $(DESTDIR)$(libdir)
	install -p -m 644 $(SRCINCLUDE) $(DESTDIR)$(includedir)

$(PYNLIB) : $(PYNSRC) $(PYNHDR)
	rm -f -r build/
	rm -f $(LIBDIR)/PyUNILIB.*.so
	$(PYNEXE) setup.py build_ext --inplace

reftest: $(TESTS_OUT)
	@echo "*WARNING* : replacing reference output in test directory"
	cp -p $(TESTS_OUT) $(TESTDIR)/

lib$(NAME).a: $(OBJS)
	ar crs $@ $^
	$(RM) $(OBJS)

lib$(NAME).so: $(SRCS)
	$(FC) $(FFLAGS) -shared -o $@ -fPIC $^

lib$(NAME).dll: $(SRCS)
	$(FC) $(FFLAGS) -shared -o $@ $^

# ---- generic rules for objects

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# ---- generic rules for tests

%.exe: %.o $(LIBS)
	$(LD) $(LDFLAGS) $< $(LDLIBS) -o $@

%.out: %.exe
	./$< > $@
	@diff -EZbB --strip-trailing-cr $@ $(TESTDIR)/$@ > $*.diff || echo "  *** COMPARISON FAILED, see "$*.diff

%.nml: $(TESTDIR)/%.nml
	cp -p $< ./$@
