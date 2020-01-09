FC := gfortran
SRCDIR := ../src
BINDIR := ../bin
MODDIR := ../mod

all: test_legolas

libPFUnit.a: $(BINDIR)/*.o
	$(AR) -r $@ $?

ifeq (nagfor,$(findstring nagfor,$(FC)))
  FFLAGS += -fpp
endif

%.o : %.F90
	$(FC) -c $(FFLAGS) -I$(MODDIR) -I$(BINDIR) -I$(SRCDIR) $<

LATEST_PFUNIT_DIR := $(lastword $(shell echo $(wildcard $(PFUNIT_DIR)/PFUNIT-4.*) | xargs -n1 | sort -V))
include $(LATEST_PFUNIT_DIR)/include/PFUNIT.mk

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)

test_legolas: libPFUnit.a

test_legolas_TESTS := core_tests.pf
test_legolas_OTHER_LIBRARIES := -L. -lPFUnit -lblas -llapack
$(eval $(call make_pfunit_test,test_legolas))


clean:
	$(RM) *.o *.mod *.a *.inc *.F90 test_legolas
