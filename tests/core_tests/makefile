FC := gfortran
SRCDIR := ../../src
BINDIR := ../../bin
MODDIR := ../../mod

all:	test_legolas_core

libpFUnit.a: $(BINDIR)/*.o
	$(AR) -r $@ $?

ifeq (nagfor,$(findstring nagfor,$(FC)))
  FFLAGS += -fpp
endif
FFLAGS += -cpp

%.o : %.F90
	$(FC) -c $(FFLAGS) -I$(MODDIR) -I$(BINDIR) -I$(SRCDIR) $<

LATEST_PFUNIT_DIR := $(lastword $(shell echo $(wildcard $(PFUNIT_DIR)/PFUNIT-4.*) | xargs -n1 | sort -V))
include $(LATEST_PFUNIT_DIR)/include/PFUNIT.mk

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)

test_legolas_core: libpFUnit.a

test_legolas_core_TESTS := core_tests.pf
test_legolas_core_OTHER_LIBRARIES := -L. -lpFUnit -lblas -llapack
$(eval $(call make_pfunit_test,test_legolas_core))

clean:
	$(RM) *.o *.mod *.a *.inc *.F90 test_legolas_core
