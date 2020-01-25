FC := gfortran
SRCDIR := ../src
BINDIR := ../bin
MODDIR := ../mod

all:	list_builtins \
			test_legolas
	
list_builtins: 
	$(FC) -E -dM - < /dev/null > builtins.inc
		
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

test_legolas: libpFUnit.a

test_legolas_TESTS := core_tests.pf comparison_tests.pf
test_legolas_OTHER_LIBRARIES := -L. -lpFUnit -lblas -llapack
$(eval $(call make_pfunit_test,test_legolas))


clean:
	$(RM) *.o *.mod *.a *.inc *.F90 test_legolas
