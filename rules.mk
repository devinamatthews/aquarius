all: ctf $(DEFAULT_COMPONENTS)

ALL_COMPONENTS = libs bins bench_ao_ccsd bench_ao_ccsd_lambda bench_ao_ccsdt \
                 bench_cholesky_ccsd bench_cholesky_ccsd_lambda bench_cholesky_ccsdt

libs: phase1

bins bench_ao_ccsd bench_ao_ccsd_lambda bench_ao_ccsdt bench_cholesky_ccsd \
     bench_cholesky_ccsd_lambda bench_cholesky_ccsdt: libs phase2

bins: bench_ao_ccsd bench_ao_ccsd_lambda bench_ao_ccsdt bench_cholesky_ccsd \
      bench_cholesky_ccsd_lambda bench_cholesky_ccsdt

bindir = $(topdir)/bin
libdir = $(topdir)/lib

ALL_LIBS_LINK = $(wildcard $(topdir)/src/autocc/*.o) \
                $(wildcard $(topdir)/src/input/*.o) \
                $(wildcard $(topdir)/src/time/*.o) \
                $(wildcard $(topdir)/src/integrals/*.o) \
                $(wildcard $(topdir)/src/memory/*.o) \
                $(wildcard $(topdir)/src/tensor/*.o) \
                $(wildcard $(topdir)/src/util/*.o) \
                $(wildcard $(topdir)/src/cc/*.o) \
                $(wildcard $(topdir)/src/scf/*.o) \
                $(wildcard $(topdir)/src/symmetry/*.o) \
                $(wildcard $(topdir)/src/operator/*.o) \
                $(wildcard $(topdir)/src/task/*.o)

VPATH=$(srcdir)$(subst $(topdir),,$(shell pwd))

DEPDIR = .deps
DEPS += $(topdir)/.dummy $(addprefix $(DEPDIR)/,$(notdir $(patsubst %.o,%.Po,$(wildcard *.o))))
ALL_SUBDIRS = $(sort $(SUBDIRS) $(foreach comp,$(ALL_COMPONENTS),$(value $(addsuffix _SUBDIRS,$(comp)))))

_CPPFLAGS = $(CPPFLAGS)
_DEFS = $(DEFS)
_LDFLAGS = $(LDFLAGS) -L$(topdir)/lib
_INCLUDES = $(INCLUDES) -I. -I$(topdir) -I$(srcdir) -I$(srcdir)/src -I$(CTFDIR)/include #-I$(ELEMENTAL)/include
_CFLAGS = $(OPT) $(WARN) $(CFLAGS)
_CXXFLAGS = $(OPT) $(WARN) $(CXXFLAGS)
#_F77FLAGS = $(F77FLAGS)
#_F90FLAGS = $(F90FLAGS)
_DEPENDENCIES = $(DEPENDENCIES) Makefile $(topdir)/config.mk $(topdir)/rules.mk $(topdir)/config.h
_LIBS = $(LIBS) $(CTF_LIBS) $(ELEMENTAL_LIBS) $(BLAS_LIBS)

#F77COMPILE = $(F77) $(_INCLUDES) $(_F77FLAGS)
#F90COMPILE = $(F90) $(_INCLUDES) $(_F90FLAGS)
CCOMPILE = $(CC) $(_DEFS) $(_INCLUDES) $(_CPPFLAGS) $(_CFLAGS)
CXXCOMPILE = $(CXX) $(_DEFS) $(_INCLUDES) $(_CPPFLAGS) $(_CXXFLAGS)
ifeq ($(DEPSTYPE),normal)
    CCOMPILEDEPS = $(CCOMPILE) $(DEPFLAGS)
    CXXCOMPILEDEPS = $(CXXCOMPILE) $(DEPFLAGS)
else
    CCOMPILEDEPS = $(CCOMPILE) $(DEPFLAGS) > $(DEPDIR)/$(notdir $*).Po; $(CCOMPILE)
    CXXCOMPILEDEPS = $(CXXCOMPILE) $(DEPFLAGS) > $(DEPDIR)/$(notdir $*).Po; $(CXXCOMPILE)
endif

LINK = $(CXX) $(_CXXFLAGS) $(_LDFLAGS) -o $@
ARCHIVE = $(AR) $@

#
# Automatic library dependency generation derived from: Steve Dieters
# http://lists.gnu.org/archive/html/help-make/2010-03/msg00072.html
#
libdeps = \
$(foreach d, $(patsubst -L%,%,$(filter -L%,$(1))),\
 $(foreach l, $(patsubst -l%,%,$(filter -l%,$(1))),\
  $(if $(shell if [ -e $(d)/lib$(l).a ]; then echo "X"; fi),\
   $(d)/lib$(l).a\
  )\
 )\
) $(filter %.a,$(1))

#.NOTPARALLEL:
.PHONY: all default clean $(ALL_COMPONENTS) ctf phase1 phase2
FORCE:

ctf: FORCE
	@if [ -d ctf -a $(make_ctf) = yes ];then \
		echo "Making ctf in ctf"; \
		(cd ctf && $(MAKE) ctf); \
	fi
	
phase1: ctf
	@if [ -f config.mk -o "x$(MAKECMDGOALS)" = "xlibs" ]; then \
		for dir in $(ALL_SUBDIRS); do \
			echo "Compiling in $$dir"; \
			(cd $$dir && $(MAKE) libs); \
		done; \
	fi

phase2: phase1 ctf
	@for dir in $(ALL_SUBDIRS); do \
		echo "Linking $(MAKECMDGOALS) in $$dir"; \
		(cd $$dir && $(MAKE) $(MAKECMDGOALS)); \
	done

clean:
	rm -rf $(DEPDIR) *.o
	@if [ -d ctf -a $(make_ctf) = yes ];then \
    	echo "Cleaning in ctf"; \
    	(cd ctf && $(MAKE) clean); \
	fi; \
	for subdir in $(ALL_SUBDIRS); do \
		echo "Cleaning in $$subdir"; \
		(cd $$subdir && $(MAKE) clean); \
	done

$(bindir)/%: $(_DEPENDENCIES) $(call libdeps,$(_LIBS)) $(ALL_LIBS_LINK)
	@mkdir -p $(dir $@)
	$(LINK) $(filter %.o,$^) $(FLIBS) $(_LIBS)

$(libdir)/%: $(_DEPENDENCIES)
#	@mkdir -p $(dir $@)
#	$(ARCHIVE) $(filter %.o,$^)

%.o: %.f $(_DEPENDENCIES)
	$(F77COMPILE) -c -o $@ $<

%.o: %.f90 $(_DEPENDENCIES)
	$(F90COMPILE) -c -o $@ $<

%.o: %.c $(_DEPENDENCIES)
	@mkdir -p $(DEPDIR)
	$(CCOMPILEDEPS) -c -o $@ $<

%.o: %.cxx $(_DEPENDENCIES)
	@mkdir -p $(DEPDIR)
	$(CXXCOMPILEDEPS) -c -o $@ $<

-include $(DEPS)
