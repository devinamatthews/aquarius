all: $(DEFAULT_COMPONENTS)

ALL_COMPONENTS = aquarius slide scf autocc input time cc memory

slide: memory
cc: scf input autocc time memory
input: slide
scf: input slide autocc memory
aquarius: cc input scf slide autocc time memory

LOWER_NO_UNDERSCORE = 1
LOWER_UNDERSCORE = 2
UPPER_NO_UNDERSCORE = 3
UPPER_UNDERSCORE = 4

F77COMPILE = $(F77) $(_INCLUDES) $(_F77FLAGS)
F90COMPILE = $(F90) $(_INCLUDES) $(_F90FLAGS)
CCOMPILE = $(CC) $(_DEFS) $(_INCLUDES) $(_CPPFLAGS) $(_CFLAGS)
CXXCOMPILE = $(CXX) $(_DEFS) $(_INCLUDES) $(_CPPFLAGS) $(_CXXFLAGS)
CCOMPILEDEPS = $(CCOMPILE) $(DEPFLAGS)
CXXCOMPILEDEPS = $(CXXCOMPILE) $(DEPFLAGS)

LINK = $(CXX) $(_CXXFLAGS) $(_LDFLAGS) -o $@
ARCHIVE = $(AR) $@

bindir = $(topdir)/bin
libdir = $(topdir)/lib

DEPDIR = .deps
DEPS += $(topdir)/.dummy $(addprefix $(DEPDIR)/,$(notdir $(patsubst %.o,%.Po,$(wildcard *.o))))
ALL_SUBDIRS = $(sort $(SUBDIRS) $(foreach comp,$(ALL_COMPONENTS),$(value $(addsuffix _SUBDIRS,$(comp)))))

_CPPFLAGS = $(CPPFLAGS)
_DEFS = $(DEFS) -DFORTRAN_INTEGER_SIZE=$(FORTRAN_INTEGER_SIZE) -DF77_NAME=$(F77_NAME) -DF90_NAME=$(F90_NAME) -DTOPDIR=\"$(topdir)\"
_LDFLAGS = $(LDFLAGS) -L$(topdir)/lib
_INCLUDES = $(INCLUDES) -I. -I$(topdir) -I$(topdir)/src -I$(CYCLOPSTF)/include -I$(ELEMENTAL)/include
_CFLAGS = $(CFLAGS)
_CXXFLAGS = $(CXXFLAGS)
_F77FLAGS = $(F77FLAGS)
_F90FLAGS = $(F90FLAGS)
_DEPENDENCIES = $(DEPENDENCIES) Makefile $(topdir)/config.mk $(topdir)/rules.mk
_LIBS = $(LIBS) $(CYCLOPSTF_LIBS) $(ELEMENTAL_LIBS) $(BLAS_LIBS)

.NOTPARALLEL:
.PHONY: all default clean $(ALL_COMPONENTS)
FORCE:

$(ALL_COMPONENTS):
	@for dir in $(SUBDIRS) $($@_SUBDIRS); do \
		echo "Making $@ in $$dir"; \
		(cd $$dir && $(MAKE) $@); \
	done

clean:
	rm -rf $(DEPDIR) *.o
	@for subdir in $(ALL_SUBDIRS); do \
		echo "Making clean in $$subdir"; \
		(cd $$subdir && $(MAKE) clean); \
	done

$(bindir)/%: $(_DEPENDENCIES) 
	@mkdir -p $(dir $@)
	$(LINK) $(filter %.o,$^) $(FLIBS) $(_LIBS)

$(libdir)/%: $(_DEPENDENCIES)
	@mkdir -p $(dir $@)
	$(ARCHIVE) $(filter %.o,$^)

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
