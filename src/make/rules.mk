all: $(DEFAULT_COMPONENTS)

ALL_COMPONENTS = libs bins

bins: libs

LOWER_NO_UNDERSCORE = 1
LOWER_UNDERSCORE = 2
UPPER_NO_UNDERSCORE = 3
UPPER_UNDERSCORE = 4


bindir = ${top_dir}/bin
libdir = ${top_dir}/lib

DEPDIR = .deps
DEPS += ${top_dir}/.dummy $(addprefix $(DEPDIR)/,$(notdir $(patsubst %.o,%.Po,$(wildcard *.o))))
ALL_SUBDIRS = $(sort $(SUBDIRS) $(foreach comp,$(ALL_COMPONENTS),$(value $(addsuffix _SUBDIRS,$(comp)))))

_CPPFLAGS = $(CPPFLAGS)
_DEFS = $(DEFS) -DFORTRAN_INTEGER_SIZE=$(FORTRAN_INTEGER_SIZE) -DF77_NAME=$(F77_NAME) -DF90_NAME=$(F90_NAME) -DTOPDIR=\"${top_dir}\"
_LDFLAGS = $(LDFLAGS) -L${top_dir}/lib
_INCLUDES = $(INCLUDES) -I. -I${top_dir} -I${top_dir}/src -I$(CTF_DIR)/include -I$(ELEMENTAL)/include
_CFLAGS = $(CFLAGS)
_CXXFLAGS = $(CXXFLAGS)
_F77FLAGS = $(F77FLAGS)
_F90FLAGS = $(F90FLAGS)
_DEPENDENCIES = $(DEPENDENCIES) Makefile ${top_dir}/config.mk ${top_dir}/rules.mk
_LIBS = $(LIBS) $(CTF_LIBS) $(ELEMENTAL_LIBS) $(BLAS_LIBS)

F77COMPILE = $(F77) $(_INCLUDES) $(_F77FLAGS)
F90COMPILE = $(F90) $(_INCLUDES) $(_F90FLAGS)
CCOMPILE = $(CC) $(_DEFS) $(_INCLUDES) $(_CPPFLAGS) $(_CFLAGS)
CXXCOMPILE = $(CXX) $(_DEFS) $(_INCLUDES) $(_CPPFLAGS) $(_CXXFLAGS)
CCOMPILEDEPS = $(CCOMPILE) $(DEPFLAGS)
CXXCOMPILEDEPS = $(CXXCOMPILE) $(DEPFLAGS)

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
.PHONY: all default clean $(ALL_COMPONENTS)
FORCE:

$(ALL_COMPONENTS):
	@for dir in $(SUBDIRS) $($@_SUBDIRS); do \
    echo "top dir is"; \
    echo ${top_dir}; \
		echo "Making $@ in $$dir"; \
		(cd $$dir && $(MAKE) $@); \
	done

clean:
	rm -rf $(DEPDIR) *.o
	@for subdir in $(ALL_SUBDIRS); do \
		echo "Making clean in $$subdir"; \
		(cd $$subdir && $(MAKE) clean); \
	done

$(bindir)/%: $(_DEPENDENCIES) $(call libdeps,$(_LIBS))
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
