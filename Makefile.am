bin_PROGRAMS = $(top_builddir)/bin/aquarius
__top_builddir__bin_aquarius_SOURCES = \
	src/agora/integrals/hgp/hgp.cxx \
	src/agora/integrals/ishida/ishida.cxx \
	src/agora/integrals/ishida/kei.cxx \
	src/agora/integrals/ishida/nai.cxx \
	src/agora/integrals/ishida/ovi.cxx \
	src/agora/integrals/obara_saika/moments.cxx \
	src/agora/integrals/obara_saika/os.cxx \
	\
	src/frameworks/autocc/diagram.cxx \
	src/frameworks/autocc/fraction.cxx \
	src/frameworks/autocc/fragment.cxx \
	src/frameworks/autocc/generator.cxx \
	src/frameworks/autocc/line.cxx \
	src/frameworks/autocc/operator.cxx \
	src/frameworks/autocc/term.cxx \
	\
	src/frameworks/integrals/1eints.cxx \
	src/frameworks/integrals/2eints.cxx \
	src/frameworks/integrals/cholesky.cxx \
	src/frameworks/integrals/context.cxx \
	src/frameworks/integrals/fmgamma.cxx \
	src/frameworks/integrals/momints.cxx \
	src/frameworks/integrals/nucints.cxx \
	src/frameworks/integrals/rys.cxx \
	\
	src/frameworks/logging/logging.cxx \
	\
	src/frameworks/molecule/basis.cxx \
	src/frameworks/molecule/center.cxx \
	src/frameworks/molecule/element.cxx \
	src/frameworks/molecule/molecule.cxx \
	src/frameworks/molecule/shell.cxx \
	\
	src/frameworks/symmetry/symmetry.cxx \
	\
	src/frameworks/task/config.cxx \
	src/frameworks/task/task.cxx \
	\
	src/frameworks/tensor/implementations/ctf_tensor.cxx \
	\
	src/frameworks/time/time.cxx \
	\
	src/frameworks/util/distributed.cxx \
	\
	src/main/main.cxx \
	\
	src/tasks/jellium/jellium.cxx
	
VPATH += $(srcdir)

AM_CPPFLAGS = -I$(srcdir)/src -I$(srcdir)/external/lawrap \
	-I$(srcdir)/external/marray/include -I$(srcdir)/external/mpiwrap/include \
	-I$(srcdir)/external/stl_ext/include @ctf_INCLUDES@
__top_builddir__bin_aquarius_LDADD = @ctf_LIBS@
	
if HAVE_ELEMENTAL
AM_CPPFLAGS += @elemental_INCLUDES@
__top_builddir__bin_aquarius_LDADD += @elemental_LIBS@
__top_builddir__bin_aquarius_SOURCES += \
	src/scf/uhf_elemental.cxx \
	src/cc/tda_elemental.cxx \
	src/cc/rhftda_elemental.cxx
endif
	
if HAVE_LIBINT2
AM_CPPFLAGS += @libint2_INCLUDES@
__top_builddir__bin_aquarius_LDADD += @libint2_LIBS@
__top_builddir__bin_aquarius_SOURCES += src/integrals/libint2eints.cxx
endif

__top_builddir__bin_aquarius_LDADD += $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

if CTF_IS_LOCAL
$(PROGRAMS): external/ctf/lib/libctf.a

.PHONY: ALWAYS
external/ctf/lib/libctf.a: ALWAYS
	$(MAKE) -C external/ctf
endif
