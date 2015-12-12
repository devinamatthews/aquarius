bin_PROGRAMS = $(top_builddir)/bin/aquarius
__top_builddir__bin_aquarius_SOURCES = \
	src/autocc/autocc.cxx \
	src/autocc/diagram.cxx \
	src/autocc/fraction.cxx \
	src/autocc/fragment.cxx \
	src/autocc/generator.cxx \
	src/autocc/line.cxx \
	src/autocc/operator.cxx \
	src/autocc/term.cxx \
	\
	src/cc/1edensity.cxx \
	src/cc/2edensity.cxx \
	src/cc/ccd.cxx \
	src/cc/ccsd_density.cxx \
	src/cc/ccsdt_density.cxx \
	src/cc/ccsdtq_density.cxx \
	src/cc/ccsdtq_1a_density.cxx \
	src/cc/ccsdtq_1b_density.cxx \
	src/cc/ccsdtq_3_density.cxx \
	src/cc/cc4_density.cxx \
	src/cc/ccsd_t.cxx \
	src/cc/ccsd_t_l.cxx \
	src/cc/ccsd_t_n_opt.cxx \
	src/cc/ccsd_t_n.cxx \
	src/cc/e_ccsd_t_n.cxx \
	src/cc/ccsd_tq_n_opt.cxx \
	src/cc/ccsd_tq_n.cxx \
	src/cc/ccsd.cxx \
	src/cc/ccsdipgf.cxx \
	src/cc/ccsdtipgf.cxx \
	src/cc/ccsdt_q.cxx \
	src/cc/ccsdt_q_l.cxx \
	src/cc/ccsdt_q_n_opt.cxx \
	src/cc/ccsdt_q_n.cxx \
	src/cc/ccsdt.cxx \
	src/cc/ccsdtq.cxx \
	src/cc/ccsdtq_1a.cxx \
	src/cc/ccsdtq_1b.cxx \
	src/cc/ccsdtq_3.cxx \
	src/cc/cc4.cxx \
	src/cc/eomeeccsd.cxx \
	src/cc/eomeeccsdt.cxx \
	src/cc/lambdaccsd.cxx \
	src/cc/lambdaccsdt.cxx \
	src/cc/lambdaccsdt_q.cxx \
	src/cc/lambdaccsdtq.cxx \
	src/cc/lambdaccsdtq_1a.cxx \
	src/cc/lambdaccsdtq_1b.cxx \
	src/cc/lambdaccsdtq_3.cxx \
	src/cc/lambdacc4.cxx \
	src/cc/lccd.cxx \
	src/cc/mp3.cxx \
	src/cc/mp4dq.cxx \
	src/cc/perturbedccsd.cxx \
	src/cc/perturbedlambdaccsd.cxx \
	src/cc/rhfccsd.cxx \
	src/cc/tda_local.cxx \
	src/cc/rhftda_local.cxx \
	src/cc/rhfeomeeccsd.cxx \
	\
	src/input/basis.cxx \
	src/input/config.cxx \
	src/input/molecule.cxx \
	\
	src/integrals/1eints.cxx \
	src/integrals/2eints.cxx \
	src/integrals/cfour1eints.cxx \
	src/integrals/cfour2eints.cxx \
	src/integrals/center.cxx \
	src/integrals/context.cxx \
	src/integrals/element.cxx \
	src/integrals/fmgamma.cxx \
	src/integrals/kei.cxx \
	src/integrals/nai.cxx \
	src/integrals/os.cxx \
	src/integrals/ovi.cxx \
	src/integrals/shell.cxx \
	\
	src/jellium/jellium.cxx \
	\
	src/main/main.cxx \
	\
	src/operator/2eoperator.cxx \
	src/operator/aomoints.cxx \
	src/operator/rhfaomoints.cxx \
	src/operator/moints.cxx \
	src/operator/sparseaomoints.cxx \
	src/operator/sparserhfaomoints.cxx \
	src/operator/fcidump.cxx \
	\
	src/scf/aouhf.cxx \
	src/scf/cfourscf.cxx \
	src/scf/uhf_local.cxx \
	src/scf/uhf.cxx \
	\
	src/symmetry/symmetry.cxx \
	\
	src/task/task.cxx \
	\
	src/tensor/ctf_tensor.cxx \
	src/tensor/spinorbital_tensor.cxx \
	src/tensor/symblocked_tensor.cxx \
	\
	src/time/time.cxx
	
VPATH += $(srcdir)

AM_CPPFLAGS = -I$(srcdir)/src @ctf_INCLUDES@ @marray_INCLUDES@ @mpiwrap_INCLUDES@ @lawrap_INCLUDES@
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
$(PROGRAMS): src/external/ctf/lib/libctf.a

.PHONY: ALWAYS
src/external/ctf/lib/libctf.a: ALWAYS
	$(MAKE) -C src/external/ctf
endif
