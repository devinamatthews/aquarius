include ../../config.mk
include ../../rules.mk

libs: $(libdir)/libintegrals.a
$(libdir)/libintegrals.a: 1eints.o 2eints.o blasx.o center.o cholesky.o \
                          context.o element.o elementdata.o fmgamma.o \
                          hrr.o ishprim.o keiprim.o momprim.o \
                          naiprim.o osinv.o osprim.o oviprim.o rys.o shell.o vrr.o \
                          cfour1eints.o cfour2eints.o
