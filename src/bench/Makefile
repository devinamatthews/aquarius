include ../../config.mk
include ../make/make.in
include ../make/rules.mk

bins: $(bindir)/bench-cholesky-ccsd-lambda
$(bindir)/bench-cholesky-ccsd-lambda: cholesky-ccsd-lambda.o \
                                      $(libdir)/libautocc.a $(libdir)/libinput.a \
                                      $(libdir)/libtime.a $(libdir)/libslide.a \
                                      $(libdir)/libmemory.a $(libdir)/libtensor.a \
                                      $(libdir)/libutil.a

bins: $(bindir)/bench-cholesky-ccsd
$(bindir)/bench-cholesky-ccsd: cholesky-ccsd.o \
							   $(libdir)/libautocc.a $(libdir)/libinput.a \
                               $(libdir)/libtime.a $(libdir)/libslide.a \
                               $(libdir)/libmemory.a $(libdir)/libtensor.a \
                               $(libdir)/libutil.a

bins: $(bindir)/bench-cholesky-ccsdt
$(bindir)/bench-cholesky-ccsdt: cholesky-ccsdt.o \
                                $(libdir)/libautocc.a $(libdir)/libinput.a \
                                $(libdir)/libtime.a $(libdir)/libslide.a \
                                $(libdir)/libmemory.a $(libdir)/libtensor.a \
                                $(libdir)/libutil.a

bins: $(bindir)/bench-ao-ccsd-lambda
$(bindir)/bench-ao-ccsd-lambda: ao-ccsd-lambda.o \
                                $(libdir)/libautocc.a $(libdir)/libinput.a \
                                $(libdir)/libtime.a $(libdir)/libslide.a \
                                $(libdir)/libmemory.a $(libdir)/libtensor.a \
                                $(libdir)/libutil.a
  
bins: $(bindir)/bench-ao-ccsd
$(bindir)/bench-ao-ccsd: ao-ccsd.o \
						 $(libdir)/libautocc.a $(libdir)/libinput.a \
                         $(libdir)/libtime.a $(libdir)/libslide.a \
                         $(libdir)/libmemory.a $(libdir)/libtensor.a \
                         $(libdir)/libutil.a

bins: $(bindir)/bench-ao-ccsdt
$(bindir)/bench-ao-ccsdt: ao-ccsdt.o \
                          $(libdir)/libautocc.a $(libdir)/libinput.a \
                          $(libdir)/libtime.a $(libdir)/libslide.a \
                          $(libdir)/libmemory.a $(libdir)/libtensor.a \
                          $(libdir)/libutil.a

LIBS += -lautocc -linput -ltime -lslide -ltensor -lmemory -lutil
