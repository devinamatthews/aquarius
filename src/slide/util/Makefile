include ../../../config.mk
include ../../make/make.in
include ../../make/rules.mk

libs: $(libdir)/libslide.a
$(libdir)/libslide.a: blasx.o memory.o symmetry.o transform.o symdata.o \
                      elementdata.o

INCLUDES += -I$(topdir)/src/slide -I$(topdir)/src/slide/include
