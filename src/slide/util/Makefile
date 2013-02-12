include ../../../config.mk
include ../../../rules.mk

slide: $(libdir)/libslide.a
$(libdir)/libslide.a: blasx.o memory.o symmetry.o transform.o symdata.o \
                      elementdata.o

INCLUDES += -I$(topdir)/src/slide -I$(topdir)/src/slide/include
