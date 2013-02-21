include ../../../config.mk
include ../../../rules.mk

libs: $(libdir)/libslide.a
$(libdir)/libslide.a: center.o context.o element.o integrals.o shell.o slide.o

INCLUDES += -I$(topdir)/src/slide -I$(topdir)/src/slide/include
