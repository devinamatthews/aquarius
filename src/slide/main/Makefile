include ../../../config.mk
include ../../make/make.in
include ../../make/rules.mk

libs: $(libdir)/libslide.a
$(libdir)/libslide.a: center.o context.o element.o integrals.o shell.o slide.o

INCLUDES += -I$(topdir)/src/slide -I$(topdir)/src/slide/include
