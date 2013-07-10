include ../../../config.mk
include ../../make/make.in
include ../../make/rules.mk

libs: $(libdir)/libslide.a
$(libdir)/libslide.a: fmgamma.o hrr.o ish.o kei.o nai.o os.o osinv.o rys.o \
                      vrr.o ovi.o moment.o

INCLUDES += -I$(topdir)/src/slide -I$(topdir)/src/slide/include
