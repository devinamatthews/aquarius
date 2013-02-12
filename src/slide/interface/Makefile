include ../../../config.mk
include ../../../rules.mk

slide: $(libdir)/libslide.a ../slide.mod
$(libdir)/libslide.a: slide_cxx.o slide_f90_mod.o slide_f90.o

../slide.mod: slide_f90_mod.o
	$(CP) slide.mod ..

INCLUDES += -I$(topdir)/src/slide -I$(topdir)/src/slide/include
