include ../../config.mk
include ../../rules.mk

bins: $(bindir)/aquarius
$(bindir)/aquarius: main.o $(libdir)/libcc.a $(libdir)/libscf.a \
                    $(libdir)/libautocc.a $(libdir)/libinput.a \
                    $(libdir)/libtime.a $(libdir)/libslide.a \
                    $(libdir)/libmemory.a $(libdir)/libtensor.a \
                    $(libdir)/libutil.a

LIBS += -lcc -lscf -lautocc -linput -ltime -lslide -ltensor -lmemory -lutil
