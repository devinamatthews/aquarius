include ../../config.mk
include ../../rules.mk

aquarius: $(bindir)/aquarius
$(bindir)/aquarius: main.o $(libdir)/libcc.a $(libdir)/libscf.a \
                    $(libdir)/libautocc.a $(libdir)/libinout.a \
                    $(libdir)/libtime.a $(libdir)/libslide.a \
                    $(libdir)/libmemory.a

LIBS += -lcc -lscf -lautocc -linput -ltime -lslide -lmemory
