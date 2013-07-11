# Aquarius
#
# see README for more detailed instructions
#
# type make to build Aquarius
# type make examples to build all examples
# type make test
#   Devin please write something that does a binary CC functionality test -Edgar

HNAME	:= $(shell hostname | cut -d . -f 1)

.PHONY: examples
all $(MAKECMDGOALS):
	@if [ ! -f config.mk ] ;  then \
	  echo -n "topdir = " > config.mk; \
	  pwd >> config.mk; \
	  echo >> config.mk; \
	  if [ $(shell hostname | grep 'edison\|hopper\|carver' ) ] ;  then \
	    echo 'Hostname recognized as a NERSC machine, using pre-made config.mk file'; \
	    cat mkfiles/config.mk.nersc >> config.mk;   \
	  elif [ $(shell hostname | grep 'surveyor\|intrepid\|challenger\|udawn' ) ] ;  then \
	    echo 'Hostname recognized as a BG/P machine, using pre-made config.mk file'; \
	    cat mkfiles/config.mk.bgp >> config.mk;   \
	  elif [ $(shell hostname | grep 'vesta\|mira\|cetus\|seq' ) ] ;  then \
	    cat mkfiles/config.mk.bgq >> config.mk;   \
	  else \
	    echo 'Hostname not recognized: assuming linux, specialize config.mk if necessary'; \
	    cat mkfiles/config.mk.linux >> config.mk;   \
	  fi; \
	fi; \
	cd src/make; \
	$(MAKE) $@;
