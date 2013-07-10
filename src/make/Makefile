include make.in
include ../../config.mk
SUBDIRS = ${top_dir}/src 
include rules.mk

clean: clean_bin clean_lib

clean_bin:
	rm -rf $(bindir)/*

clean_lib:
	rm -rf $(libdir)/*
