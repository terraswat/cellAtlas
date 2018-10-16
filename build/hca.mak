# chrisw

THIS_DIR=$(shell pwd)

DATA_DIR=data
LIB_DIR=lib

MARROW_H5=$(DATA_DIR)/immune_census/ica_bone_marrow_h5.h5
CORD_BLOOD_H5=$(DATA_DIR)/immune_census/ica_bone_marrow_h5.h5

TARGETS=

test:marrow.tsv

marrow.tsv:
	h5dump --noindex --width=0 --output=1.tmp $(MARROW_H5) ;
	\

all: $(TARGETS)

clean_all: clean_targets clean_tmp

clean_targets:
	rm -f $(TARGETS) ;

clean_tmp:
	rm -f $(wildcard *.tmp) ;

