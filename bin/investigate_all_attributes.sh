#!/bin/bash
# Go through all the attributes in the clustered scanpy objects and
# print a summary of them so we can decide which ones should go into the
# attribute enrichment analysis.
ATTRI_SCR=/projects/sysbio/users/cellAtlas/bin/attributes/attribute_investigate.py
CLUSTERDIR=/projects/sysbio/users/cellAtlas/data/cluster
OUTPUTDIR=/projects/sysbio/users/cellAtlas/data/attr.descriptions

source /projects/sysbio/users/cellAtlas/env/bin/activate

for filename in ${CLUSTERDIR}/*; do
	OUT=${OUTPUTDIR}/$(basename $filename)
	OUT=${OUT%.*}.attr.desc
	#echo $OUT
	python3 $ATTRI_SCR -i $filename -o $OUT
done
