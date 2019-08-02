#!/bin/bash
# Go through all the clustered scanpy objects and
# attribute enrichment analysis.
# If each directory for the output is not made you need to uncomment the
# mkdir command in the for loop below.
ATTR_FILTERED_SCR=/projects/sysbio/users/cellAtlas/bin/attributes/filter.py
CLUSTERDIR=/projects/sysbio/users/cellAtlas/data/cluster
OUTPUTDIR=/projects/sysbio/users/cellAtlas/data/attr.filtered
KEEPDIR=/projects/sysbio/users/cellAtlas/data/attr.keep

source /projects/sysbio/users/cellAtlas/env/bin/activate

for filename in ${CLUSTERDIR}/*; do
	echo $(basename $filename)
	OUT=${OUTPUTDIR}/$(basename $filename)
	OUT=${OUT%.*}.attr.filtered.tab

        KEEP=${KEEPDIR}/$(basename $filename)
        KEEP=${KEEP%.*}.attr.keep
        #echo $KEEP
        #echo $OUT
	python3 $ATTR_FILTERED_SCR -i $filename -o $OUT -k $KEEP
done
