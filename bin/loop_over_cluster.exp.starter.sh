#!/bin/bash
ATTR_FILTERED_SCR=/projects/sysbio/users/cellAtlas/bin/attributes/filter.py
CLUSTERDIR=/projects/sysbio/users/cellAtlas/data/cluster.exp
OUTPUTDIR=/projects/sysbio/users/cellAtlas/data/cluster.exp.less_1000

for dir in ${CLUSTERDIR}/*; do
	echo $dir
	for filename in ${dir}/*; do
		echo $(basename $filename)
		echo $(wc -l $filename)
	done
	#OUT=${OUTPUTDIR}/$(basename $filename)
	#OUT=${OUT%.*}.attr.filtered.tab

        #KEEP=${KEEPDIR}/$(basename $filename)
        #KEEP=${KEEP%.*}.attr.keep
        #echo $KEEP
        #echo $OUT
	#python3 $ATTR_FILTERED_SCR -i $filename -o $OUT -k $KEEP
done
