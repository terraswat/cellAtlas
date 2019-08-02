#!/bin/bash
# Go through all the clustered scanpy objects and
# attribute enrichment analysis.
# If each directory for the output is not made you need to uncomment the
# mkdir command in the for loop below.
CLUSTER_EXP_SCR=/projects/sysbio/users/cellAtlas/cluster/scanpy_clustered_to_traj_exp.py
CLUSTERDIR=/projects/sysbio/users/cellAtlas/data/cluster
OUTPUTDIR=/projects/sysbio/users/cellAtlas/data/cluster.exp

source /projects/sysbio/users/cellAtlas/env/bin/activate

for filename in ${CLUSTERDIR}/*; do
	echo $(basename $filename)
	OUT=${OUTPUTDIR}/$(basename $filename)
	OUT=${OUT%.*}
	#mkdir $OUT
	python3 $CLUSTER_EXP_SCR -i $filename -o $OUT
done
