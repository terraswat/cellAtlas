#!/bin/bash

# set a list of scanpy object file names (can contain only one filename)
basedirectory="/projects/sysbio/users/cellAtlas/data/scanpyObj/"
objects=("ucsc_human_cortex.h5ad" "quakeBrainGeo1.h5ad")
listOfFiles=( "${basedirectory}${objects[@]/%/}" )

# use all the files in a directory containing scanpy objects
#listOfFiles=(/projects/sysbio/users/cellAtlas/data/scanpyObj/*)



# loop through the list and execute for each scanpy object:
for INFILE in "${listOfFiles[@]}"
do
    SCANPYOBJECT=$(basename "$INFILE" .h5ad)
    LOGFILE="/projects/sysbio/users/cellAtlas/data/scanpyObj/${SCANPYOBJECT}_raw_log.h5ad"
    CLUSTERFILE="/projects/sysbio/users/cellAtlas/data/cluster/${SCANPYOBJECT}_clustered.h5ad"
    
    echo "clustering $SCANPYOBJECT"
    #echo "infile $INFILE"
    #echo "LOGFILE $LOGFILE"
    #echo "clusterfile $CLUSTERFILE"
    
    # run the scanpyToCluster.py script
    python ./scanpyToCluster.py -scanpy_obj $INFILE -output_log_transformed_data $LOGFILE -output_cluster_data $CLUSTERFILE

done



