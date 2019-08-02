#!/bin/bash


# the directory of the clustered scanpy object
clusterdir="/projects/sysbio/users/cellAtlas/data/cluster/"
#clusterListOfFiles=( "${clusterdirectory}${objects[@]/%/}.h5ad" )

# the directory of the scanpy object containing the normalized full expression matrix
expressiondir="/projects/sysbio/users/cellAtlas/data/scanpyObj/logarithmized_raw_data/"
#ExpressionListOfFiles=( "${expressiondirectory}${objects[@]/%/}.h5ad" )

# output directory for rank files per cluster
outdir="/projects/sysbio/users/cellAtlas/data/umap.coords/"




# loop through the list and execute for each scanpy object:
for fn in `ls /projects/sysbio/users/cellAtlas/data/cluster/ |grep .h5ad`;
do
    echo "------------------------------------------"
    echo "exporting UMAP coords for $fn"
    #echo "clusterfile $CLUSTERFILE"

    # run the exportSampleUMAPCoords.py script
    python /projects/sysbio/users/cellAtlas/cluster/exportSampleUMAPCoords.py -clustered_scanpy_obj $fn -output_sample_umap_coords $outdir


done
