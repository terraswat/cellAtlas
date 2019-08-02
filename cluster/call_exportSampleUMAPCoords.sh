#!/bin/bash

# set a list of scanpy object names (without file extension)
objects=("GSE84465_GBM_DarmanisQuake" "EMTAB3929_preimplantation_embryos" "GSE81076-GPL18573_pancreas_Grun" "GSE81076-GPL16791_pancreas_Grun" "GSE57872_GBM_PatelBernstein" "SRP073808_embryonicStemCells_Koh" "GSE64016_H1andFUCCI_Leng" "GSE63818_primordial_germ_cells_Guo" "GSE52529-GPL16791_myoblasts_Trapnell" "GSE79102_myeloproloferative_disease_Kiselev")
#objects=("quakeBrainGeo1" "GSE84465_GBM_DarmanisQuake" "GSE63818_primordial_germ_cells_Guo")
#objects=("ica_bone_marrow")
#objects=("ica_cord_blood")
#objects=("mouse_tabulaMuris_droplet" "mouse_tabulaMuris_facs")
#objects=("GSE52529-GPL11154_myoblasts_Trapnell")

# the directory of the clustered scanpy object
clusterdir="/projects/sysbio/users/cellAtlas/data/cluster/"
#clusterListOfFiles=( "${clusterdirectory}${objects[@]/%/}.h5ad" )

# the directory of the scanpy object containing the normalized full expression matrix
expressiondir="/projects/sysbio/users/cellAtlas/data/scanpyObj/logarithmized_raw_data/"
#ExpressionListOfFiles=( "${expressiondirectory}${objects[@]/%/}.h5ad" )

# output directory for rank files per cluster
outdir="/projects/sysbio/users/cellAtlas/data/umap.coords/"




# loop through the list and execute for each scanpy object:
for obj in "${objects[@]}"
do
    CLUSTERFILE="${clusterdir}${obj}_clustered.h5ad"
    
    echo "------------------------------------------"
    echo "exporting UMAP coords for $obj"
    #echo "clusterfile $CLUSTERFILE"

    # run the exportSampleUMAPCoords.py script
    python exportSampleUMAPCoords.py -clustered_scanpy_obj $CLUSTERFILE -output_sample_umap_coords $outdir


done



