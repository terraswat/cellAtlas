#!/bin/bash

# set a list of scanpy object file names (can contain only one filename)
basedirectory="/projects/sysbio/users/cellAtlas/data/scanpyObj/"
objects=("10xGenomics_pbmc_8k_4k")
#objects=("10xGenomics_t_4k" "10xGenomics_t_3k")
#objects=("kriegstein6k")
#objects=("10xGenomics_pbmc8k" "10xGenomics_t_3k_4k_aggregate")
#objects=("GSE94820_DCnMono.discovery_PBMCs_Broad" "GSE94820_deeper.characterization_PBMCs_Broad")
#objects=("mouse_tabulaMuris_droplet_Marrow" "mouse_tabulaMuris_facs_Marrow" "mouse_tabulaMuris_droplet_Limb_Muscle" "mouse_tabulaMuris_facs_Limb_Muscle" "mouse_tabulaMuris_droplet_Heart_and_Aorta" "mouse_tabulaMuris_facs_Heart" "mouse_tabulaMuris_droplet_Tongue" "mouse_tabulaMuris_facs_Tongue" "mouse_tabulaMuris_droplet_Bladder" "mouse_tabulaMuris_droplet_Spleen" "mouse_tabulaMuris_droplet_Thymus" "mouse_tabulaMuris_droplet_Kidney" "mouse_tabulaMuris_droplet_Liver" "mouse_tabulaMuris_droplet_Lung" "mouse_tabulaMuris_droplet_Mammary_Gland" "mouse_tabulaMuris_droplet_Trachea" "mouse_tabulaMuris_facs_Bladder" "mouse_tabulaMuris_facs_Brain_Myeloid" "mouse_tabulaMuris_facs_Brain_Non-Myeloid" "mouse_tabulaMuris_facs_Fat" "mouse_tabulaMuris_facs_Kidney" "mouse_tabulaMuris_facs_Large_Intestine" "mouse_tabulaMuris_facs_Liver" "mouse_tabulaMuris_facs_Lung" "mouse_tabulaMuris_facs_Mammary_Gland" "mouse_tabulaMuris_facs_Pancreas" "mouse_tabulaMuris_facs_Skin" "mouse_tabulaMuris_facs_Spleen" "mouse_tabulaMuris_facs_Thymus" "mouse_tabulaMuris_facs_Trachea")
#objects=("GSE52529-GPL11154_myoblasts_Trapnell")
#objects=("SRP073808_embryonicStemCells_Koh" "GSE64016_H1andFUCCI_Leng" "GSE63818_primordial_germ_cells_Guo" "GSE52529-GPL16791_myoblasts_Trapnell" "GSE52529-GPL11154_myoblasts_Trapnell" "GSE79102_myeloproloferative_disease_Kiselev")
#objects=("GSE84465_GBM_DarmanisQuake" "EMTAB3929_preimplantation_embryos" "GSE81076-GPL18573_pancreas_Grun" "GSE81076-GPL16791_pancreas_Grun" "GSE57872_GBM_PatelBernstein" "SRP073808_embryonicStemCells_Koh" "GSE64016_H1andFUCCI_Leng" "GSE63818_primordial_germ_cells_Guo" "GSE52529-GPL16791_myoblasts_Trapnell" "GSE52529-GPL11154_myoblasts_Trapnell" "GSE79102_myeloproloferative_disease_Kiselev")
#objects=("ica_bone_marrow.h5ad")
#objects=("quakeBrainGeo1.h5ad")
#objects=("GSE63818_primordial_germ_cells_Guo.h5ad")
#listOfFiles=( "${basedirectory}${objects[@]/%/}" )

# use all the files in a directory containing scanpy objects
#listOfFiles=(/projects/sysbio/users/cellAtlas/data/scanpyObj/*)



# loop through the list and execute for each scanpy object:
for obj in "${objects[@]}"
do
    INFILE="${basedirectory}${obj}.h5ad"
    LOGFILE="/projects/sysbio/users/cellAtlas/data/scanpyObj/logarithmized_raw_data/${obj}_raw_log.h5ad"
    CLUSTERFILE="/projects/sysbio/users/cellAtlas/data/cluster/${obj}_clustered.h5ad"
    
    echo "------------------------------------------"
    echo "start working on data set: $obj"
    #echo "infile $INFILE"
    #echo "LOGFILE $LOGFILE"
    #echo "clusterfile $CLUSTERFILE"
    
    # run the scanpyToCluster.py script
    python ./scanpyToCluster.py -scanpy_obj $INFILE -output_log_transformed_data $LOGFILE -output_cluster_data $CLUSTERFILE
    #echo "python ./scanpyToCluster.py -scanpy_obj $INFILE -output_log_transformed_data $LOGFILE -output_cluster_data $CLUSTERFILE"
done



