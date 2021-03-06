#!/bin/bash

# set a list of scanpy object names (without file extension)
objects=("GSE94820_DCnMono_discovery_PBMCs_Broad" "GSE94820_deeper_characterization_PBMCs_Broad")
#objects=("kriegstein6k" "10xGenomics_pbmc4k" "10xGenomics_pbmc8k" "10xGenomics_t_3k_4k_aggregate" "GSE94820_DCnMono_discovery_PBMCs_Broad" "GSE94820_deeper_characterization_PBMCs_Broad")
#objects=("mouse_tabulaMuris_droplet_Marrow" "mouse_tabulaMuris_facs_Marrow" "mouse_tabulaMuris_droplet_Limb_Muscle" "mouse_tabulaMuris_facs_Limb_Muscle" "mouse_tabulaMuris_droplet_Heart_and_Aorta" "mouse_tabulaMuris_facs_Heart" "mouse_tabulaMuris_droplet_Tongue" "mouse_tabulaMuris_facs_Tongue" "mouse_tabulaMuris_droplet_Bladder" "mouse_tabulaMuris_droplet_Spleen" "mouse_tabulaMuris_droplet_Thymus" "mouse_tabulaMuris_droplet_Kidney" "mouse_tabulaMuris_droplet_Liver" "mouse_tabulaMuris_droplet_Lung" "mouse_tabulaMuris_droplet_Mammary_Gland" "mouse_tabulaMuris_droplet_Trachea" "mouse_tabulaMuris_facs_Bladder" "mouse_tabulaMuris_facs_Brain_Myeloid" "mouse_tabulaMuris_facs_Brain_Non-Myeloid" "mouse_tabulaMuris_facs_Fat" "mouse_tabulaMuris_facs_Kidney" "mouse_tabulaMuris_facs_Large_Intestine" "mouse_tabulaMuris_facs_Liver" "mouse_tabulaMuris_facs_Lung" "mouse_tabulaMuris_facs_Mammary_Gland" "mouse_tabulaMuris_facs_Pancreas" "mouse_tabulaMuris_facs_Skin" "mouse_tabulaMuris_facs_Spleen" "mouse_tabulaMuris_facs_Thymus" "mouse_tabulaMuris_facs_Trachea")
#objects=("GSE84465_GBM_DarmanisQuake" "EMTAB3929_preimplantation_embryos" "GSE81076-GPL18573_pancreas_Grun" "GSE81076-GPL16791_pancreas_Grun" "GSE57872_GBM_PatelBernstein" "SRP073808_embryonicStemCells_Koh" "GSE64016_H1andFUCCI_Leng" "GSE63818_primordial_germ_cells_Guo" "GSE52529-GPL16791_myoblasts_Trapnell" "GSE52529-GPL11154_myoblasts_Trapnell" "GSE79102_myeloproloferative_disease_Kiselev")
#objects=("ica_bone_marrow")
#objects=("ica_cord_blood")
#objects=("mouse_tabulaMuris_droplet" "mouse_tabulaMuris_facs")

# the directory of the clustered scanpy object
clusterdir="/projects/sysbio/users/cellAtlas/data/cluster/"
#clusterListOfFiles=( "${clusterdirectory}${objects[@]/%/}.h5ad" )

# output directory for rank files per cluster
rnkdir="/projects/sysbio/users/cellAtlas/data/cluster.rnks/"




# loop through the list and execute for each scanpy object:
for obj in "${objects[@]}"
do
    CLUSTERFILE="${clusterdir}${obj}_clustered.h5ad"
    
    echo "------------------------------------------"
    echo "exporting rnk files for $obj"
    
    
    # make directory to write rank file per cluster for the data set
    OUTPUTDIR="${rnkdir}${obj}/"
    mkdir -p $OUTPUTDIR
    
    #echo "clusterfile $CLUSTERFILE"
    #echo "output dir $OUTPUTDIR"
    
    # run the exportClusterRankings.py script
    python exportClusterRankings.py -clustered_scanpy_obj $CLUSTERFILE -output_cluster_rankings $rnkdir


done



