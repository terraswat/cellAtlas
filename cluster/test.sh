python exportClusterRankings.py -clustered_scanpy_obj /projects/sysbio/users/cellAtlas/data/cluster/quakeBrainGeo1_clustered.h5ad -output_cluster_rankings /projects/sysbio/users/cellAtlas/data/cluster.rnks/

python exportClusterExpressionMatrix.py -clustered_scanpy_obj /projects/sysbio/users/cellAtlas/data/cluster/quakeBrainGeo1_clustered.h5ad -raw_log_scanpy_obj /projects/sysbio/users/cellAtlas/data/scanpyObj/logarithmized_raw_data/quakeBrainGeo1_raw_log.h5ad -output_cluster_expression /projects/sysbio/users/cellAtlas/data/cluster.exp/

python exportSampleExpressionMatrix.py -clustered_scanpy_obj /projects/sysbio/users/cellAtlas/data/cluster/quakeBrainGeo1_clustered.h5ad -raw_log_scanpy_obj /projects/sysbio/users/cellAtlas/data/scanpyObj/logarithmized_raw_data/quakeBrainGeo1_raw_log.h5ad -output_sample_expression /projects/sysbio/users/cellAtlas/data/sample.exp/

python exportSampletSNECoords.py -clustered_scanpy_obj /projects/sysbio/users/cellAtlas/data/cluster/quakeBrainGeo1_clustered.h5ad -output_sample_tSNE_coords /projects/sysbio/users/cellAtlas/data/tsne.coords/

python exportSampleUMAPCoords.py -clustered_scanpy_obj /projects/sysbio/users/cellAtlas/data/cluster/quakeBrainGeo1_clustered.h5ad -output_sample_umap_coords /projects/sysbio/users/cellAtlas/data/umap.coords/