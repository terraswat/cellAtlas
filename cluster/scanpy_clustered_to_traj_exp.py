"""
Makes a tab delimited expression file [cell, gene] per cluster
of scanpy object.

Does gene mean dispersion filtering and then magic on each cluster of a
clustered scanpy object.

example usage:
    scanpy_clustered_to_traj_exp.py -i path/to/scanpy_obj -o path/to/directory/to/fill
"""
import argparse
import sys
import scanpy.api as sc
import pandas as pd
import os

def parse_args():

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i',"--in_file", type=str, help=""
                        )

    parser.add_argument('-o',"--out_dir", type=str, help="",
                        )

    opts = parser.parse_args()
    in_file, out_dir = opts.in_file, opts.out_dir

    return in_file, out_dir


def read_data(in_file):
    return sc.read(in_file)

def data_transform_and_output(out_dir, adata):
    
    MIN_CLUSTER_SIZE = 33
    
    counts = adata.obs.louvain.value_counts()
    print("size of the clusters, less than %d will be filtered out"%MIN_CLUSTER_SIZE)
    print(counts)
    clusters = counts.index[(counts > MIN_CLUSTER_SIZE).tolist()]
    #print(clusters)
    
    for cluster in clusters:
         samples = adata.obs.index[adata.obs.louvain == cluster].tolist()
         subset = adata[samples, :]
          
         subset_filtered = sc.pp.filter_genes_dispersion(subset.X, n_top_genes=3500)
         subset = subset[ :, subset_filtered.gene_subset]
         print("number of genes is %d"%subset.shape[1])
         print("computing magic for cluster_%s"%cluster)
         sc.pp.magic(subset, n_pca=min(subset.X.shape), name_list='all_genes')
         df = pd.DataFrame(
             data = subset.X,
             columns = subset.var.index.tolist(),
             index = subset.obs.index.tolist()
         )  
         cluster_id = "cluster_%s"%cluster
         out_path = os.path.join(out_dir, cluster_id)
         df.to_csv(out_path, sep="\t")


def main():

    in_file, out_dir = parse_args()

    data = read_data(in_file)

    data_transform_and_output(out_dir, data)


if __name__ == "__main__":
    sys.exit(main())
