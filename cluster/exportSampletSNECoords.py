import argparse
import scanpy.api as sc
from scanpyLibrary import *


def commandline_args():
    """
    Command line arguments
    """
    parser = argparse.ArgumentParser(description = "This program takes a clustered scanpy object and exports each clusters expression matrix in tsv format. Rows are cells and genes are columns.")

    parser.add_argument('-clustered_scanpy_obj', type = str, nargs = 1, help = 'The path to a sample\'s clustered scanpy object.')

    parser.add_argument('-output_sample_tSNE_coords', type = str, nargs = 1, help = 'The path to export a sample\'s log transformed expression matrix.')

    return parser.parse_args()

def create_sample_dictionary(sample_path):
    """
    Creates a sample dictionary to store the scanpy anndata object

    :param sample_path: A string that is the path to the sample's scanpy h5ad object
    return:
    sample_dict: A dictionary where the key is the sample name and the value is the scanpy anndata object
    """
    sample_dict = {'_'.join(sample_path.split('/')[-1].split('.')[0].split('_')[:-1]) : sc.read(sample_path)}
    return sample_dict

def generate_sample_tSNE_coords_output_file(sample_cluster_dict, output_dir):
    """
    :param sample_clustered_dict:
    :param output_dir:
    return:
    """
    for sample_key, adata in sample_cluster_dict.items():
        tsne_coord = pd.DataFrame(adata.obsm.X_tsne, index=adata.obs.index.tolist())
        tsne_coord.columns = ['x', 'y']
        tsne_coord.to_csv(output_dir + sample_key + "_tSNE_coords.tsv", sep='\t')

def main():

    args = commandline_args()
    sample_cluster_dict = create_sample_dictionary(args.clustered_scanpy_obj[0])
    generate_sample_tSNE_coords_output_file(sample_cluster_dict, args.output_sample_tSNE_coords[0])

if __name__ == '__main__':
    main()
