
import argparse
import scanpy.api as sc
from scanpyLibrary import *


def commandline_args():
    """
    Command line arguments
    """
    parser = argparse.ArgumentParser(description = "This program takes a scanpy object to preprocess, compute principle component analysis, and cluster.")

    parser.add_argument('-scanpy_obj', type = str, nargs = 1, help = 'The path to the parent directory containing a sample\'s raw scanpy object.')

    parser.add_argument('-output_log_transformed_data', type = str, nargs = 1, help = 'The path to export a sample\'s raw log transformed scanpy object.')

    parser.add_argument('-output_cluster_data', type = str, nargs = 1, help = 'The path to export a sample\'s clustered scanpy object.')

    return parser.parse_args()

def create_sample_dictionary(sample_path):
    """
    Creates a sample dictionary to store the scanpy anndata object

    :param sample_path: A string that is the path to the sample's scanpy h5ad object
    return:
    sample_dict: A dictionary where the key is the sample name and the value is the scanpy anndata object
    """
    sample_dict = {sample_path.split('/')[-1].split('.')[0] : sc.read(sample_path)}
    return sample_dict

def main():

    args = commandline_args()
    sample_dict = create_sample_dictionary(args.scanpy_obj[0])

    # Discuss if this should be a command line argument
    gene_id_conversion_file_path = '/projects/sysbio/users/cellAtlas/data/primary/hugo_to_ensembl.tsv'

    print("...mapping genes and removing duplicate genes...")
    gene_mapped_sample_dict = Gene_Mapping(sample_dict, gene_id_conversion_file_path).sample_gene_mapping()

    print("...filter cells and genes...")
    filter_cell_and_gene_sample_dict = Filter_Cells_And_Genes(gene_mapped_sample_dict).sample_filtering()

    print("...preprocessing data...")
    preprocessed_sample_dict = Basic_Preprocessing(filter_cell_and_gene_sample_dict, args.output_log_transformed_data[0]).sample_preprocessing()

    print("...dimensionality reduction...")
    dimension_reduced_sample_dict = Dimensionality_Reduction(preprocessed_sample_dict).sample_dimension_reduction()
    
    print("...clustering...")
    Cluster(dimension_reduced_sample_dict, args.output_cluster_data[0]).sample_clustering()
    
    print("---Done---")

if __name__ == '__main__':
    main()
