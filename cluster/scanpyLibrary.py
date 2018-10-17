import numpy as np
import pandas as pd
import scanpy.api as sc
from sklearn.preprocessing import StandardScaler

class Gene_Mapping(object):

    def __init__(self, sample_dict, gene_id_conversion_file):

        self.sample_dict = sample_dict
        self.gene_id_conversion_file = gene_id_conversion_file

    def generate_gene_mapping_dictionary(self):
        """
        Create an ensembl to Hugo gene conversion dictionary.

        :param gene_id_conversion_file: A string that is the path to the file containing the ensembl to hugo gene mappings.

        return:
        ensembl2symbol: A dictionary where the keys are the ensembl id's and the values are the hugo id's
        """

        symbol2ensembl = pd.read_csv(self.gene_id_conversion_file, sep='\t', index_col=1)['ensembl']

        symbol2ensembl = symbol2ensembl[symbol2ensembl.notnull()]

        ensembl2symbol = dict(zip(symbol2ensembl.values,symbol2ensembl.index.values))

        return ensembl2symbol


    def sample_adata_ensembl_to_hugo_conversion(self, adata, ensembl2symbol):
        """
        Converts a sample's ensembl gene space to hugo in adata.var and create unique gene names
        for multiple Hugo mappings

        :param adata: A sample's scanpy anndata object
        :param ensembl2symbol: A dictionary where the keys are the ensembl id's and the values are the hugo id's

        return: Updated scanpy adata.var to unique Hugo name space 
        """
        adata.var = adata.var.rename(index=ensembl2symbol)
        adata.var_names_make_unique()

    def sample_gene_mapping(self):
        """
        """
        ensembl2symbol = self.generate_gene_mapping_dictionary()
        for sample_key, adata in self.sample_dict.items():
            self.sample_adata_ensembl_to_hugo_conversion(adata, ensembl2symbol)
        return self.sample_dict

class Filter_Cells_And_Genes(object):

    def __init__(self, sample_dict):

        self.sample_dict = sample_dict

    def filter_cells(self, adata):
        """
        Filter cell outliers based on number of genes expressed. A cell
        must have a minimum number of three genes expressed to pass filtering.

        This method adds the columns 'n_genes', which is the number of genes expressed
        per cell.

        References
        https://scanpy.readthedocs.io/en/latest/api/scanpy.api.pp.filter_cells.html?highlight=filter%20cells

        :param adata: A sample's scanpy anndata object

        return: None
        """
        sc.pp.filter_cells(adata, min_genes=3)


    def filter_genes(self, adata):
        """
        Filter genes based on number of cells or counts. A gene must be expressed in a minimum of two hundred
        cells to pass filtering

        This method adds the columns 'n_cells', which is the number of cells that express a gene, to adata.var.

        References
        https://scanpy.readthedocs.io/en/latest/api/scanpy.api.pp.filter_genes.html?highlight=filter%20genes

        :param adata: A sample's scanpy anndata object
        return: None
        """
        sc.pp.filter_genes(adata, min_cells=1)


    def compute_mitochondria_percentage(self, adata):
        """
        Adds the column 'percent_mito' to the adata.obs column.

        The 'percent_mito' is computed by summing the mitochondria genes expression values and dividing by the
        sum of all genes expressed in the cell.

        :param adata: A sample's scanpy anndata object
        return: None
        """
        mito_genes = [gene for gene in adata.var_names if gene.startswith('MT.') or gene.startswith('MT-')]
        adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)

    def compute_n_counts(self, adata):
        """
        Adds the column 'n_counts' to the adata.obs column.

        The 'n_counts' is computed by summing the expression of all genes for a given cell.

        :param adata: A sample's scanpy anndata object
        return: None
        """
        adata.obs['n_counts'] =  np.sum(adata.X, axis=1)


    def sample_filtering(self):

        for sample_key, adata in self.sample_dict.items():

            self.filter_cells(adata)

            self.filter_genes(adata)

            self.compute_mitochondria_percentage(adata)

            self.compute_n_counts(adata)

        return self.sample_dict


class Basic_Preprocessing(object):

    def __init__(self, sample_dict, output_dir):

        self.sample_dict = sample_dict
        self.output_dir = output_dir

    def normalize(self, adata):
        """
        Normalize each cell by total counts over all genes, so that every cell
        has the same total count after normalization. Each cell has a total
        count equal to the median of the 'counts_per_cell', which is the sum
        of a cell's expression data, after normalization.

        Adds the column 'n_counts' to the adata.obs column. The 'n_counts' is
        computed by summing the expression of all genes for a given cell.

        References
        https://scanpy.readthedocs.io/en/latest/api/scanpy.api.pp.normalize_per_cell.html?highlight=normalize_per_cell

        :param adata: A sample's scanpy anndata object
        return: A sample's updated scanpy anndata object with normalized version of the original
        """
        sc.pp.normalize_per_cell(adata)

    def filter_genes_dispersion(self, sample_key, adata):
        """
        Extract highly variable genes.

        Calculates the average expression and dispersion for each gene, places these genes into bins,
        and then calculates a z-score for dispersion within each bin. This helps control for the relationship
        between variability and average expression.

        References
        https://scanpy.readthedocs.io/en/latest/api/scanpy.api.pp.filter_genes_dispersion.html?highlight=filter_genes_dispersion
        https://satijalab.org/seurat/pbmc3k_tutorial.html

        Macosko
        Volume 161, ISSUE 5, P1202-1214, May 21, 2015
        Highly Parallel Genome-wide Expression Profiling of Individual Cells Using Nanoliter Droplets

        :param sample_key: A string that is the sample's dictionary key to access self.sample_dict
        :param adata: A sample's normalized scanpy anndata object
        return: A sample's updated scanpy anndata object with low varying genes removed

        """

        filter_result = sc.pp.filter_genes_dispersion(adata.X, flavor='seurat') # Select highly variable genes
        self.sample_dict[sample_key] = adata[:, filter_result.gene_subset] # subset the genes and set the filtered anndata object as the dictionary value

    def log_transform(self, adata):
        """
        Logarithmize the data matrix plus one and save the logarithmized raw data for differential
        expression testing and plotting.

        References
        https://scanpy.readthedocs.io/en/latest/api/scanpy.api.pp.log1p.html?highlight=log

        :param adata: A sample's normalized and low variable gene filtered scanpy anndata object
        return: A sample's updated log transformed scanpy anndata obect
        """
        sc.pp.log1p(adata, copy=True).write(self.output_dir)

    def regress_out_unwanted_variation(self, adata):
        """
        Regress out unwanted sources of variation. Uses linear regression. This is inspired by Seurat’s regressOut function in R.

        References
        https://scanpy.readthedocs.io/en/latest/api/scanpy.api.pp.regress_out.html?highlight=regress_out
        https://satijalab.org/seurat/pbmc3k_tutorial.html

        :param adata: A sample's log transformed and low variable gene filtered scanpy anndata object
        return: A sample's updated scanpy anndata object with the corrected data matrix
        """
        sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])

    def sample_preprocessing(self):

        for sample_key, adata in self.sample_dict.items():

            self.normalize(adata)
            self.filter_genes_dispersion(sample_key, adata)
            self.normalize(adata) # renormalize after filtering
            self.log_transform(adata)
            self.regress_out_unwanted_variation(adata)

        return self.sample_dict

class Dimensionality_Reduction(object):

    def __init__(self, sample_dict):
        self.sample_dict = sample_dict

    def scale_data(self, adata):
        """
        Standardize expression features by removing the mean and scaling to unit variance.

        References
        http://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html

        :param adata: A sample's scanpy anndata object
        return: A sample's scanpy anndata matrix scaled
        """
        adata = StandardScaler().fit_transform(adata.X)

    def compute_principle_component_analysis(self, adata):
        """
        Computes PCA coordinates, loadings and variance decomposition. Uses the implementation of scikit-learn.

        References
        https://scanpy.readthedocs.io/en/latest/api/scanpy.api.tl.pca.html?highlight=pca

        :param adata: A sample's scaled scanpy anndata object

        return: Updates adata with the following fields.

        X_pca (.obsm) – PCA representation of data.
        PCs (.varm) – The principal components containing the loadings.
        variance_ratio (.uns[‘pca’]) – Ratio of explained variance.
        variance (.uns[‘pca’]) – Explained variance, equivalent to the eigenvalues of the covariance matrix.
        """
        sc.tl.pca(adata, n_comps=100)

    def sample_dimension_reduction(self):
        for sample_key, adata in self.sample_dict.items():
            self.scale_data(adata)
            self.compute_principle_component_analysis(adata)
        return self.sample_dict

class Cluster(object):

    def __init__(self, sample_dict, output_dir):

        self.sample_dict = sample_dict
        self.output_dir = output_dir

    def compute_tSNE(self, adata):
        """
        Compute t-distributed stochastic neighborhood embedding. The scanpy implementation uses scikit-learn.

        :param adata: A sample's dimension reduced scanpy object
        return: Updates adata with the following fields.
        X_tsne (np.ndarray (adata.obs, dtype float)) – tSNE coordinates of data.
        """
        sc.tl.tsne(adata)

    def compute_neighborhood_graph(self, adata):
        """
        Compute a neighborhood graph of observations.

        References
        https://scanpy.readthedocs.io/en/latest/api/scanpy.api.pp.neighbors.html?highlight=neighbors

        :param adata: A sample's dimension reduced scanpy object
        return: Updates adata with the following fields.
        connectivities (sparse matrix (.uns[‘neighbors’], dtype float32)) – Weighted adjacency matrix of the neighborhood graph of data points. Weights should be interpreted as connectivities.

        distances (sparse matrix (.uns[‘neighbors’], dtype float32)) – Instead of decaying weights, this stores distances for each pair of neighbors.
        """
        sc.pp.neighbors(adata)

    def compute_umap(self, adata):
        """
        Embed the neighborhood graph using UMAP.

        References
        https://scanpy.readthedocs.io/en/latest/api/scanpy.api.tl.umap.html?highlight=umap

        :param adata: A sample's dimension reduced scanpy object
        return: Updates adata with the following fields.
        X_umap (adata.obsm) – UMAP coordinates of data.
        """
        sc.tl.umap(adata)

    def compute_louvain_clustering(self, adata):
        """
        Cluster cells using the Louvain algorithm.

        References
        https://scanpy.readthedocs.io/en/latest/api/scanpy.api.tl.louvain.html?highlight=louvain

        :param adata: A sample's dimension reduced scanpy object
        return: Updates adata with the following fields.
        None – By default (copy=False), updates adata with the following fields:

        louvain : pandas.Series (adata.obs, dtype category)

            Array of dim (number of samples) that stores the subgroup id (‘0’, ‘1’, …) for each cell.

        AnnData – When copy=True is set, a copy of adata with those fields is returned.
        """
        sc.tl.louvain(adata)

    def compute_louvain_cluster_gene_rankings(self, adata):
        """
        Rank genes for louvain clusters.

        References
        https://scanpy.readthedocs.io/en/latest/api/scanpy.api.tl.rank_genes_groups.html?highlight=rank_genes_groups

        :param adata: A sample's dimension reduced scanpy object
        return: Updates adata with the following fields.
        names (structured np.ndarray (.uns[‘rank_genes_groups’])) – Structured array to be indexed by group id storing the gene names. Ordered according to scores.

        scores (structured np.ndarray (.uns[‘rank_genes_groups’])) – Structured array to be indexed by group id storing the score for each gene for each group. Ordered according to scores.

        logfoldchanges (structured np.ndarray (.uns[‘rank_genes_groups’])) – Structured array to be indexed by group id storing the log2 fold change for each gene for each group. Ordered according to scores. Only provided if method is ‘t-test’ like.
        """
        sc.tl.rank_genes_groups(adata, 'louvain')

    def sample_clustering(self):

        for sample_key, adata in self.sample_dict.items():
            self.compute_tSNE(adata)
            self.compute_neighborhood_graph(adata)
            self.compute_umap(adata)
            self.compute_louvain_clustering(adata)
            self.compute_louvain_cluster_gene_rankings(adata)

            # This will need to be updated if the script is planning to be scaled
            adata.write(self.output_dir)
