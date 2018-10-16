#!/usr/bin/env python
# coding: utf-8

# In[2]:


import glob
import pandas as pd
import h5py
import scanpy.api as sc
import collections
import scipy.sparse as sp_sparse
import numpy as np
#import tables


# In[12]:


# run once for each platform ID
platform_id = "GPL18573"
platform_id = "GPL16791"


# In[13]:


# read the expression matrix
exprMatrix = sc.read("/projects/sysbio/users/cellAtlas/data/human/GSE81076_pancreas_Grun/GSE81076-"+platform_id+"_genelevel_counts.tsv")

# transpose into genes as columns (vars), samples / cells as rows (obs)
exprMatrix = exprMatrix.T
print(exprMatrix.shape)


# In[14]:


exprMatrix.obs.index
exprMatrix.var.index


# In[15]:


# read the metadata
meta_df = pd.read_csv("/projects/sysbio/users/cellAtlas/data/human/GSE81076_pancreas_Grun/GSE81076-"+platform_id+"_metadata.tsv"
                     ,sep = "\t", header=0, index_col=0)

print(meta_df.shape)

# make sure the index is the same in metadata and expression matrix
print(set(meta_df.index == exprMatrix.obs.index))


# In[16]:


# save metadata in scanpy object's obs dataframe
for col in meta_df.columns:
    exprMatrix.obs[col] = meta_df[col]
#print(exprMatrix.obs)
#print(exprMatrix.obs.index)
print(exprMatrix)


# In[17]:


# write scanpy object to file
exprMatrix.write("/projects/sysbio/users/cellAtlas/scanpyObjects/GSE81076-"+platform_id+"_pancreas_Grun.h5ad")


# In[ ]:




