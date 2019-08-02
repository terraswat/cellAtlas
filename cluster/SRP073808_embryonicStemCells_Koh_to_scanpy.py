#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import pandas as pd
import h5py
import scanpy.api as sc
import collections
import scipy.sparse as sp_sparse
import numpy as np
#import tables


# In[2]:


# read the expression matrix
exprMatrix = sc.read("/projects/sysbio/users/cellAtlas/data/primary/human/SRP073808_embryonicStemCells_Koh/SRP073808_genelevel_counts.tsv")

# transpose into genes as columns (vars), samples / cells as rows (obs)
exprMatrix = exprMatrix.T
print(exprMatrix.shape)


# In[5]:


exprMatrix.obs.index
#exprMatrix.var.index


# In[10]:


# read the metadata
meta_df = pd.read_csv("/projects/sysbio/users/cellAtlas/data/primary/human/SRP073808_embryonicStemCells_Koh/SraRunTable.txt"
                     ,sep = "\t", header=0, index_col=11)

print(meta_df.shape)
#print(len(set(meta_df.index)))
#print(meta_df)

# some samples are missing in the expression matrix
meta_df = meta_df.loc[exprMatrix.obs.index]
print(meta_df.shape)

# make sure the index is the same in metadata and expression matrix
print(set(meta_df.index == exprMatrix.obs.index))


# In[11]:


# save metadata in scanpy object's obs dataframe
for col in meta_df.columns:
    exprMatrix.obs[col] = meta_df[col]
#print(exprMatrix.obs)
#print(exprMatrix.obs.index)
print(exprMatrix)


# In[8]:


# write scanpy object to file
exprMatrix.write("/projects/sysbio/users/cellAtlas/data/scanpyObj/SRP073808_embryonicStemCells_Koh.h5ad")


# In[ ]:




