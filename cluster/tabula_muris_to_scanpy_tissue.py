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
import tables


# In[40]:


tm_ds = "droplet"
tm_ds = "facs"


# In[41]:


# read the h5ad file
tm_droplet = sc.read("/projects/sysbio/users/cellAtlas/data/primary/mouse/tabula_muris/TM_"+tm_ds+"_mat.h5ad")
print(tm_droplet.X.shape)


# In[19]:


# read the meta data file
meta_df = pd.read_csv("/projects/sysbio/users/cellAtlas/data/primary/mouse/tabula_muris/TM_"+tm_ds+"_metadata.csv"
                     ,sep = ",", header=0, index_col=0,low_memory=False)
print(meta_df.shape)
# make sure the index is the same in metadata and expression matrix
print(set(meta_df.index == tm_droplet.obs.index))


# In[42]:


for col in meta_df.columns:
    tm_droplet.obs[col] = meta_df[col]
print(tm_droplet)


# In[43]:


tissue_set = set(tm_droplet.obs["tissue"])
for t in tissue_set:
    cells_in_tissue = tm_droplet.obs.index[tm_droplet.obs["tissue"]==t]
    cells_in_tissue_index = [tm_droplet.obs.index.tolist().index(x) for x in cells_in_tissue]
    print(t)
    #print(cells_in_tissue)
    print(len(cells_in_tissue_index))
    tm_droplet_tissue = tm_droplet[cells_in_tissue_index,]
    #print(tm_droplet_tissue)
    #print(set(tm_droplet_tissue.obs["tissue"]))
    
    tm_droplet_tissue.write("/projects/sysbio/users/cellAtlas/data/scanpyObj/mouse_tabulaMuris_"+tm_ds+"_"+t+".h5ad")


# In[ ]:




