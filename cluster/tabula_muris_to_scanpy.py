#!/usr/bin/env python
# coding: utf-8

# In[17]:


import glob
import pandas as pd
import h5py
import scanpy.api as sc
import collections
import scipy.sparse as sp_sparse
import numpy as np
#import tables


# ## droplet data

# In[18]:


# read the h5ad file
tm_droplet = sc.read("/projects/sysbio/users/cellAtlas/data/mouse/tabula_muris/TM_droplet_mat.h5ad")
print(tm_droplet.X.shape)


# In[19]:


# read the meta data file
meta_df = pd.read_csv("/projects/sysbio/users/cellAtlas/data/mouse/tabula_muris/TM_droplet_metadata.csv"
                     ,sep = ",", header=0, index_col=0,low_memory=False)
print(meta_df.shape)
# make sure the index is the same in metadata and expression matrix
print(set(meta_df.index == tm_droplet.obs.index))


# In[20]:


# save metadata in scanpy object's obs dataframe
for col in meta_df.columns:
    tm_droplet.obs[col] = meta_df[col]
#print(tm_droplet.obs)
#print(tm_droplet.obs.index)
print(tm_droplet)


# In[21]:


# write scanpy object to file
tm_droplet.write("/projects/sysbio/users/cellAtlas/scanpyObjects/mouse_tabulaMuris_droplet.h5ad")


# ## facs data

# In[22]:


# read the h5ad file
tm_facs = sc.read("/projects/sysbio/users/cellAtlas/data/mouse/tabula_muris/TM_facs_mat.h5ad")
print(tm_facs.X.shape)


# In[23]:


# read the meta data file
meta_df = pd.read_csv("/projects/sysbio/users/cellAtlas/data/mouse/tabula_muris/TM_facs_metadata.csv"
                     ,sep = ",", header=0, index_col=0,low_memory=False)
print(meta_df.shape)
# make sure the index is the same in metadata and expression matrix
print(set(meta_df.index == tm_facs.obs.index))


# In[24]:


# save metadata in scanpy object's obs dataframe
for col in meta_df.columns:
    tm_facs.obs[col] = meta_df[col]
#print(tm_facs.obs)
#print(tm_facs.obs.index)
print(tm_facs)


# In[25]:


# write scanpy object to file
tm_facs.write("/projects/sysbio/users/cellAtlas/scanpyObjects/mouse_tabulaMuris_facs.h5ad")


# ## check out gene ortholog mapping

# In[28]:


ortho_map = pd.read_csv("/projects/sysbio/users/cellAtlas/data/mouse_human_orthologs.txt"
                     ,sep = "\t", header=0)
print(ortho_map.shape) # (149657, 7)
print(ortho_map.columns)
# Index(['Gene stable ID', 'Transcript stable ID', 'Gene name',
#       'Human gene stable ID', 'Human gene name',
#       'Human protein or transcript stable ID',
#       'Human orthology confidence [0 low, 1 high]'],
#      dtype='object')


# In[36]:


mouse_genes = set(tm_facs.var.index).intersection(set(ortho_map['Gene name']))
print(len(mouse_genes)) #20954
print(len(tm_facs.var.index)) #23433


# In[44]:


ortho_map_overlap = ortho_map.loc[ortho_map["Gene name"].isin(mouse_genes)]
print(ortho_map_overlap.shape) #(99537, 7)
print(len(set(ortho_map_overlap["Gene name"]))) #20954
print(len(set(ortho_map_overlap["Human gene name"]))) #16311
print(len(set(ortho_map_overlap["Human gene stable ID"]))) #16335


# In[ ]:




