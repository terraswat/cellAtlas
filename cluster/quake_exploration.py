
# coding: utf-8

# In[1]:


import scanpy.api as sc
import pandas as pd


# In[2]:


quake = sc.read('/projects/sysbio/users/cellAtlas/data/human/quake_BrainGeo1/sc000JAJ.tsv').T


# In[19]:


quake_meta = pd.read_csv('/projects/sysbio/users/cellAtlas/data/human/joint/meta.tsv', sep='\t')


# In[16]:


for c in quake_meta.columns.values.tolist():
    quake.obs[c] = list(quake_meta[c])


# In[ ]:


quake.write('/projects/sysbio/users/cellAtlas/scanpyObjects/quakeBrainGeo1.h5ad')

