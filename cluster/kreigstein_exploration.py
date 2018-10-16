
# coding: utf-8

# In[6]:


import scanpy.api as sc
import pandas as pd


# In[23]:


kreigstein = sc.read('/projects/sysbio/users/cellAtlas/data/human/kriegstein_jointCirmBrain1/sc002MIV.tsv').T


# In[24]:


kreigstein_meta = pd.read_csv('/projects/sysbio/users/cellAtlas/data/human/kriegstein_jointCirmBrain1/meta.tsv', sep='\t')


# In[27]:


for c in kreigstein_meta.columns.values.tolist():
    kreigstein.obs[c] = list(kreigstein_meta[c])


# In[29]:


kreigstein.write('/projects/sysbio/users/cellAtlas/scanpyObjects/kriegstein_jointCirmBrain1.h5ad')

