
# coding: utf-8

# In[1]:


import scanpy.api as sc
import pandas as pd


# In[2]:


hematopoietic = sc.read('/projects/sysbio/users/cellAtlas/data/human/hematopoietic/GSE79331_normalized_minus2ndcol.tab').T


# In[4]:


hematopoietic.write('/projects/sysbio/users/cellAtlas/scanpyObjects/hematopoietic.h5ad')

