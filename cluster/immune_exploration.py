
# coding: utf-8

# In[1]:


import pandas as pd
import scanpy.api as sc


# In[ ]:


ica_cord_blood = sc.read_10x_h5(glob.glob('/projects/sysbio/users/cellAtlas/data/human/immune_census/ica_cord_blood_h5.h5')[0], "GRCh38")
ica_cord_blood.write('/projects/sysbio/users/cellAtlas/scanpyObjects/ica_cord_blood.h5ad')


# In[ ]:


ica_bone_marrow = sc.read_10x_h5(glob.glob('/projects/sysbio/users/cellAtlas/data/human/immune_census/ica_bone_marrow_h5.h5')[0], "GRCh38")
ica_bone_marrow.write('/projects/sysbio/users/cellAtlas/scanpyObjects/ica_bone_marrow.h5ad')

