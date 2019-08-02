#!/usr/bin/env python
# coding: utf-8

# In[168]:


import glob
import pandas as pd
import h5py
import scanpy.api as sc
import collections
import scipy.sparse as sp_sparse
import numpy as np
import tables


# In[169]:


# pick one data set here, but run code through for both
pick_one_dataset = "ica_cord_blood"
pick_one_dataset = "ica_bone_marrow"

scanpy_object = sc.read_10x_h5(glob.glob('/projects/sysbio/users/cellAtlas/data/primary/human/immune_census/'+pick_one_dataset+'_h5.h5')[0], "GRCh38")

print(scanpy_object.shape)


# # Explore the sample IDs

# In[170]:


# Donor ID is encoded in sample IDs
donor_id = [y.split("Manton")[1] for y in [x.split("_")[0] for x in scanpy_object.obs.index]]
#print(set(donor_id))


# In[171]:


# Sequencer is the same for all samples
sequencer = [x.split("_")[1] for x in scanpy_object.obs.index]
#print(set(sequencer))


# In[172]:


# Biomaterial IDis encoded in sample IDs

# split the barcode part of sample IDs
barcode = [x.split("_")[2] for x in scanpy_object.obs.index]
before_barcode = [x.split("-")[0] for x in barcode]
#after_barcode = [x.split("-")[2] for x in barcode]
#barcode = [x.split("-")[1] for x in barcode]

# connect donor ID with number before barcode to get the biomaterial ID
biomaterial_id = [y+"_"+x for x,y in zip(donor_id,before_barcode)]


# # Read metadata tables
# 

# In[173]:


# read the metadata 
# get rid of second line containing descriptions of the columns
# get rid of any column containing only NA (aka. empty columns)
# make column names unique (some are named the same in different meta data files, but contain different values)

meta_donor = pd.read_csv("/projects/sysbio/users/cellAtlas/data/primary/human/immune_census/metadata_tables/donor_organism.txt"
                     ,sep = "\t", header=0, index_col=0)
meta_donor.drop(meta_donor.index[0], inplace=True)
meta_donor.dropna(axis=1, how='all', inplace=True)
meta_donor.columns = ["donor."+x for x in meta_donor.columns]
meta_donor.index.name = "donor."+meta_donor.index.name

meta_specimen = pd.read_csv("/projects/sysbio/users/cellAtlas/data/primary/human/immune_census/metadata_tables/specimen_from_organism.txt"
                     ,sep = "\t", header=0, index_col=0)
meta_specimen.drop(meta_specimen.index[0], inplace=True)
meta_specimen.dropna(axis=1, how='all', inplace=True)
meta_specimen.columns = ["specimen."+x for x in meta_specimen.columns]
meta_specimen.index.name = "specimen."+meta_specimen.index.name

meta_cellsuspension = pd.read_csv("/projects/sysbio/users/cellAtlas/data/primary/human/immune_census/metadata_tables/cell_suspension.txt"
                     ,sep = "\t", header=0, index_col=0)
meta_cellsuspension.drop(meta_cellsuspension.index[0], inplace=True)
meta_cellsuspension.dropna(axis=1, how='all', inplace=True)
meta_cellsuspension.columns = ["cellsuspension."+x for x in meta_cellsuspension.columns]
meta_cellsuspension.index.name = "cellsuspension."+meta_cellsuspension.index.name

meta_sequencing = pd.read_csv("/projects/sysbio/users/cellAtlas/data/primary/human/immune_census/metadata_tables/sequencing_process.txt"
                     ,sep = "\t", header=0, index_col=0)
meta_sequencing.drop(meta_sequencing.index[0], inplace=True)
meta_sequencing.dropna(axis=1, how='all', inplace=True)
meta_sequencing.columns = ["sequencing."+x for x in meta_sequencing.columns]
meta_sequencing.index.name = "sequencing."+meta_sequencing.index.name


# In[174]:


# collect all metadata in a dataframe later to be added to the scanpy object obs
metadata_df = pd.DataFrame(index=scanpy_object.obs.index)

# set the donor id as a column in the metadata dataframe
metadata_df["donor_id"] = donor_id

# set the biomaterial id as a column in the metadata dataframe
metadata_df["biomaterial_id"] = biomaterial_id

print(metadata_df.shape)


# In[175]:


# meta_donor
# expand each column to fit the samples and add to metadata dataframe
for col in meta_donor.columns:
    expanded_list = [meta_donor.loc[row,col] for row in donor_id]
    metadata_df[col] = expanded_list
print(metadata_df.shape)


# In[176]:


# meta_specimen
# expand each column to fit the samples and add to metadata dataframe
for col in meta_specimen.columns:
    expanded_list = [meta_specimen.loc[row,col] for row in biomaterial_id]
    metadata_df[col] = expanded_list
print(metadata_df.shape)


# In[177]:


# meta_cellsuspension & meta_sequencing

# make index match biomaterial ID
meta_cellsuspension.index = [x.split("_")[0]+"_"+x.split("_")[1] for x in meta_cellsuspension.index]

# include the information from meta_sequencing in meta_cellsuspension by mapping the sequencing process id
sequencing_process_id = [x.split("||")[1] for x in meta_cellsuspension["cellsuspension.process_ids"]]
meta_cellsuspension["cellsuspension.seq_process_id"] = [x.split("||")[1] for x in meta_cellsuspension["cellsuspension.process_ids"]]
for col in meta_sequencing.columns:
    expanded_list = [meta_sequencing.loc[row,col] for row in sequencing_process_id]
    meta_cellsuspension[col] = expanded_list
print(meta_cellsuspension.shape)

# expand each column to fit the samples and add to metadata dataframe
for col in meta_cellsuspension.columns:
    expanded_list = [meta_cellsuspension.loc[row,col] for row in biomaterial_id]
    metadata_df[col] = expanded_list
print(metadata_df.shape)


# # Connect metadata to scanpy object

# In[178]:


# save metadata in scanpy object's obs dataframe
for col in metadata_df.columns:
    scanpy_object.obs[col] = metadata_df[col]
print(scanpy_object)


# In[179]:


# write to file
scanpy_object.write("/projects/sysbio/users/cellAtlas/data/scanpyObj/"+pick_one_dataset+".h5ad")


# In[ ]:




