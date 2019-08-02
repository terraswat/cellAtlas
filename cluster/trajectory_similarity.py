
import os
import pandas as pd
from collections import defaultdict
import numpy as np

path='/projects/sysbio/users/cellAtlas/data/trajectories/sim_all/'
dfs=[]
final_df=pd.DataFrame()
indecies = set()
for fn in os.listdir(path):
    df=pd.read_csv(path+fn, sep='\t', index_col=0)
    indecies_ = set(df.index)
    indecies = indecies.union(indecies_)

final_df = pd.DataFrame(index=indecies)

for fn in os.listdir(path):
    df=pd.read_csv(path+fn, sep='\t', index_col=0)
    if df.shape[0]>30:
        final_df[df.columns[0]]=df[df.columns[0]]

print(final_df.shape)
print(final_df.isna().sum().sum())
print(final_df.max().max())
final_df = final_df.fillna(final_df.max().max())
#intr=set(final_df.index).intersection(final_df.columns)
print(final_df.isnull().values.any())
#final_df=final_df.loc[intr,intr]
#final_df=final_df.iloc[0:10,0:10]
final_df.to_csv(path+'trajectory_dist.tab',sep='\t')
#final_df=final_df.drop([c for c in final_df.columns.values.tolist() if 'WISH' in c],axis=1)
#final_df=final_df.drop([r for r in final_df.index.values.tolist() if 'WISH' in r],axis=0)

n= final_df.shape[0]

#upper_tr=np.triu_indices(n,1)
lower_tr=np.tril_indices(n,-1)

tri = final_df.values[lower_tr] * final_df.T.values[lower_tr]

#final_df.values[upper_tr] = tri
final_df.values[lower_tr] = tri
final_df=final_df.T
final_df.values[lower_tr]=tri

np.fill_diagonal(final_df.values, 0)
assert np.allclose(final_df, final_df.T)
final_df.to_csv(path+'trajectory_dist_sym.tab', sep='\t')
