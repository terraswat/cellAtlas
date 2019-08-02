
import  matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
from sklearn.preprocessing import LabelEncoder
encoder = LabelEncoder()
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-attr', '--attributes_fn', default=None, required=True, help='tab delim attr files')
parser.add_argument('-feature', '--feature_fn', default=None, required=True, help='tab delim feature file')
args = parser.parse_args()
'''
width=20
height=10
plt.figure(figsize=(width,height))  #units are always inches
########----------------------Canvas--------------------##########

panel_width=18
panel_height=9

panel=plt.axes([0.01,0.1,panel_width,panel_height]) #this adds axis labels
'''
#feature=args.feature_fn
#attributes=args.attributes_fn
features='/projects/sysbio/users/cellAtlas/data/example.formats/sim.tab'
attributes='/projects/sysbio/users/cellAtlas/data/example.formats/attr.tab'

features_df=pd.read_csv(features, sep='\t', index_col=0)
attr_df=pd.read_csv(attributes, sep='\t', index_col=0)
attr_df['biosample_cell_type']=encoder.fit_transform(attr_df['biosample_cell_type'].values.tolist())

louvain=set(attr_df.louvain.values.tolist())
cell_type=set(attr_df.biosample_cell_type.values.tolist())

colors1=sns.color_palette("RdBu_r", len(louvain)).as_hex()
colors2=sns.color_palette("BrBG", len(cell_type)).as_hex()

louvain_colors=dict(zip(louvain,colors1[:len(louvain)]))
cell_type_colors=dict(zip(cell_type,colors2[:len(cell_type)]))

for i in attr_df.index.values.tolist():
    louvain=attr_df.loc[i,'louvain'] #louvain of sample
    cell_type=attr_df.loc[i,'biosample_cell_type'] #cell type
    
    
    l_color=louvain_colors[louvain]
    c_color=cell_type_colors[cell_type]
    attr_df.loc[i,'louvain_colors']=l_color
    attr_df.loc[i,'cell_type_colors']=c_color
#print(attr_df.head())
    
#print(len(louvain)) #8 clusters
#print(len(cell_type)) #9 cell)_type
cmap = sns.diverging_palette(220, 20, as_cmap=True)
g = sns.clustermap(features_df, cmap=cmap,row_cluster=True, col_cluster=False,row_colors=attr_df[['louvain_colors','cell_type_colors']],linewidths=0, xticklabels=False, yticklabels=False)
#for some reason, it cant take two legends....
#legend_louvain = [mpatches.Patch(color=c, label=l) for c,l in attr_df[['louvain_colors','louvain']].drop_duplicates().values]
#l2=g.ax_heatmap.legend(loc='upper left',bbox_to_anchor=(0.05,1.3),handles=legend_louvain,frameon=True)
#l2.set_title(title='louvain cluster',prop={'size':10})
attr_df['biosample_cell_type']=encoder.inverse_transform(attr_df['biosample_cell_type'].values.tolist())
legend_cell_type = [mpatches.Patch(color=k, label=v) for k,v in attr_df[['cell_type_colors','biosample_cell_type']].drop_duplicates().values]
l1=g.ax_heatmap.legend(loc='upper left',bbox_to_anchor=(1.01,0.6),handles=legend_cell_type,frameon=True)
l1.set_title(title='cell type',prop={'size':10})

for label in attr_df.louvain.unique():
    g.ax_col_dendrogram.bar(0, 0, color=louvain_colors[label],
                            label=label, linewidth=0)
g.ax_col_dendrogram.legend(bbox_to_anchor=(1.22,0.01), ncol=1).set_title(title='louvain cluster',prop={'size':10})

# Adjust the postion of the main colorbar for the heatmap
g.cax.set_position([.01, .2, .03, .45])
g.savefig('/projects/sysbio/users/cellAtlas/plots/test/IAheatmap.png')


