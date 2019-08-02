import seaborn as sns
import pandas as pd
import matplotlib as mpl
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist,squareform
import numpy as np
import argparse
import pandas as pd
from matplotlib.colors import ListedColormap
from sklearn.preprocessing import LabelEncoder
from matplotlib.pyplot import cm
from scipy.cluster.hierarchy import fcluster
encoder1= LabelEncoder()

parser = argparse.ArgumentParser()
parser.add_argument('-attr', '--attributes_fn', default=None, required=True, help='tab delim attr files')
parser.add_argument('-feature', '--feature_fn', default=None, required=True, help='tab delim feature file')
parser.add_argument('-o', '--out', default=None, required=True, help='output for png')
args = parser.parse_args()

feature=args.feature_fn
attributes=args.attributes_fn
#features='/projects/sysbio/users/cellAtlas/data/example.formats/sim.tab'
#attributes='/projects/sysbio/users/cellAtlas/data/example.formats/attr.tab'

features_df=pd.read_csv(feature, sep='\t', index_col=0)
condensed_features=pdist(features_df.values) #condensed distance matrix
attr_df=pd.read_csv(attributes, sep='\t', index_col=0)
attr_df=attr_df.loc[features_df.index.values.tolist()]
#print(attr_df.shape)
#print(attr_df.head())
attr_df['algorithm']=encoder1.fit_transform(attr_df['algorithm'].values.tolist())
encoder2=LabelEncoder()
attr_df['Dataset']=encoder2.fit_transform(attr_df['Dataset'].values.tolist())


cluster=list(set(attr_df.cluster.values.tolist()))
algorithm=list(set(attr_df.algorithm.values.tolist()))
dataset=list(set(attr_df.Dataset.values.tolist()))

#colors1=sns.color_palette("RdBu_r", len(cluster)).as_hex()
colors2=sns.color_palette("BrBG", len(algorithm)).as_hex()
#colors3=['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6']
colors3=sns.color_palette('deep',len(dataset)).as_hex()

#cluster_colors=ListedColormap(colors1)
algorithm_colors=ListedColormap(colors2)
dataset_colors=ListedColormap(colors3)

#c_colors=dict(zip(cluster,colors1))
a_colors=dict(zip(algorithm,colors2))
d_colors=dict(zip(dataset,colors3))

for i in attr_df.index.values.tolist():
    #louvain=attr_df.loc[i,'cluster'] #louvain of sample
    algorithm=attr_df.loc[i,'algorithm'] #algorithm
    dataset=attr_df.loc[i,'Dataset']
        
    #l_color=c_colors[louvain]
    a_color=a_colors[algorithm]
    dataset_color=d_colors[dataset]

    #attr_df.loc[i,'cluster_colors']=((l_color))
    attr_df.loc[i,'algorithm_colors']=((a_color))
    attr_df.loc[i,'dataset_colors']=((dataset_color))
 
# Calculate the distance between each sample
Z = hierarchy.linkage(condensed_features, 'ward')
panel_width=10/20
panel_height=9/10
plt.figure(figsize=(20,10))  #units are always inches
cmap = cm.rainbow(np.linspace(0, 1, 6))
hierarchy.set_link_color_palette([mpl.colors.rgb2hex(rgb[:3]) for rgb in cmap])
cluster_cutoff=.2*max(Z[:,2]) #cutoff that defines clusters
cluster_assignments = fcluster(Z,cluster_cutoff,'distance')

panel=plt.axes([0.053,0.04,panel_width,panel_height]) #this adds axis labels)
panel.plot([cluster_cutoff,cluster_cutoff],[0,1530],linewidth=2,zorder=3,color='r')
d=hierarchy.dendrogram(Z, labels=features_df.index, leaf_rotation=45, orientation="left",distance_sort=False,color_threshold=cluster_cutoff,get_leaves=True)
panel.tick_params(axis='both',which='both', bottom=True, right=False, left=False, top=False, labelbottom=True,labelleft=False,labeltop=False,labelright=False)

#print('len of leaves',len(d['color_list']))

trajectory_annotation=attr_df.loc[list(d['ivl'][::-1])]
trajectory_annotation.index.names = ['']
attr_df['Dataset']=encoder2.inverse_transform(attr_df['Dataset'].values.tolist())
attr_df['algorithm']=encoder1.inverse_transform(attr_df['algorithm'].values.tolist())

trajectory_annotation_legend=attr_df.loc[list(d['ivl'][::-1])]
#print(trajectory_annotation.head(10))
#panel1=plt.axes([0.55,0.04,2/20,9/10]) #this adds axis labels)
#heatmap1=sns.heatmap(trajectory_annotation[['cluster']],annot=False,cbar=False,yticklabels=False,cmap=cluster_colors)
#panel1.tick_params(axis='both',which='both', bottom=True, right=False, left=False, top=False, labelbottom=True,labelleft=False,labeltop=False,labelright=False)
#legend_cell_type = [mpatches.Patch(color=k, label=v) for k,v in attr_df[['cluster_colors','cluster']].drop_duplicates().values]
#panel1.legend(loc='upper left',bbox_to_anchor=(2.01,0.9),handles=legend_cell_type,frameon=True).set_title(title='cluster',prop={'size':10})


#attr_df['biosample_cell_type']=encoder.inverse_transform(attr_df['biosample_cell_type'].values.tolist())
panel2=plt.axes([0.65,0.04,2/20,9/10]) #this adds axis labels)
heatmap2=sns.heatmap(trajectory_annotation[['algorithm']],annot=False,cbar=False,yticklabels=False,cmap=algorithm_colors)
panel2.tick_params(axis='both',which='both', bottom=True, right=False, left=False, top=False, labelbottom=True,labelleft=False,labeltop=False,labelright=False)
legend_cell_type = [mpatches.Patch(color=k, label=v) for k,v in trajectory_annotation_legend[['algorithm_colors','algorithm']].drop_duplicates().values]
panel2.legend(loc='upper left',bbox_to_anchor=(2.35,0.75),handles=legend_cell_type,frameon=True).set_title(title='Algorithm',prop={'size':10})

#plotting dataset
panel3=plt.axes([0.75,0.04,2/20,9/10]) #this adds axis labels)
heatmap3=sns.heatmap(trajectory_annotation[['Dataset']],annot=False,cbar=False,yticklabels=False,cmap=dataset_colors)
panel3.tick_params(axis='both',which='both', bottom=True, right=False, left=False, top=False, labelbottom=True,labelleft=False,labeltop=False,labelright=False)
legend_dataset = [mpatches.Patch(color=k, label=v) for k,v in trajectory_annotation_legend[['dataset_colors','Dataset']].drop_duplicates().values]
panel3.legend(loc='upper left',bbox_to_anchor=(1.03,0.65),handles=legend_dataset,frameon=True).set_title(title='Dataset',prop={'size':10})
plt.savefig(args.out,dpi=300)

#folowing block gets the color of the leafs in the dendogram and makes a dict of {color-->trajectories within the same color}. Currently not used
'''
def get_cluster_classes(den, label='ivl'):
    from collections import defaultdict
    cluster_idxs = defaultdict(list)
    for c, pi in zip(den['color_list'], den['icoord']):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))

    cluster_classes = defaultdict(list)
    for c, l in cluster_idxs.items():
        i_l = [den[label][i] for i in l]
        cluster_classes[c] = i_l

    return cluster_classes
#for cluster, items in get_cluster_classes(d).items():
#    print(cluster, items)
'''

#print(cluster_assignments)
cluster_assignments_filtered=[] #clusters with less than 5 samples are labeled as ungrouped
cluster_count=dict()
#counting how many samples are within cluster, and excluding clusters that have less than 5 samples
for cluster in cluster_assignments:
    cluster_count[cluster]=cluster_count.get(cluster,0)+1
for i in range(len(list(cluster_assignments))):
    cluster=cluster_assignments[i]
    if cluster_count[cluster]<5:
        cluster_assignments[i]=0
    
#print(cluster_count)
#print(cluster_assignments)
encoder4=LabelEncoder()
family_output = pd.DataFrame({'dataset':features_df.index, 'trajectory_families':cluster_assignments})
family_output['trajectory_families']=encoder4.fit_transform(family_output['trajectory_families'].values.tolist())
family_output=family_output.set_index('dataset')
family_output=family_output.loc[trajectory_annotation.index]
family_output.index.names = ['']
families=sorted(family_output.trajectory_families.unique())
colors1=sns.color_palette("coolwarm",len(families)).as_hex()
family_colors=ListedColormap(colors1)
f_colors=dict(zip(families,colors1))

for i in family_output.index.values.tolist():
    family=family_output.loc[i,'trajectory_families'] #
    family_color=f_colors[family]
    family_output.loc[i,'family_colors']=family_color
#print(family_output)

panel4=plt.axes([0.555,0.04,2/20,9/10]) #this adds axis labels)
heatmap4=sns.heatmap(family_output[['trajectory_families']],annot=False,cbar=False,yticklabels=False,cmap=family_colors) #when i set the cmap to string colormap, it works, but when i set it to family_colors it doesnt. WEIRD they have the same colors!!!:: from duncan:  probably an encoding thing.

panel4.tick_params(axis='both',which='both', bottom=True, right=False, left=False, top=False, labelbottom=True,labelleft=False,labeltop=False,labelright=False)
family_output['trajectory_families']=encoder4.inverse_transform(family_output['trajectory_families'].values.tolist())

legend_dataset = [mpatches.Patch(color=k, label=v) for k,v in family_output.sort_values(by=['trajectory_families'])[['family_colors','trajectory_families']].drop_duplicates().values]
panel4.legend(loc='upper left',bbox_to_anchor=(3.2,1),title='Trajectory Families',prop={'size':10},handles=legend_dataset)

plt.savefig(args.out,dpi=300)

#cluster_output.to_csv('PRELIMINARYtrajectoryClustering.tsv',sep='\t')

