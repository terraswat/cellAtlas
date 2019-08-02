
import os
with open('/projects/sysbio/users/cellAtlas/data/trajectories/trajectory_meta.tsv','w') as f:
    f.write('index\tDataset\talgorithm\tcluster\n')
    for fn in os.listdir('/projects/sysbio/users/cellAtlas/data/trajectories/sim_all'):
        if fn!='--monocle' and fn!='trajectory_dist.tab' and fn!='trajectory_dist_sym.tab':
            fname=fn.split('.')
            left=fname[0]
            right=fname[1]
            cluster_num=int(left[-1])
            left=left.replace('_cluster_%d'%cluster_num,'')
            
            dataset=left.replace('_cluster_','')
            for c in dataset:
                if c.isdigit():
                    dataset=dataset.replace(c,'')
            print(dataset)
            algorithm=right.replace('tsv--','') 
            
            f.write('{}\t{}\t{}\t{}\n'.format(fn,dataset,algorithm,cluster_num))
