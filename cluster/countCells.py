import scanpy.api as sc
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--INPUT', default=None, required=True, help='input of scnapy object absolute path')
args = parser.parse_args()


obj=sc.read(args.INPUT)
matrix=obj.X
n_cells=matrix.shape[0]

abs_path=os.path.abspath(args.INPUT)
print('{}\t{} number of cells'.format(abs_path,n_cells))

