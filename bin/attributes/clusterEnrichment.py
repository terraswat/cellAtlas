"""
Takes in an attribute file that has a "louvain" cluster attribute.

"""
import argparse
import sys
import os
import pandas as pd
from oneByAllAttrsStats import oneByAllStats
from data_type_dict import data_type_dictionary

def parse_args():

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i',"--in_file", type=str, help=""
                        )

    parser.add_argument('-o',"--out_dir", type=str, help="",
                        )

    opts = parser.parse_args()
    in_file, out_dir = opts.in_file, opts.out_dir

    return in_file, out_dir


def read_data(in_file):
    df = pd.read_table(in_file, index_col=0)
    cluster_attr = df["louvain"].astype("category")
    cluster_dummy_df = pd.get_dummies(cluster_attr)
    df = df.drop("louvain", axis=1)
    return cluster_dummy_df, df


def data_transform_and_output(attrdf, clusterdf, out_dir):
    data_type_dict = data_type_dictionary(attrdf)
    
    for column in clusterdf:
        oneByAll = oneByAllStats(attrdf, data_type_dict, clusterdf[column], "bin")
        oneByAll.name = "cluster_%d"%column
        out_file = os.path.join(out_dir, "cluster_%d"%column)
        oneByAll.to_csv(out_file, sep="\t")


def main():

    in_file, out_dir = parse_args()

    cluster_dummy_df, attrdf = read_data(in_file)

    new_data = data_transform_and_output(attrdf, cluster_dummy_df, out_dir)


if __name__ == "__main__":
    sys.exit(main())
