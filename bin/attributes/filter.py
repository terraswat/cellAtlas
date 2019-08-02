"""
Remove attributes from a .tab file.

example usage:
    python filter.py -i scanpyObj.h5ad -o attr.tab -k keep.list

keep.list is a single column file with an attirbute name per line.

**must have at least two columns that you want to keep
"""
import argparse
import sys
import pandas as pd
import numpy as np
import scanpy.api as sc


def parse_args():

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i',"--in_file", type=str, help=""
                        )


    parser.add_argument('-k',"--keep_file", type=str, help=""
                        )

    parser.add_argument('-o',"--out_file", type=str, help="",
                        default="full_sim.euc.tab"
                        )

    opts = parser.parse_args()
    in_file, out_file, keep_file = opts.in_file, opts.out_file, opts.keep_file

    return in_file, out_file, keep_file


def read_data(in_file, keep_file):
    adf = sc.read(in_file).obs
    keep_list = np.loadtxt(keep_file, dtype=str, delimiter="\t").tolist() 
    return adf, keep_list


def data_transform(data_obj):
    adf, keep_list = data_obj
    return adf[keep_list]

def write_data(out_file, data_obj):
    data_obj.to_csv(out_file, sep="\t")


def main():

    in_file, out_file, keep_file = parse_args()

    data = read_data(in_file, keep_file)

    new_data = data_transform(data)

    write_data(out_file, new_data)


if __name__ == "__main__":
    sys.exit(main())
