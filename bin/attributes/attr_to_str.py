"""
Turn a column into a string data type by prepending a specified string
to it.

Input is a attribute.tab file, columns are attribute names and rows
are cell ids/ sample ids. And a prepend.tab, where the first column is
attribute names that need the data type changed and the second column
is the string to prepend to the previous values.

e.g. if specified in prepend.tab file my_attr\tcluster 
my_attr -> my_attr
1       -> cluster_1
2       -> cluster_2
3       -> cluster_3

example usage:
    python tostr.py -i attr.tab -o attr_strd.tab -p prepend.tab
"""
import argparse
import sys
import pandas as pd
import numpy as np


def parse_args():

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i',"--in_file", type=str, help=""
                        )

    parser.add_argument('-p',"--prepend_file", type=str, help=""
                        )

    parser.add_argument('-o',"--out_file", type=str, help="",
                        )

    opts = parser.parse_args()
    in_file, out_file, prepend_file = opts.in_file, opts.out_file, opts.prepend_file

    return in_file, out_file, prepend_file


def read_data(in_file, prepend_file):
    adf = pd.read_table(in_file, index_col=0)
    prepend_df = pd.read_table(prepend_file, header=None) 
    prepend_dict = dict(zip(prepend_df[0], prepend_df[1]))
    return adf, prepend_dict


def data_transform(data_obj):
    adf, prepend_dict = data_obj
    def prependit(x):
        if np.isnan(x):
            return x
        else:
            return prepender + "_" + str(int(x))

    for key in prepend_dict.keys():
        prepender = prepend_dict[key]

        a = list(map(prependit, adf[key].tolist()))
        adf[key] = a
         
    return adf


def write_data(out_file, data_obj):
    data_obj.to_csv(out_file, sep="\t")


def main():

    in_file, out_file, rename_file = parse_args()

    data = read_data(in_file, rename_file)

    new_data = data_transform(data)

    write_data(out_file, new_data)


if __name__ == "__main__":
    sys.exit(main())
