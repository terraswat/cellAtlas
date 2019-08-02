"""
Script renames attributes in a .tab file.

example usage:
    python rename.py -i attr.tab -o renamed_attr.tab -r rename.tab

rename.tab is a two column file, where the first column is the original
attribute name and the second column is the new attribute name
"""
import argparse
import sys
import pandas as pd


def parse_args():

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i',"--in_file", type=str, help=""
                        )


    parser.add_argument('-r',"--rename_file", type=str, help=""
                        )

    parser.add_argument('-o',"--out_file", type=str, help="",
                        default="full_sim.euc.tab"
                        )

    opts = parser.parse_args()
    in_file, out_file, rename_file = opts.in_file, opts.out_file, opts.rename_file

    return in_file, out_file, rename_file


def read_data(in_file, rename_file):
    adf = pd.read_table(in_file, index_col=0)
    rename_df = pd.read_table(rename_file, header=None) 
    rename_dict = dict(zip(rename_df[0], rename_df[1]))
    return adf, rename_dict


def data_transform(data_obj):
    adf, rename_dict = data_obj
    return adf.rename(rename_dict, axis=1)

def write_data(out_file, data_obj):
    data_obj.to_csv(out_file, sep="\t")


def main():

    in_file, out_file, rename_file = parse_args()

    data = read_data(in_file, rename_file)

    new_data = data_transform(data)

    write_data(out_file, new_data)


if __name__ == "__main__":
    sys.exit(main())
